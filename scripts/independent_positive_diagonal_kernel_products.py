"""Independent latent-space raw/error-contrast products for a positive diagonal.

Unlike the frozen uniform-only reader, retain every D contraction. This is
numerical audit algebra, not basis selection, covariance fitting or precision
calibration. Entity kernels use latent overlaps, without dense component kernels.
"""
import numpy as np
from scipy import sparse, linalg


class PositiveDiagonalKernelProducts:
    def __init__(self, labels, incidence, diagonal):
        labels = np.asarray(labels)
        if labels.ndim != 1 or not len(labels): raise ValueError('Nonempty component labels required')
        self.n = len(labels); self.diagonal = np.asarray(diagonal, dtype=float)
        if self.diagonal.shape != (self.n,) or not np.isfinite(self.diagonal).all() or np.any(self.diagonal <= 0):
            raise ValueError('Strictly positive finite residual diagonal required')
        if any(not isinstance(k, str) or not k for k in incidence) or {'residual','species'} & set(incidence):
            raise ValueError('Distinct nonempty nonreserved incidence names required')
        _, codes = np.unique(labels, return_inverse=True)
        self.components = len(set(codes)); self.names = ['residual', *incidence, 'species']
        self.operators = []; self.widths = []
        for value in incidence.values():
            z = sparse.csr_matrix(value, dtype=float, copy=True)
            if z.shape[0] != self.n or not np.isfinite(z.data).all(): raise ValueError('Invalid entity operator')
            z.sum_duplicates(); z.eliminate_zeros(); z.sort_indices()
            z = z[:, np.unique(z.indices)].tocsr()
            # An entity is required to belong to one declared component.
            present = {}
            rr, cc = z.nonzero()
            for row, col in zip(rr, cc):
                if int(col) in present and present[int(col)] != int(codes[row]):
                    raise ValueError('Shared entity crosses declared components')
                present[int(col)] = int(codes[row])
            self.operators.append(z); self.widths.append(z.shape[1])
        k = len(self.names); self.base = np.zeros((k,k)); self.cross = {}
        self.base[0,0] = float(self.diagonal @ self.diagonal)
        for i,a in enumerate(self.operators):
            energy = np.asarray(a.multiply(a).sum(axis=1)).ravel()
            self.base[0,i+1] = float(self.diagonal @ energy)
            for j,b in enumerate(self.operators):
                cross = (a.T @ b).tocsr(); self.cross[i,j] = cross
                self.base[i+1,j+1] = float(cross.multiply(cross).sum())
        self.base[1:,0] = self.base[0,1:]

    def tree(self, factor):
        return PositiveDiagonalTreeProducts(self, factor)


class PositiveDiagonalTreeProducts:
    def __init__(self, bank, factor):
        self.bank = bank; self.factor = np.asarray(factor, dtype=float)
        if self.factor.ndim != 2 or self.factor.shape[0] != bank.n or not np.isfinite(self.factor).all():
            raise ValueError('Finite n-row species factor required')
        self.raw = bank.base.copy(); self.species_core = self.factor.T @ self.factor
        self.raw[-1,-1] = float(np.sum(self.species_core * self.species_core))
        self.raw[0,-1] = self.raw[-1,0] = float(bank.diagonal @ np.sum(self.factor * self.factor, axis=1))
        for i,z in enumerate(bank.operators,1):
            latent = z.T @ self.factor; self.raw[i,-1] = float(np.sum(latent * latent))
        self.raw[-1,1:-1] = self.raw[1:-1,-1]

    def project(self, design):
        bank = self.bank; n = bank.n; x = np.asarray(design, dtype=float)
        if x.ndim != 2 or x.shape[0] != n or not 0 < x.shape[1] < n or not np.isfinite(x).all():
            raise ValueError('Finite active design with residual dimensions required')
        scaled = np.asarray(x, dtype=np.longdouble)
        maximum = np.max(abs(scaled), axis=0)
        if np.any(maximum == 0): raise ValueError('Inactive zero design column')
        scaled /= maximum; scaled /= np.sqrt(np.sum(scaled * scaled, axis=0))
        scaled = np.asarray(scaled, dtype=float)
        singular = linalg.svd(scaled, compute_uv=False, lapack_driver='gesvd')
        if singular[-1] <= 10 * max(scaled.shape) * np.finfo(float).eps * singular[0]:
            raise ValueError('Design rank deficiency or boundary requires review')
        q,_,_ = linalg.qr(scaled, mode='economic', pivoting=True)
        k = len(bank.names); p = x.shape[1]
        core = np.zeros((k,p,p)); correction = np.zeros((k,k))
        dq = bank.diagonal[:,None] * q
        core[0] = q.T @ dq; correction[0,0] = float(np.sum(dq * dq))
        species = q.T @ self.factor
        core[-1] = species @ species.T
        correction[0,-1] = correction[-1,0] = float(np.sum((dq.T @ self.factor) * species))
        coefficients = [np.asarray(z.T @ q).T for z in bank.operators]
        for i,a in enumerate(coefficients,1):
            core[i] = a @ a.T
            correction[0,i] = correction[i,0] = float(np.sum(np.asarray(bank.operators[i-1].T @ dq).T * a))
        for i,a in enumerate(coefficients,1):
            for j,b in enumerate(coefficients,1):
                correction[i,j] = float(np.sum((bank.cross[i-1,j-1].T @ a.T).T * b))
            image = bank.operators[i-1] @ a.T
            correction[i,-1] = float(np.sum((image.T @ self.factor) * species))
        correction[-1,1:-1] = correction[1:-1,-1]
        correction[-1,-1] = float(np.sum(self.species_core * (species.T @ species)))
        small = np.asarray([[float(np.sum(a*b)) for b in core] for a in core])
        projected = self.raw - 2 * correction + small
        dimension = max(n, self.factor.shape[1], p, *bank.widths, 1)
        envelope = 64 * np.finfo(float).eps * dimension
        condition = float(singular[0] / singular[-1])
        return dict(raw=self.raw.copy(), projected=projected,
            raw_error=envelope * abs(self.raw),
            projected_error=envelope * condition * (abs(self.raw) + 2*abs(correction) + abs(small)),
            normalized_design_condition_number=condition,
            diagonal_core=core[0].copy(), diagonal_image_energy=correction[0,0],
            scope='Fresh positive-diagonal latent-overlap raw/REML products. No uniform substitution, inherited numerical envelope, covariance selection, clipping or variance fit.')
