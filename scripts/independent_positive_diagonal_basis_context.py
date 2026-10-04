"""Reuse independently derived latent products, recomputing each D correction.

Long-double normalization/pivoted QR/gesvd and latent overlaps deliberately
differ from the component/image-product producer. No rank-based basis changes,
uniform substitutions or inherited numerical audits/envelopes are used.
"""
import numpy as np
from scipy import linalg

from independent_positive_diagonal_kernel_products import PositiveDiagonalKernelProducts


class IndependentPositiveDiagonalBasisContext:
    def __init__(self, labels, incidence, factor):
        self.bank = PositiveDiagonalKernelProducts(labels, incidence, np.ones(len(labels)))
        self.factor = np.array(factor, dtype=float, copy=True)
        self.tree = self.bank.tree(self.factor)
        self.factor.flags.writeable = False
        self.diagonal_energies = np.asarray([
            *[np.asarray(z.multiply(z).sum(axis=1)).ravel() for z in self.bank.operators],
            np.sum(self.factor * self.factor, axis=1)])
        self.diagonal_energies.flags.writeable = False

    def design(self, design):
        return IndependentPositiveDiagonalDesignContext(self, design)


class IndependentPositiveDiagonalDesignContext:
    def __init__(self, parent, design):
        self.parent = parent; bank = parent.bank; x = np.asarray(design, dtype=float)
        n = bank.n
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
        self.q = linalg.qr(scaled, mode='economic', pivoting=True)[0]
        self.q.flags.writeable = False
        self.p = x.shape[1]; self.condition = float(singular[0] / singular[-1])
        self.envelope = 64 * np.finfo(float).eps * max(n, parent.factor.shape[1], self.p, *bank.widths, 1)
        self.species = self.q.T @ parent.factor
        self.coefficients = [np.asarray(z.T @ self.q).T for z in bank.operators]
        k = len(bank.names); core = np.zeros((k, self.p, self.p)); correction = np.zeros((k, k))
        core[-1] = self.species @ self.species.T
        for i, a in enumerate(self.coefficients, 1):
            core[i] = a @ a.T
            for j, b in enumerate(self.coefficients, 1):
                correction[i, j] = float(np.sum((bank.cross[i - 1, j - 1].T @ a.T).T * b))
            image = bank.operators[i - 1] @ a.T
            correction[i, -1] = float(np.sum((image.T @ parent.factor) * self.species))
        correction[-1, 1:-1] = correction[1:-1, -1]
        correction[-1, -1] = float(np.sum(parent.tree.species_core * (self.species.T @ self.species)))
        self.core = core; self.correction = correction
        self.small = np.asarray([[float(np.sum(a * b)) for b in core] for a in core])
        for a in [self.core, self.correction, self.small, self.species, *self.coefficients]:
            a.flags.writeable = False

    def audit(self, diagonal):
        bank = self.parent.bank; d = np.asarray(diagonal, dtype=float)
        if d.shape != (bank.n,) or not np.isfinite(d).all() or np.any(d <= 0):
            raise ValueError('Strictly positive finite residual diagonal required')
        raw = self.parent.tree.raw.copy()
        raw[0, 0] = float(d @ d)
        raw[0, 1:] = self.parent.diagonal_energies @ d; raw[1:, 0] = raw[0, 1:]
        dq = d[:, None] * self.q; dcore = self.q.T @ dq
        correction = self.correction.copy(); small = self.small.copy()
        correction[0, 0] = float(np.sum(dq * dq))
        for i, a in enumerate(self.coefficients, 1):
            correction[0, i] = float(np.sum(np.asarray(bank.operators[i - 1].T @ dq).T * a))
        correction[0, -1] = float(np.sum((dq.T @ self.parent.factor) * self.species))
        correction[1:, 0] = correction[0, 1:]
        small[0, 0] = float(np.sum(dcore * dcore))
        for j in range(1, len(bank.names)): small[0, j] = float(np.sum(dcore * self.core[j]))
        small[1:, 0] = small[0, 1:]
        projected = raw - 2 * correction + small
        return dict(raw=raw, projected=projected, raw_error=self.envelope * abs(raw),
            projected_error=self.envelope * self.condition * (abs(raw) + 2 * abs(correction) + abs(small)),
            normalized_design_condition_number=self.condition, diagonal_core=dcore.copy(),
            diagonal_image_energy=correction[0, 0],
            scope='Fresh general-D latent error-contrast products; unchanged nonresidual overlaps shared across supplied diagonals, every residual contraction and bound rebuilt. No numerical audit inherited, basis deletion or covariance fit.')
