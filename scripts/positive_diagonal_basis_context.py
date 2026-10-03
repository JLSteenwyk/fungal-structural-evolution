"""Reuse fresh nonresidual products across diagonals, without inherited audits.

All kernels, rank diagnostics and error bounds are newly computed. No saved
uniform qualification or error envelope is consumed. Retaining/removing any
basis requires a separate named source certificate; this module never does it.
"""
import numpy as np
from scipy.linalg import qr

from covariance_basis_audit import _diagnostics
from covariance_basis_context import ComponentKernelProducts


class PositiveDiagonalBasisContext:
    def __init__(self, labels, incidence, factor):
        # Ones are a disposable placeholder for the residual row only. All
        # nonresidual products are unchanged by the eventual supplied D.
        self.bank = ComponentKernelProducts(labels, incidence, np.ones(len(labels)))
        self.factor = np.array(factor, dtype=float, copy=True)
        self.factor.flags.writeable = False
        self.tree = self.bank.tree(self.factor)
        self.names = list(self.bank.names)
        self.diagonal_images = np.column_stack([
            *[np.asarray(z.multiply(z).sum(axis=1)).ravel() for z in self.bank.operators.values()],
            np.sum(self.factor * self.factor, axis=1)])
        self.diagonal_images.flags.writeable = False

    def design(self, design):
        return PositiveDiagonalDesignContext(self, design)


class PositiveDiagonalDesignContext:
    def __init__(self, parent, design):
        self.parent = parent; bank = parent.bank
        x = np.asarray(design, dtype=float)
        if x.ndim != 2 or x.shape[0] != bank.n or not 0 < x.shape[1] < bank.n or not np.isfinite(x).all():
            raise ValueError('Finite active design with residual dimensions required')
        maximum = np.max(abs(x), axis=0)
        if np.any(maximum == 0): raise ValueError('Inactive zero design columns require separate accounting')
        scaled = x / maximum; scaled /= np.sqrt(np.sum(scaled * scaled, axis=0))
        singular = np.linalg.svd(scaled, compute_uv=False)
        if singular[-1] <= 10 * max(x.shape) * np.finfo(float).eps * singular[0]:
            raise ValueError('Design rank or boundary requires review')
        self.q = qr(scaled, mode='economic')[0]; self.q.flags.writeable = False
        self.columns = x.shape[1]; self.condition = float(singular[0] / singular[-1])
        images = [z @ (z.T @ self.q) for z in bank.operators.values()]
        images.append(parent.factor @ (parent.factor.T @ self.q))
        self.images = np.asarray(images); self.images.flags.writeable = False
        flat = self.images.reshape(len(images), -1)
        self.products = flat @ flat.T
        self.cores = np.asarray([self.q.T @ value for value in images])
        self.cores.flags.writeable = False
        small = self.cores.reshape(len(images), -1)
        self.core_products = small @ small.T
        dimension = max(bank.n, parent.factor.shape[1], self.columns,
                        *(z.shape[1] for z in bank.operators.values()), 1)
        self.envelope = 64 * np.finfo(float).eps * dimension

    def audit(self, diagonal):
        parent = self.parent; bank = parent.bank
        d = np.asarray(diagonal, dtype=float)
        if d.shape != (bank.n,) or not np.isfinite(d).all() or np.any(d <= 0):
            raise ValueError('Strictly positive finite residual diagonal required')
        raw = parent.tree.raw.copy()
        raw[0,0] = d @ d
        raw[0,1:] = d @ parent.diagonal_images; raw[1:,0] = raw[0,1:]
        k = len(parent.names); pair = np.empty((k,k)); small = np.empty((k,k))
        pair[1:,1:] = self.products; small[1:,1:] = self.core_products
        dq = d[:,None] * self.q; core = self.q.T @ dq
        pair[0,0] = np.sum(dq * dq)
        pair[0,1:] = dq.ravel() @ self.images.reshape(k-1,-1).T
        pair[1:,0] = pair[0,1:]
        small[0,0] = np.sum(core * core)
        small[0,1:] = core.ravel() @ self.cores.reshape(k-1,-1).T
        small[1:,0] = small[0,1:]
        projected = raw - 2 * pair + small
        raw_error = self.envelope * abs(raw)
        projected_error = self.envelope * self.condition * (abs(raw) + 2*abs(pair) + abs(small))
        return dict(kernel_names=list(parent.names), records=bank.n, fixed_effect_columns=self.columns,
            residual_dimension=bank.n-self.columns, family_components=len(bank.parts),
            normalized_design_condition_number=self.condition,
            raw_gram=raw.tolist(), projected_gram=projected.tolist(),
            raw_roundoff_envelope=raw_error.tolist(), projected_roundoff_envelope=projected_error.tolist(),
            exactly_zero_incidence_names=[name for name,z in bank.operators.items() if not z.nnz],
            raw_diagnostics=_diagnostics(raw,raw_error,parent.names),
            reml_diagnostics=_diagnostics(projected,projected_error,parent.names),
            covariance_model_selected=False,
            scope='Fresh component and tree products shared only across supplied residual diagonals; fixed-design images/cores reused without changing any kernel. Each D Gram and error envelope is recomputed. No saved uniform audit/envelope inherited, automatic basis deletion, clipping, fit or inferential calibration.')
