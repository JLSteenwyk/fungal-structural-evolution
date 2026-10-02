"""Audit covariance bases before REML; do not select or fit a variance model.

For H=I-QQ', REML observes H K H rather than the unprojected kernel K.
All Frobenius products use component-local entity kernels and low-rank
projections. No full n-by-n covariance or residual projector is allocated.
Cancellation and rank boundaries are review states, never repaired by clipping.
"""
import numpy as np
from scipy import sparse
from scipy.linalg import qr


def _diagnostics(gram, error_bound, names):
    diagonal = np.diag(gram)
    if not np.isfinite(gram).all() or np.any(diagonal < -np.diag(error_bound)):
        raise ArithmeticError('Kernel Gram lost finite nonnegative squared norms')
    resolved = diagonal > np.diag(error_bound)
    active = np.flatnonzero(resolved)
    unresolved = np.flatnonzero(~resolved)
    if len(active):
        normalization = np.sqrt(np.outer(diagonal[active], diagonal[active]))
        normalized = gram[np.ix_(active, active)] / normalization
        normalized_error = error_bound[np.ix_(active, active)] / normalization
        singular = np.linalg.svd(normalized, compute_uv=False)
        roundoff = float(np.linalg.norm(normalized_error, ord=2))
        if np.linalg.eigvalsh((normalized + normalized.T) / 2)[0] < -roundoff:
            raise ArithmeticError('Kernel Gram lost positive semidefiniteness beyond its review bound')
        tolerance = max(len(active) * np.finfo(float).eps * singular[0], roundoff)
        ranks = [int(np.count_nonzero(singular > tolerance * multiplier))
                 for multiplier in (.1, 1., 10.)]
    else:
        singular = np.empty(0); roundoff = tolerance = 0.; ranks = [0, 0, 0]
    status = ('unresolved_kernel_norm_requires_review' if len(unresolved) else
              'rank_boundary_requires_review' if ranks[0] != ranks[2] else
              'dependent_covariance_bases_require_review' if ranks[1] < len(active) else
              'numerically_independent_covariance_bases')
    return dict(disposition=status, resolved_kernel_indices=active.tolist(),
                unresolved_kernel_names=[names[i] for i in unresolved],
                normalized_singular_values=singular.tolist(),
                propagated_normalized_roundoff_bound=roundoff, rank_tolerance=float(tolerance),
                rank_at_one_tenth_tolerance=ranks[0], rank=ranks[1],
                rank_at_ten_times_tolerance=ranks[2])


def audit_covariance_bases(block_labels, incidence, species_factor,
                           residual_diagonal, design):
    """Return raw/projected Grams, numerical ranks and explicit review states.

    Every incidence column must be contained in one declared component.
    The residual diagonal is a supplied positive covariance basis; it need
    not be uniform. Names, zero bases and dependent bases remain in the audit.
    No returned rank authorizes arbitrary basis deletion: that could change
    the nonnegative variance cone and the scientific model.
    """
    labels = np.asarray(block_labels)
    if labels.ndim != 1 or not len(labels):
        raise ValueError('Nonempty one-dimensional block labels required')
    n = len(labels); _, codes = np.unique(labels, return_inverse=True)
    diagonal = np.asarray(residual_diagonal, dtype=float)
    factor = np.asarray(species_factor, dtype=float)
    x = np.asarray(design, dtype=float)
    if diagonal.shape != (n,) or not np.isfinite(diagonal).all() or np.any(diagonal <= 0):
        raise ValueError('Positive finite residual basis required')
    if factor.ndim != 2 or factor.shape[0] != n or not np.isfinite(factor).all():
        raise ValueError('Finite species factor with n rows required')
    if x.ndim != 2 or x.shape[0] != n or not 0 < x.shape[1] < n or not np.isfinite(x).all():
        raise ValueError('Finite fixed design with positive residual dimension required')
    if any(not isinstance(name, str) or not name for name in incidence) or {'residual', 'species'} & set(incidence):
        raise ValueError('Distinct nonempty entity names must not use reserved basis names')
    maximum = np.max(abs(x), axis=0)
    if np.any(maximum == 0):
        raise ValueError('Audit requires the active design; retain zero terms separately')
    scaled = x / maximum
    scaled /= np.sqrt(np.sum(scaled * scaled, axis=0))
    singular = np.linalg.svd(scaled, compute_uv=False)
    cutoff = max(x.shape) * np.finfo(float).eps * singular[0]
    if singular[-1] <= 10 * cutoff:
        raise ValueError('Fixed design is rank deficient or at its rank review boundary')
    q, _ = qr(scaled, mode='economic')
    operators = {}
    for name, original in incidence.items():
        z = sparse.csr_matrix(original, dtype=float, copy=True)
        if z.shape[0] != n or not np.isfinite(z.data).all():
            raise ValueError('Invalid incidence operator')
        z.sum_duplicates(); z.eliminate_zeros(); z.sort_indices()
        rows, cols = z.nonzero()
        low = np.full(z.shape[1], n, dtype=np.int64)
        high = np.full(z.shape[1], -1, dtype=np.int64)
        np.minimum.at(low, cols, codes[rows]); np.maximum.at(high, cols, codes[rows])
        present = high >= 0
        if np.any(low[present] != high[present]):
            raise ValueError('Shared entity crosses declared components')
        operators[name] = z
    names = ['residual', *operators, 'species']; k = len(names)
    raw = np.zeros((k, k)); raw[0, 0] = diagonal @ diagonal
    order = np.argsort(codes, kind='stable')
    parts = np.split(order, np.flatnonzero(np.diff(codes[order])) + 1)
    for part in parts:
        local = []
        for i, z in enumerate(operators.values(), 1):
            zs = z[part]
            # Empty global columns must never become a dense latent matrix.
            zs = zs[:, np.unique(zs.indices)]
            kernel = (zs @ zs.T).toarray(); local.append(kernel)
            raw[0, i] += diagonal[part] @ np.diag(kernel)
            cross = zs.T @ factor[part]
            raw[i, -1] += np.sum(cross * cross)
        if local:
            flat = np.asarray(local).reshape(len(local), -1)
            raw[1:-1, 1:-1] += flat @ flat.T
    raw[-1, -1] = np.sum((factor.T @ factor)**2)
    raw[0, -1] = diagonal @ np.sum(factor * factor, axis=1)
    raw[1:, 0] = raw[0, 1:]; raw[-1, 1:-1] = raw[1:-1, -1]
    images = [diagonal[:, None] * q]
    images.extend(z @ (z.T @ q) for z in operators.values())
    images.append(factor @ (factor.T @ q))
    images = np.asarray(images)
    image_products = images.reshape(k, -1) @ images.reshape(k, -1).T
    cores = np.asarray([q.T @ value for value in images]).reshape(k, -1)
    core_products = cores @ cores.T
    projected = raw - 2 * image_products + core_products
    # A conservative rounding envelope for the subtraction and contractions.
    # It is a numerical review bound, not statistical uncertainty or clipping.
    condition = float(singular[0] / singular[-1])
    factor_bound = 64 * np.finfo(float).eps * max(
        n, factor.shape[1], x.shape[1], *(z.shape[1] for z in operators.values()), 1)
    raw_error = factor_bound * np.abs(raw)
    projected_error = factor_bound * condition * (abs(raw) + 2 * abs(image_products) + abs(core_products))
    return dict(kernel_names=names, records=n, fixed_effect_columns=x.shape[1],
                residual_dimension=n-x.shape[1], family_components=len(parts),
                normalized_design_condition_number=condition,
                raw_gram=raw.tolist(), projected_gram=projected.tolist(),
                raw_roundoff_envelope=raw_error.tolist(),
                projected_roundoff_envelope=projected_error.tolist(),
                exactly_zero_incidence_names=[name for name, z in operators.items() if not z.nnz],
                raw_diagnostics=_diagnostics(raw, raw_error, names),
                reml_diagnostics=_diagnostics(projected, projected_error, names),
                covariance_model_selected=False,
                scope='Numerical structural identifiability audit only. No variance fitting, '
                      'basis removal, clipping, scientific adequacy or calibration. '
                      'Gram conditioning squares basis conditioning; unresolved norms '
                      'and rank boundaries require review.')
