"""Exact covariance contractions in four absolute variance coordinates.

The blocked information calculation avoids stored n-by-n matrices but still
has quadratic work. It is a validation/reference path, not a full-grid launch.
No inverse information, KR degrees of freedom, or calibrated intervals produced.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from matched_mixed_covariance import MatchedCovariance


COMPONENTS = ('residual', 'background', 'family', 'species')


def covariance_information(background, family, factor, design, variances, block_size=64):
    if isinstance(block_size, bool) or not isinstance(block_size, (int, np.integer)) or block_size < 1:
        raise ValueError('block_size must be a positive integer')
    values = np.asarray(variances, dtype=float)
    if values.shape != (4,):
        raise ValueError('Expected four absolute variances')
    cov = MatchedCovariance(background, family, factor, *values)
    x = np.asarray(design, dtype=float)
    if (x.ndim != 2 or x.shape[0] != cov.n or not 0 < x.shape[1] < cov.n
            or not np.isfinite(x).all() or np.linalg.matrix_rank(x) != x.shape[1]):
        raise ValueError('Expected finite full-rank fixed design and positive residual degrees of freedom')
    t = cov.solve(x)
    phi = cho_solve(cho_factor(x.T @ t, lower=True), np.eye(x.shape[1]))

    def project(a):
        return cov.solve(a) - t @ (phi @ (t.T @ a))

    def kernel(i, a):
        if i == 0:
            return a
        if i in (1, 2):
            labels, groups = (cov.bg, cov.ng) if i == 1 else (cov.fam, cov.nf)
            return cov._group_sum(a, labels, groups)[labels]
        return cov.factor @ (cov.factor.T @ a)

    o = [kernel(i, t) for i in range(4)]
    first = np.array([-t.T @ a for a in o])
    inverse_o = [cov.solve(a) for a in o]
    second = np.array([[a.T @ b for b in inverse_o] for a in o])
    # G_i = U_i U_i'. Sum tr((P U_i)' G_j (P U_i)) in column blocks.
    # This includes the residual identity kernel; zero fitted variances do not
    # remove their derivative kernels or make components automatically estimable.
    information = np.zeros((4, 4))
    sizes = (cov.n, cov.ng, cov.nf, cov.factor.shape[1])
    for i, size in enumerate(sizes):
        for start in range(0, size, block_size):
            stop = min(size, start + block_size)
            if i == 3:
                basis = cov.factor[:, start:stop]
            else:
                labels = np.arange(cov.n) if i == 0 else (cov.bg if i == 1 else cov.fam)
                basis = (labels[:, None] == np.arange(start, stop)[None, :]).astype(float)
            pu = project(basis)
            for j in range(i, 4):
                information[i, j] += .5 * np.sum(pu * kernel(j, pu))
        information[i+1:, i] = information[i, i+1:]
    if not all(np.isfinite(a).all() for a in [phi, first, second, information]):
        raise FloatingPointError('Nonfinite covariance contractions')
    diagonal = np.diag(information)
    if np.any(diagonal < 0):
        raise FloatingPointError('Negative information diagonal')
    scale = np.sqrt(diagonal)
    denominator = np.outer(scale, scale)
    normalized = np.divide(information, denominator, out=np.zeros_like(information), where=denominator > 0)
    eigenvalues = np.linalg.eigvalsh(normalized)
    # Diagnostic numerical rank, not a boundary correction or license to invert.
    tolerance = 4 * np.finfo(float).eps * max(float(np.max(abs(eigenvalues))), 1.)
    return dict(components=COMPONENTS, absolute_variances=values,
                conditional_beta_covariance=phi, first_contractions=first,
                second_contractions=second, expected_reml_information=information,
                normalized_information_eigenvalues=eigenvalues,
                numerical_information_rank=int(np.sum(eigenvalues > tolerance)),
                rank_tolerance=tolerance, exact_zero_components=values == 0)
