"""Candidate KR covariance adjustment for linear covariance kernels.

Numerical formula only: no reference distribution, interval or coverage claim.
Unlike pbkrtest's fallback, singular information is explicitly unqualified.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve


def adjust_covariance(contractions):
    phi = np.asarray(contractions['conditional_beta_covariance'], dtype=float)
    first = np.asarray(contractions['first_contractions'], dtype=float)
    second = np.asarray(contractions['second_contractions'], dtype=float)
    information = np.asarray(contractions['expected_reml_information'], dtype=float)
    if phi.ndim != 2 or phi.shape[0] != phi.shape[1]:
        raise ValueError('Expected square fixed-effect covariance')
    p = phi.shape[0]
    if first.shape != (4, p, p) or second.shape != (4, 4, p, p) or information.shape != (4, 4):
        raise ValueError('Expected four absolute-variance kernel contractions')
    if not all(np.isfinite(a).all() for a in [phi, first, second, information]):
        raise ValueError('Nonfinite contractions')
    np.testing.assert_allclose(information, information.T, rtol=1e-10, atol=1e-12)
    scale = np.sqrt(np.maximum(np.diag(information), 0))
    if np.any(scale == 0):
        return dict(status='information_requires_review', reason='nonpositive_information_diagonal')
    normalized = information / np.outer(scale, scale)
    eigenvalues = np.linalg.eigvalsh(normalized)
    # This numerical eligibility threshold is not a statistical regularity test.
    if eigenvalues[0] <= 1e-10 * eigenvalues[-1]:
        return dict(status='information_requires_review', reason='singular_or_ill_conditioned_information')
    w = cho_solve(cho_factor(normalized, lower=True), np.eye(4)) / np.outer(scale, scale)
    u = sum(w[i, j] * (second[i, j] - first[i] @ phi @ first[j])
            for i in range(4) for j in range(4))
    adjusted = phi + 2 * phi @ u @ phi
    adjusted = (adjusted + adjusted.T) / 2
    if not np.isfinite(adjusted).all() or np.linalg.eigvalsh(adjusted)[0] <= 0:
        return dict(status='adjusted_covariance_requires_review', reason='nonpositive_or_nonfinite_covariance')
    return dict(status='candidate_adjustment_pending_statistical_validation',
                adjusted_beta_covariance=adjusted, variance_parameter_covariance=w,
                exact_zero_components=np.asarray(contractions['exact_zero_components'], dtype=bool),
                normalized_information_eigenvalues=eigenvalues)
