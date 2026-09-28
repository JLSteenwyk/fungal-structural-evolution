"""Rank-one KR moment calculation; candidate intervals, not calibrated inference."""
import numpy as np
from scipy.stats import t
from matched_kr_covariance import adjust_covariance


def scalar_interval(contractions, beta, contrast, alpha=.05):
    phi = np.asarray(contractions['conditional_beta_covariance'], dtype=float)
    beta, contrast = np.asarray(beta, float), np.asarray(contrast, float)
    if (beta.shape != (len(phi),) or contrast.shape != beta.shape
            or not np.isfinite(beta).all() or not np.isfinite(contrast).all()
            or not np.any(contrast) or not np.isfinite(alpha) or not 0 < alpha < 1):
        raise ValueError('Expected finite coefficient/contrast vectors, nonzero contrast and 0<alpha<1')
    adjusted = adjust_covariance(contractions)
    if adjusted['status'] != 'candidate_adjustment_pending_statistical_validation':
        return adjusted
    s = float(contrast @ phi @ contrast)
    h = phi @ contrast
    # Rank-one Theta implies A1=A2, hence denominator df=2/A2 and F scale=1.
    derivatives = np.array([h @ a @ h for a in contractions['first_contractions']])
    a2 = float(derivatives @ adjusted['variance_parameter_covariance'] @ derivatives / s**2)
    if not np.isfinite(a2) or a2 <= 0:
        return dict(status='scalar_moments_require_review', reason='nonpositive_or_nonfinite_A2')
    df = 2/a2
    variance = float(contrast @ adjusted['adjusted_beta_covariance'] @ contrast)
    estimate = float(contrast @ beta)
    half = float(t.isf(alpha/2, df)*np.sqrt(variance))
    if not np.isfinite([df, variance, half]).all() or variance <= 0:
        return dict(status='scalar_interval_requires_review', reason='nonfinite_or_nonpositive_result')
    return dict(status='candidate_interval_pending_coverage_validation', estimate=estimate,
                adjusted_variance=variance, denominator_df=df, numerator_df=1,
                f_scaling=1., A2=a2, nominal_level=1-alpha,
                lower=estimate-half, upper=estimate+half,
                exact_zero_components=adjusted['exact_zero_components'].tolist())
