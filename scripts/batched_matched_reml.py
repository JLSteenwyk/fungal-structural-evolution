"""Shared-covariance REML evaluations for multiple responses, without fitting ratios."""
import numpy as np
from scipy.linalg import cho_factor, cho_solve


def evaluate(covariance, design, responses):
    x, y = np.asarray(design, dtype=float), np.asarray(responses, dtype=float)
    if x.ndim != 2 or y.ndim != 2 or x.shape[0] != covariance.n or y.shape[0] != covariance.n or not y.shape[1]:
        raise ValueError('Expected matching observation-by-column matrices')
    n, p = x.shape
    if not p or n <= p or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Finite data and positive residual degrees of freedom required')
    if np.linalg.matrix_rank(x) != p:
        raise ValueError('Rank-deficient design')
    information = x.T @ covariance.solve(x)
    information = (information+information.T)/2
    chol = cho_factor(information, lower=True)
    beta = cho_solve(chol, x.T @ covariance.solve(y))
    residual = y-x@beta
    quadratic = np.sum(residual*covariance.solve(residual), axis=0)
    valid = np.isfinite(quadratic) & (quadratic > 0)
    scale = np.full(y.shape[1], np.nan)
    scale[valid] = quadratic[valid]/(n-p)
    objective = np.full(y.shape[1], np.nan)
    determinant = 2*np.log(np.diag(chol[0])).sum()
    objective[valid] = .5*(covariance.logdet+determinant+(n-p)*(1+np.log(2*np.pi*scale[valid])))
    return dict(valid=valid, beta=beta, profiled_scale=scale,
                negative_profiled_reml=objective, residual_quadratic=quadratic,
                conditional_beta_covariance=scale[:,None,None]*cho_solve(chol,np.eye(p))[None,:,:],
                residual_degrees_of_freedom=n-p)
