#!/usr/bin/env python3
"""Whole-protein adapter for the validated ordinary-ML covariance optimizer."""
import numpy as np
from fit_matched_ml import fit_variance_ratios
from matched_mixed_covariance import MatchedCovariance
from cached_matched_ml import profiled_ml


def fit_input(background, family, factor, matrix, columns, maximum_ratio=10000., maxiter=300):
    values = np.asarray(matrix, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(columns) or len(columns) < 2:
        raise ValueError('Expected outcome, intercept and optional active covariates')
    if columns[1] != 'intercept' or len(set(columns)) != len(columns):
        raise ValueError('Invalid column identities')
    if not np.isfinite(values).all() or not np.all(values[:, 1] == 1):
        raise ValueError('Nonfinite values or invalid intercept')
    if len(values) <= values.shape[1] - 1:
        raise ValueError('No positive residual degrees of freedom')
    covariates = values[:, 2:]
    if np.any(np.ptp(covariates, axis=0) <= 1e-12):
        raise ValueError('Inventory must remove inactive covariates before fitting')
    scales = np.std(covariates, axis=0)
    design = np.column_stack([values[:, 1], covariates / scales])
    response = values[:, 0]
    fitted = fit_variance_ratios(background, family, factor, design, response,
                                 maximum_ratio=maximum_ratio, maxiter=maxiter)
    maximum_error = 0.
    evaluated = {}
    for candidate in fitted['candidates']:
        key = tuple(candidate['theta'])
        if key not in evaluated:
            evaluated[key] = profiled_ml(MatchedCovariance(background, family, factor, 1., *np.expm1(key)), design, response)
        direct = evaluated[key]['negative_profiled_ml']
        np.testing.assert_allclose(direct, candidate['objective'], rtol=1e-9, atol=1e-7)
        maximum_error = max(maximum_error, abs(direct - candidate['objective']))
    direct = profiled_ml(MatchedCovariance(background, family, factor, 1., *fitted['ratios']), design, response)
    for key in ['beta', 'profiled_scale', 'residual_quadratic', 'conditional_beta_covariance']:
        np.testing.assert_allclose(fitted[key], direct[key], rtol=1e-7, atol=1e-8)
    conversion = np.r_[1., 1. / scales]
    fitted.update(likelihood='ordinary_gaussian_ml', records=len(values), columns=list(columns),
                  covariate_scales=scales.tolist(),
                  raw_unit_beta=(np.asarray(fitted['beta']) * conversion).tolist(),
                  raw_unit_conditional_beta_covariance=(np.asarray(fitted['conditional_beta_covariance']) * np.outer(conversion, conversion)).tolist(),
                  direct_candidate_readback=dict(candidates_checked=len(fitted['candidates']), distinct_parameter_vectors=len(evaluated), maximum_objective_error=maximum_error),
                  scope='Whole-protein conditional working-model candidate. All optimizer starts and review flags retained; no calibrated interval, model-comparison decision, causal duplication claim or biological acceleration inference.')
    return fitted
