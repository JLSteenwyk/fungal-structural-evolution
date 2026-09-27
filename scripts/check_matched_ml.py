#!/usr/bin/env python3
"""Check ML objectives and scores against explicit dense Gaussian calculations."""
import itertools
import json
from pathlib import Path
import numpy as np
from cached_matched_ml import CachedMatchedML
from matched_ml_gradient import evaluate_gradient
from screen_duplication_domain_alignment_coverage import sha


def dense(x, y, kernels, ratios):
    n = len(y)
    v = np.eye(n) + sum(r * k for r, k in zip(ratios, kernels))
    inverse = np.linalg.solve(v, np.eye(n))
    information = x.T @ inverse @ x
    beta = np.linalg.solve(information, x.T @ inverse @ y)
    residual = y - x @ beta
    score_residual = inverse @ residual
    q = float(residual @ score_residual)
    sign, determinant = np.linalg.slogdet(v)
    assert sign == 1
    return dict(beta=beta, profiled_scale=q/n, residual_quadratic=q,
                conditional_beta_covariance=q/n*np.linalg.solve(information, np.eye(x.shape[1])),
                negative_profiled_ml=.5*(determinant+n*(1+np.log(2*np.pi*q/n))),
                gradient=np.array([.5*(np.trace(inverse @ k)-n*(score_residual @ k @ score_residual)/q) for k in kernels]))


def main():
    rng = np.random.default_rng(429671)
    n = 60
    background = np.repeat(np.arange(15), 4)
    family = background // 3
    z = (background[:, None] == np.arange(15)).astype(float)
    g = (family[:, None] == np.arange(5)).astype(float)
    a, b = rng.uniform(size=(2, n))
    x = np.column_stack([np.ones(n), a-b, rng.normal(size=(n, 3)), a*a-b*b, a**3-b**3])
    y = rng.normal(size=n)
    raw = rng.normal(size=(n, 4))
    maximum_objective = maximum_gradient = 0.
    cases = invariance = nesting = 0
    for factor in [np.zeros((n, 0)), raw, np.column_stack([raw, raw[:, 0]])]:
        kernels = [z @ z.T, g @ g.T, factor @ factor.T]
        for ratios in itertools.product([0., .3, 10000.], repeat=3):
            previous = None
            for columns in [5, 6, 7]:
                design = x[:, :columns]
                cache = CachedMatchedML(background, family, factor, design, y)
                actual = cache.evaluate(*ratios)
                expected = dense(design, y, kernels, ratios)
                for key in ['beta', 'profiled_scale', 'residual_quadratic', 'conditional_beta_covariance', 'negative_profiled_ml']:
                    np.testing.assert_allclose(actual[key], expected[key], rtol=2e-7, atol=2e-7)
                assert actual['variance_profile_denominator'] == n
                assert actual['residual_degrees_of_freedom'] == n-columns
                gradient = evaluate_gradient(cache, ratios)
                np.testing.assert_allclose(gradient['negative_profiled_ml'], expected['negative_profiled_ml'], rtol=2e-7, atol=2e-7)
                np.testing.assert_allclose(gradient['ratio_gradient'], expected['gradient'], rtol=2e-7, atol=2e-7)
                np.testing.assert_allclose(gradient['log1p_ratio_gradient'], expected['gradient']*(1+np.array(ratios)), rtol=2e-7, atol=2e-7)
                maximum_objective = max(maximum_objective, abs(actual['negative_profiled_ml']-expected['negative_profiled_ml']))
                maximum_gradient = max(maximum_gradient, float(np.max(abs(gradient['ratio_gradient']-expected['gradient']))))
                # At fixed covariance, adding a column cannot worsen profiled ML.
                if previous is not None:
                    assert actual['negative_profiled_ml'] <= previous + 2e-7
                    nesting += 1
                previous = actual['negative_profiled_ml']
                transform = np.diag(np.linspace(.5, 2., columns))
                transform[0, 1:] = .4
                changed = CachedMatchedML(background, family, factor, design @ transform, y).evaluate(*ratios)
                np.testing.assert_allclose(changed['negative_profiled_ml'], actual['negative_profiled_ml'], atol=2e-7, rtol=2e-7)
                np.testing.assert_allclose(transform @ changed['beta'], actual['beta'], atol=2e-7, rtol=2e-7)
                invariance += 1
                cases += 1
    # Large fitted mean forces direct residual calculation instead of cancellation.
    design = x[:, :5]
    large = design @ np.array([1e7, -3e6, 2e6, 4e6, 1e6]) + rng.normal(size=n)
    actual = CachedMatchedML(background, family, raw, design, large).evaluate(.2, .5, .3)
    assert actual['evaluation_method'] == 'direct_residual_fallback'
    expected = dense(design, large, [z@z.T, g@g.T, raw@raw.T], [.2, .5, .3])
    for key in ['profiled_scale', 'negative_profiled_ml', 'conditional_beta_covariance']:
        np.testing.assert_allclose(actual[key], expected[key], rtol=2e-7, atol=2e-7)
    try:
        CachedMatchedML(background, family, raw, np.column_stack([design, design[:, 0]]), y)
    except ValueError:
        pass
    else:
        raise AssertionError('Rank-deficient design accepted')
    result = dict(status='passed_matched_profiled_ml_dense_checks', cases=cases,
                  ratio_derivatives_checked=3*cases, fixed_design_reparameterizations=invariance,
                  fixed_covariance_nested_comparisons=nesting, cancellation_fallback_cases=1,
                  maximum_objective_error=maximum_objective, maximum_ratio_gradient_error=maximum_gradient,
                  pins={p:sha(p) for p in ['scripts/cached_matched_ml.py', 'scripts/matched_ml_gradient.py', 'scripts/matched_mixed_covariance.py', __file__]},
                  scope='Synthetic full boundary grid, empty/duplicated species factors and nested polynomial designs. Dense Gaussian objective and trace score checks, parameterization invariance and direct-residual fallback. No optimizer, fitted biological effects, calibrated likelihood-ratio null or nonlinear joint-support claim.')
    Path('metadata/matched_ml_checks_20260927.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
