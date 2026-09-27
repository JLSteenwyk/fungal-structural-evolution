#!/usr/bin/env python3
"""Validate ML fitting against dense candidates and a closed-form variance case."""
import json
from pathlib import Path
import numpy as np
from fit_matched_ml import fit_variance_ratios
from check_matched_ml import dense
from screen_duplication_domain_alignment_coverage import sha


def main():
    rng = np.random.default_rng(630951)
    n, groups, size = 60, 15, 4
    bg = np.repeat(np.arange(groups), size)
    family = bg // 3
    z = (bg[:, None] == np.arange(groups)).astype(float)
    g = (family[:, None] == np.arange(5)).astype(float)
    a, b = rng.uniform(size=(2, n))
    x = np.column_stack([np.ones(n), a-b, rng.normal(size=(n, 3)), a*a-b*b, a**3-b**3])
    raw = rng.normal(size=(n, 3))
    rows = []
    candidates_checked = 0
    for factor in [np.zeros((n, 0)), raw]:
        kernels = [z@z.T, g@g.T, factor@factor.T]
        y = rng.normal(size=n) + z@rng.normal(size=groups) + g@rng.normal(size=5)
        if factor.shape[1]:
            y += factor@rng.normal(size=factor.shape[1])
        previous = None
        for columns in [5, 6, 7]:
            result = fit_variance_ratios(bg, family, factor, x[:, :columns], y)
            assert result['species_kernel_is_zero'] == (factor.shape[1] == 0)
            assert not result['component_identifiability_assessed']
            assert len(result['candidates']) == 22
            assert len({tuple(r['active']) for r in result['candidates']}) == 8
            for candidate in result['candidates']:
                expected = dense(x[:, :columns], y, kernels, np.expm1(candidate['theta']))
                np.testing.assert_allclose(candidate['objective'], expected['negative_profiled_ml'], rtol=1e-8, atol=1e-8)
                candidates_checked += 1
            final = dense(x[:, :columns], y, kernels, result['ratios'])
            for key in ['negative_profiled_ml', 'beta', 'profiled_scale', 'conditional_beta_covariance']:
                np.testing.assert_allclose(result[key], final[key], rtol=1e-8, atol=1e-8)
            theta = np.array(result['log1p_ratios'])
            expected_gradient = final['gradient']*(1+np.array(result['ratios']))
            np.testing.assert_allclose(result['analytic_log1p_gradient'], expected_gradient, rtol=1e-7, atol=1e-7)
            projected = np.where(theta <= 1e-7, np.minimum(expected_gradient, 0), expected_gradient)
            projected = np.where(theta >= np.log1p(result['maximum_ratio'])-1e-7, np.maximum(projected, 0), projected)
            np.testing.assert_allclose(result['projected_gradient'], projected, atol=1e-7, rtol=1e-7)
            assert result['status'] == 'ml_candidate_passed_numerical_optimization_checks'
            if previous is not None:
                assert result['negative_profiled_ml'] <= previous + 1e-6
            previous = result['negative_profiled_ml']
            rows.append(dict(factor_columns=factor.shape[1], fixed_columns=columns, status=result['status'], ratios=result['ratios'], objective=previous, wall_seconds=result['wall_seconds']))
    # Balanced intercept-only random-intercept ML: independent closed-form optimum.
    y = 3*z@rng.normal(size=groups) + rng.normal(size=n)
    means = y.reshape(groups, size).mean(axis=1)
    within = float(np.sum((y-np.repeat(means, size))**2))
    between = float(size*np.sum((means-means.mean())**2))
    residual_variance = within/(groups*(size-1))
    between_eigenvalue = between/groups
    ratio = (between_eigenvalue/residual_variance-1)/size
    assert ratio > 0
    result = fit_variance_ratios(bg, np.zeros(n, dtype=int), np.zeros((n, 0)), np.ones((n, 1)), y)
    assert result['status'] == 'ml_candidate_passed_numerical_optimization_checks'
    np.testing.assert_allclose(result['ratios'][0], ratio, rtol=2e-5, atol=1e-7)
    np.testing.assert_allclose(result['ratios'][1], 0, atol=1e-7)
    np.testing.assert_allclose(result['profiled_scale'], residual_variance, rtol=2e-5, atol=1e-7)
    np.testing.assert_allclose(result['beta'][0], y.mean(), rtol=1e-8, atol=1e-8)
    # Exact zero group means put the random-intercept optimum on the zero face.
    centered = y - np.repeat(means, size)
    zero = fit_variance_ratios(bg, np.zeros(n, dtype=int), np.zeros((n, 0)), np.ones((n, 1)), centered)
    np.testing.assert_allclose(zero['ratios'][:2], [0, 0], atol=1e-7)
    np.testing.assert_allclose(zero['profiled_scale'], np.mean(centered**2), rtol=1e-8, atol=1e-8)
    assert zero['status'] == 'ml_candidate_passed_numerical_optimization_checks'
    strong = 100*z@rng.normal(size=groups) + rng.normal(size=n)
    bounded = fit_variance_ratios(bg, np.zeros(n, dtype=int), np.zeros((n, 0)), np.ones((n, 1)), strong, maximum_ratio=21.)
    assert bounded['upper_boundary'][0]
    assert bounded['status'] == 'ml_candidate_requires_optimization_review'
    # One iteration must retain failed attempts and flag the candidate for review.
    limited = fit_variance_ratios(bg, family, raw, x[:, :5], y, maxiter=1)
    assert len(limited['candidates']) == 22 and any(not r['success'] for r in limited['candidates'])
    assert limited['status'] == 'ml_candidate_requires_optimization_review'
    paths = ['scripts/fit_matched_ml.py', 'scripts/cached_matched_ml.py', 'scripts/matched_ml_gradient.py', 'scripts/matched_mixed_covariance.py', 'scripts/check_matched_ml.py', __file__]
    proof = dict(status='passed_matched_ml_optimizer_numerical_checks', synthetic_fits=rows,
                 dense_candidate_checks=candidates_checked, closed_form_random_intercept_checks=1,
                 zero_random_intercept_boundary_checks=1, upper_bound_review_checks=1,
                 iteration_limit_failure_retention_checks=1, pins={p:sha(p) for p in paths},
                 scope='Numerical optimizer checks only, not global optimality, calibrated uncertainty, model adequacy or biological evidence. All 22 candidates retained per fit; analytic score cancellation fails explicitly. No production ML fit launched.')
    Path('metadata/matched_ml_optimizer_checks_20260927.json').write_text(json.dumps(proof, indent=2)+'\n')
    print(json.dumps(proof, indent=2))


if __name__ == '__main__':
    main()
