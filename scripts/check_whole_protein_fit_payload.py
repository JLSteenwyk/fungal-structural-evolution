#!/usr/bin/env python3
"""Recompute saved whole-protein likelihood candidates and numerical decisions."""
import itertools
import numpy as np
from cached_matched_ml import CachedMatchedML
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from readback_selected_refinement_candidates import independent_fit


def check_payload(payload, background, family, factor, matrix, columns):
    assert payload['status'] in ['ml_candidate_passed_numerical_optimization_checks', 'ml_candidate_requires_optimization_review']
    values = np.asarray(matrix)
    assert np.isfinite(values).all() and np.all(values[:,1] == 1)
    assert payload['columns'] == columns and payload['likelihood'] == 'ordinary_gaussian_ml'
    assert payload['records'] == payload['variance_profile_denominator'] == len(values)
    assert payload['residual_degrees_of_freedom'] == len(values)-len(columns)+1
    scales = np.std(values[:,2:],axis=0)
    assert np.all(scales > 0)
    np.testing.assert_array_equal(payload['covariate_scales'],scales)
    design = np.column_stack([np.ones(len(values)),values[:,2:]/scales])
    response = values[:,0]
    upper = np.log1p(payload['maximum_ratio'])
    assert payload['maximum_ratio'] == 10000.
    required = {(face,start) for face in itertools.product([False,True],repeat=3) for start in ([.05,1.,20.] if any(face) else [0.])}
    observed = set()
    cache = {}
    maximum_error = 0.
    for c in payload['candidates']:
        face = tuple(c['active'])
        key = face,c['start_ratio']
        assert key in required and key not in observed
        observed.add(key)
        theta = np.array(c['theta'])
        assert theta.shape == (3,) and np.isfinite(theta).all()
        assert np.all(theta >= 0) and np.all(theta <= upper)
        assert np.all(theta[~np.asarray(face)] == 0)
        point = tuple(theta)
        if point not in cache:
            cache[point] = independent_fit(MatchedCovariance(background,family,factor,1.,*np.expm1(theta)),design,response)
        reference = cache[point]['negative_profiled_ml']
        np.testing.assert_allclose(c['objective'],reference,rtol=1e-9,atol=1e-7)
        maximum_error = max(maximum_error,abs(c['objective']-reference))
    assert observed == required and len(payload['candidates']) == 22
    best = min(payload['candidates'],key=lambda c:c['objective'])
    theta = np.asarray(best['theta'])
    np.testing.assert_array_equal(payload['log1p_ratios'],theta)
    np.testing.assert_array_equal(payload['ratios'],np.expm1(theta))
    fitted = cache[tuple(theta)]
    for key,value in fitted.items():
        np.testing.assert_allclose(payload[key],value,rtol=1e-7,atol=1e-8)
    conversion = np.r_[1.,1/scales]
    np.testing.assert_allclose(payload['raw_unit_beta'],fitted['beta']*conversion,rtol=1e-7,atol=1e-8)
    np.testing.assert_allclose(payload['raw_unit_conditional_beta_covariance'],fitted['conditional_beta_covariance']*np.outer(conversion,conversion),rtol=1e-7,atol=1e-8)
    gradient = evaluate_gradient(CachedMatchedML(background,family,factor,design,response),np.expm1(theta))['log1p_ratio_gradient']
    projected = np.where(theta<=1e-7,np.minimum(gradient,0),np.where(theta>=upper-1e-7,np.maximum(gradient,0),gradient))
    np.testing.assert_allclose(payload['analytic_log1p_gradient'],gradient,rtol=1e-9,atol=1e-9)
    np.testing.assert_allclose(payload['projected_gradient'],projected,rtol=1e-9,atol=1e-9)
    checks = dict(best_optimizer_success=best['success'],projected_gradient_pass=bool(max(abs(projected))<=1e-3),upper_bound_contact=bool(np.any(theta>=upper-1e-6)),all_full_face_starts_agree=bool(np.ptp([c['objective'] for c in payload['candidates'] if all(c['active'])])<=1e-5))
    assert checks == payload['checks']
    passed = checks['best_optimizer_success'] and checks['projected_gradient_pass'] and not checks['upper_bound_contact'] and checks['all_full_face_starts_agree']
    assert payload['status'] == ('ml_candidate_passed_numerical_optimization_checks' if passed else 'ml_candidate_requires_optimization_review')
    np.testing.assert_array_equal(payload['zero_boundary'],theta<=1e-7)
    np.testing.assert_array_equal(payload['upper_boundary'],theta>=upper-1e-6)
    assert payload['species_kernel_is_zero'] == (not np.any(factor))
    assert payload['component_identifiability_assessed'] is False
    assert payload['component_order'] == ['background','family_component','species']
    assert payload['direct_candidate_readback']['candidates_checked'] == 22
    assert payload['direct_candidate_readback']['distinct_parameter_vectors'] == len(cache)
    return dict(candidates_checked=22,maximum_independent_objective_error=maximum_error,status=payload['status'])
