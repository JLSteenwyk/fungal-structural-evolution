"""Independent actual-D likelihood/start/curvature/search candidate readback.

Uses frozen optimizer and comparison tolerances, retaining all failed searches
and reviews. Closed whole-grid source provenance is a separate requirement.
"""
import numpy as np
from full_expanded_model_design_sources import array_digest
from weighted_shared_entity_candidate import READY, validate_source, backend_guard, candidate as original_candidate
from independent_positive_diagonal_basis_context import IndependentPositiveDiagonalBasisContext
from readback_full_covariance_qualification import numeric as audit_qualification
from independent_shared_entity_likelihood import audit_candidate, audit_local_curvature
from independent_shared_entity_likelihood_fast import ComponentSpectralLikelihood
from independent_shared_entity_optimizer import compare_independent_optimizer


def numeric(source,plan,rows,identity,matrix,response,exported,operators=None,source_audit=None,route=None,diagonal=None):
    diagonal=validate_source(identity,source_audit,route,diagonal)
    assert np.asarray(response).shape==(len(rows),) and np.isfinite(response).all()
    assert identity['response_sha256']==array_digest(response,'<f8')
    assert operators is not None
    assert plan['independent_backend']=='component_spectral_v1'
    assert all(exported[k]==v for k,v in identity.items())
    expected_keys=set(identity)|{'disposition','fit','numerical_attempted'}
    if exported['disposition']=='shared_entity_fit_error_requires_review':expected_keys|={'error_type','error_message'}
    assert set(exported)==expected_keys
    fit=exported['fit'];source_ready=identity['source_combined_disposition']==READY
    if not source_ready:
        assert exported==dict(**identity,disposition=identity['source_combined_disposition'],fit=None,numerical_attempted=False)
        return dict(disposition='source_review_or_exclusion_retained',scientific_eligibility=False)
    assert exported['numerical_attempted'] is True
    if fit is None or 'variance_ratios' not in fit:
        # Preserve real failures without permitting an invented failure to hide
        # an original source row. Deterministic production replay verifies the
        # claimed error; it never independently promotes the failed candidate.
        replay=original_candidate(source,plan,rows,identity,matrix,response,operators,source_audit,route,diagonal)
        assert replay==exported
        return dict(disposition='original_numerical_failure_reproduced_requires_review',scientific_eligibility=False)
    assert exported['disposition']==fit['status'] and fit['method']==identity['method']
    if any(start['disposition']=='failed_search_requires_review' for start in fit['starts']):
        assert original_candidate(source,plan,rows,identity,matrix,response,operators,source_audit,route,diagonal)==exported
    fresh=backend_guard(source,rows,matrix,operators,source_audit,diagonal)
    assert fit['source_covariance_qualification']==fresh
    latent=IndependentPositiveDiagonalBasisContext(source['labels'][rows],operators,source['factors'][identity['tree']][rows]).design(matrix).audit(diagonal)
    conservative,_=audit_qualification(fresh,latent,source_audit['retained_kernel_names'],len(rows))
    assert conservative is False
    assert len(fit['starts'])==3 and fit['selected_start'] in range(3)
    independent=ComponentSpectralLikelihood(source['labels'][rows],operators,source['factors'][identity['tree']][rows],
        diagonal,matrix,response,column_batch=plan['independent_audit']['column_batch'])
    replay=audit_candidate(independent,fit,**plan['independent_audit']['replay'])
    assert replay['objective_absolute_error']<=plan['independent_audit']['replay']['objective_atol']
    assert replay['coefficient_and_conditional_covariance_comparison_passed'] and replay['profiled_scale_comparison_passed'] and replay['variance_components_comparison_passed']
    norms=np.asarray(fit['kernel_normalization']);upper=np.log1p(fit['maximum_scaled_variance'])
    assert fit['maximum_scaled_variance']==plan['optimizer']['maximum_scaled_variance']
    tolerance=64*np.finfo(float).eps*max(1.,upper);finite=[];start_checks=[]
    middle=min(1.,fit['maximum_scaled_variance']/2);rng=np.random.default_rng(20261002)
    initial_points=[np.zeros(len(norms)),np.full(len(norms),np.log1p(middle)),np.log1p(middle*rng.uniform(.1,1.8,len(norms)))]
    for number,start in enumerate(fit['starts']):
        assert start['start']==number and 0<start['search_evaluations']<=plan['optimizer']['max_evaluations']
        np.testing.assert_allclose(start['initial_coordinates'],initial_points[number],rtol=1e-13,atol=1e-14)
        if 'final_coordinates' not in start:
            assert start['disposition']=='failed_search_requires_review'
            assert start['error_type'] in ['ValueError','ArithmeticError','RuntimeError','LinAlgError']
            continue
        assert start['final_replay_evaluations']==1 and 0<=start['iterations']<=plan['optimizer']['max_iterations']
        point=np.asarray(start['final_coordinates'],dtype=float)
        assert point.shape==norms.shape and np.isfinite(point).all() and np.all((point>=0)&(point<=upper))
        value=independent.evaluate(np.expm1(point)/norms,fit['method'])
        error=abs(value['negative_profiled_likelihood']-start['objective'])
        assert error<=plan['independent_audit']['replay']['objective_atol']
        gradient=value['gradient']*np.exp(point)/norms;low=point<=tolerance;high=point>=upper-tolerance
        gradient[low]=np.minimum(gradient[low],0);gradient[high]=np.maximum(gradient[high],0)
        kkt=float(np.max(abs(gradient),initial=0))
        np.testing.assert_allclose(kkt,start['maximum_projected_gradient'],rtol=1e-3,atol=plan['independent_audit']['replay']['gradient_atol'])
        assert start['lower_boundary_indices']==np.flatnonzero(low).tolist() and start['upper_boundary_indices']==np.flatnonzero(high).tolist()
        expected='candidate_pending_independent_audit' if start['optimizer_success'] and start['maximum_projected_gradient']<=plan['optimizer']['gradient_tolerance'] and not np.any(high) else 'optimizer_or_variance_boundary_requires_review'
        assert start['disposition']==expected
        finite.append((start['objective'],number));start_checks.append(dict(start=number,objective_absolute_error=float(error),independent_projected_gradient=kkt))
    assert finite and min(finite)[1]==fit['selected_start']
    selected=fit['starts'][fit['selected_start']]
    np.testing.assert_allclose(selected['final_coordinates'],np.log1p(np.asarray(fit['variance_ratios'])*norms),rtol=1e-13,atol=1e-14)
    objectives=[f[0] for f in finite];spread=max(objectives)-min(objectives)
    np.testing.assert_allclose(fit['objective_spread'],spread,rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(fit['negative_profiled_likelihood'],min(objectives),rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(fit['multistart_objective_tolerance'],1e-7*max(1.,abs(min(objectives))),rtol=1e-13,atol=0)
    expected_status='optimized_shared_entity_candidate_pending_independent_audit' if all(s['disposition']=='candidate_pending_independent_audit' for s in fit['starts']) and spread<=fit['multistart_objective_tolerance'] else 'shared_entity_optimizer_candidate_requires_review'
    assert fit['status']==expected_status
    result=dict(disposition='producer_candidate_requires_review_independent_numeric_replay_passed',candidate_replay=replay,
        original_start_replays=start_checks,scientific_eligibility=False)
    if fit['status']=='optimized_shared_entity_candidate_pending_independent_audit':
        try:
            curvature=audit_local_curvature(independent,fit,**plan['independent_audit']['curvature'])
            search=compare_independent_optimizer(independent,fit,**plan['independent_audit']['optimizer'])
            result.update(local_curvature=curvature,independent_search=search)
            if curvature['status']=='independent_shared_entity_local_minimum_check_passed_pending_global_and_inferential_audits' and search['status']=='independent_multistart_searches_agree_pending_inferential_calibration':
                result['disposition']='independently_audited_working_candidate_pending_inferential_calibration'
            else:result['disposition']='independent_curvature_or_search_requires_review'
        except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
            result.update(disposition='independent_numerical_precision_requires_review',error_type=type(error).__name__,error_message=str(error))
    return result

