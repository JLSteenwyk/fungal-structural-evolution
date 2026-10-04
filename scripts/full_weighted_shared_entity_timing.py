"""Whole-grid timing census and deterministic qualified-input probes.

Every candidate is retained. Probes characterize computation rather than fit
an evolutionary effect; their extrapolations are conditional planning numbers.
"""
from collections import Counter
import time
import numpy as np
from full_expanded_model_design_sources import array_digest, digest
from full_weighted_shared_entity_fit_sources import cases
from weighted_shared_entity_candidate import READY, validate_source, backend_guard
from independent_positive_diagonal_basis_context import IndependentPositiveDiagonalBasisContext
from readback_full_covariance_qualification import numeric as audit_qualification
from independent_shared_entity_likelihood import audit_candidate
from independent_shared_entity_likelihood_fast import ComponentSpectralLikelihood
from shared_entity_likelihood import SharedEntityLikelihood

SCHEMA='full-four-control-qualified-input-timing-v1'
PRODUCER='complete_full_four_control_input_timing_pending_accounting_readback_v1'
READER='passed_full_four_control_input_timing_census_numeric_and_planning_readback_v1'


def groups(source,fit_plan,cohort,entries):
    selected={};counts=Counter();census=[]
    for identity,x,y,audit,route,diagonal in cases(source,fit_plan,cohort,entries):
        key=(identity['control_policy'],identity['loading_mode'],identity['tree'],identity['method'],identity['outcome'])
        group_id=digest([SCHEMA,source['fit_contract'],cohort['cohort_id'],*key])
        ready=identity['source_combined_disposition']==READY
        record=dict(candidate_id=identity['candidate_id'],identity_sha256=digest(identity),
            source_disposition=identity['source_combined_disposition'],group_id=group_id if ready else None,
            active_columns=x.shape[1],scientific_eligibility=False)
        census.append(record)
        if not ready:continue
        counts[group_id]+=1
        condition=float(audit['numerical_audit']['normalized_design_condition_number'])
        assert np.isfinite(condition) and condition>=1 and x.shape[1]<len(y)
        rank=(x.shape[1],condition,identity['candidate_id'])
        if group_id not in selected or rank>selected[group_id]['rank']:
            selected[group_id]=dict(rank=rank,identity=identity,matrix=x,response=y,
                source_audit=audit,route=route,diagonal=diagonal,key=key)
    assert len({r['candidate_id'] for r in census})==len(census)
    for row in census:
        row['benchmark_candidate_id']=selected[row['group_id']]['identity']['candidate_id'] if row['group_id'] else None
    return census,selected,counts


class TimedSpectralLikelihood(ComponentSpectralLikelihood):
    def evaluate(self,*args,**kwargs):
        started=time.perf_counter()
        try:
            self.last_result=super().evaluate(*args,**kwargs)
            return self.last_result
        finally:self.last_seconds=time.perf_counter()-started


def probe(source,fit_plan,timing_plan,rows,operators,representative,count,group_id):
    identity=representative['identity'];x=representative['matrix'];y=representative['response']
    result=dict(group_id=group_id,representative=identity,eligible_candidates=count,
        selection_active_columns=x.shape[1],selection_condition=representative['rank'][1],
        variance_points=timing_plan['scaled_variance_points'],scientific_eligibility=False)
    started=time.perf_counter()
    try:
        validate_source(identity,representative['source_audit'],representative['route'],representative['diagonal'])
        assert array_digest(y,'<f8')==identity['response_sha256']
        result['source_validation_seconds']=time.perf_counter()-started
        started=time.perf_counter()
        primary=SharedEntityLikelihood(source['labels'][rows],operators,source['factors'][identity['tree']][rows],representative['diagonal'],x,y)
        result['primary_constructor_seconds']=time.perf_counter()-started
        started=time.perf_counter()
        qualification=backend_guard(source,rows,x,operators,representative['source_audit'],representative['diagonal'])
        result['production_qualification_guard_seconds']=time.perf_counter()-started
        # Actual-D closed envelopes remain bound and unchanged.
        # The optimizer backend receives its own fresh additional guard.
        assert qualification['kernel_names']==representative['source_audit']['retained_kernel_names']
        assert all(qualification[k]['disposition']=='numerically_independent_covariance_bases' for k in ['raw_diagnostics','reml_diagnostics'])
        norms=np.asarray([float(z.multiply(z).sum())/primary.n for z in primary.incidence.values()]
            +[float(np.sum(primary.factor*primary.factor))/primary.n])/float(np.mean(primary.diagonal))
        assert np.isfinite(norms).all() and np.all(norms>0)
        started=time.perf_counter()
        fresh=backend_guard(source,rows,x,operators,representative['source_audit'],representative['diagonal'])
        latent=IndependentPositiveDiagonalBasisContext(primary.labels,operators,primary.factor).design(x).audit(primary.diagonal)
        conservative,_=audit_qualification(fresh,latent,representative['source_audit']['retained_kernel_names'],len(rows))
        assert conservative is False
        result['reader_qualification_seconds']=time.perf_counter()-started
        started=time.perf_counter()
        independent=TimedSpectralLikelihood(primary.labels,operators,primary.factor,primary.diagonal,x,y,
            column_batch=fit_plan['independent_audit']['column_batch'])
        result['independent_constructor_seconds']=time.perf_counter()-started
        result.update(kernel_normalization=norms.tolist(),parameter_names=primary.parameter_names,points=[])
        for point in timing_plan['scaled_variance_points']:
            assert 0<=point<=fit_plan['optimizer']['maximum_scaled_variance']
            ratios=np.full(len(norms),point)/norms;started=time.perf_counter()
            try:
                value=primary.evaluate(ratios,identity['method']);primary_seconds=time.perf_counter()-started
                exported=dict(value,kernel_normalization=norms,maximum_scaled_variance=fit_plan['optimizer']['maximum_scaled_variance'],
                    status='timing_probe_not_optimized',scientific_eligibility=False)
                audit=audit_candidate(independent,exported,**fit_plan['independent_audit']['replay'])
                coordinate_error=float(np.max(abs(value['gradient']-independent.last_result['gradient'])*(1+point)/norms,initial=0))
                agreement=(audit['objective_absolute_error']<=fit_plan['independent_audit']['replay']['objective_atol'] and
                    coordinate_error<=fit_plan['independent_audit']['replay']['gradient_atol'] and
                    all(audit[k] for k in ['coefficient_and_conditional_covariance_comparison_passed',
                        'profiled_scale_comparison_passed','variance_components_comparison_passed']))
                result['points'].append(dict(scaled_variance=point,
                    status='timed_likelihood_agreement_only' if agreement else 'timing_probe_numerical_agreement_requires_review',
                    primary_evaluation_seconds=primary_seconds,independent_evaluation_seconds=independent.last_seconds,
                    objective_absolute_error=audit['objective_absolute_error'],maximum_coordinate_gradient_error=coordinate_error,
                    coefficient_and_conditional_covariance_comparison_passed=audit['coefficient_and_conditional_covariance_comparison_passed'],
                    profiled_scale_comparison_passed=audit['profiled_scale_comparison_passed'],
                    variance_components_comparison_passed=audit['variance_components_comparison_passed'],
                    independent_inverse_relative_residual=audit['inverse_relative_residual'],
                    explicit_whitening_fallback_batches=independent.last_result['explicit_global_whitening_fallback_batches'],
                    independent_cached_kernel_bytes=independent.last_result['cached_component_kernel_bytes']))
            except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
                result['points'].append(dict(scaled_variance=point,status='timing_probe_precision_requires_review',
                    error_type=type(error).__name__,error_message=str(error),elapsed_seconds=time.perf_counter()-started))
        result['status']='timed_all_declared_points_agree_only' if all(
            p['status']=='timed_likelihood_agreement_only' for p in result['points']) else 'timing_group_requires_review'
    except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
        result.update(status='timing_group_construction_requires_review',error_type=type(error).__name__,error_message=str(error))
    return result


def estimate(probes,fit_plan):
    # A finite production fit with a failed search is first reproduced and
    # then receives up to four spectral replays. A fully converged fit instead
    # receives curvature and three independent searches. Keep both paths.
    primary_budget=3*(fit_plan['optimizer']['max_evaluations']+1)
    independent_search=3*(fit_plan['independent_audit']['optimizer']['max_evaluations']+1)
    summaries=[];unmeasured=measured=0
    for p in probes:
        if p['status']!='timed_all_declared_points_agree_only':unmeasured+=p['eligible_candidates'];continue
        r=len(p['parameter_names']);independent_budget=5+2+4*r+independent_search
        primary_time=max(v['primary_evaluation_seconds'] for v in p['points'])
        independent_time=max(v['independent_evaluation_seconds'] for v in p['points'])
        primary_setup=p['source_validation_seconds']+p['primary_constructor_seconds']+2*p['production_qualification_guard_seconds']
        independent_setup=p['source_validation_seconds']+p['reader_qualification_seconds']+p['independent_constructor_seconds']
        producer=primary_setup+primary_budget*primary_time
        reader=max(primary_setup+primary_budget*primary_time+independent_setup+4*independent_time,
                   independent_setup+independent_budget*independent_time)
        summaries.append(dict(group_id=p['group_id'],eligible_candidates=p['eligible_candidates'],
            primary_evaluation_budget=primary_budget,independent_evaluation_budget=independent_budget,
            conditional_budget_weighted_producer_seconds=p['eligible_candidates']*producer,
            conditional_budget_weighted_reader_seconds=p['eligible_candidates']*reader))
        measured+=p['eligible_candidates']
    return dict(measured_candidate_coverage=measured,unmeasured_review_candidate_coverage=unmeasured,
        group_planning_costs=summaries,
        conditional_budget_weighted_seconds=sum(p['conditional_budget_weighted_producer_seconds']+
            p['conditional_budget_weighted_reader_seconds'] for p in summaries),
        production_finish_eta=None,mathematical_runtime_bound=False,
        numerical_probe_agreement_is_not_fit_acceptance=True,
        scope='Maximum observed time among the declared variance points, multiplied by original explicit evaluation budgets and eligible candidate counts. '
            'Representative design/response timings do not bound all optimizer points or infer actual evaluation counts. '
            'Source validation, both fresh production guards, reader fresh-plus-latent qualification, failure reproduction and independent optimization paths are budgeted. Loading/export/SQL/closure, '
            'unmeasured review groups are excluded; all four weighting controls are retained. Independent timing readback costs are separate from future fit budgets; no project or production-finish ETA.')
