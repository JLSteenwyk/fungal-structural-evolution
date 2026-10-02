"""Deterministic multi-start optimization of qualified shared-entity variances.

Returned candidates require independent optimization/numerical audit and
scientific calibration. Finite ratio caps are explicit review boundaries.
"""
import numpy as np
from scipy.optimize import minimize
from covariance_basis_context import ComponentKernelProducts
from shared_entity_likelihood import SharedEntityLikelihood


def fit_shared_entity(likelihood,method='reml',maximum_scaled_variance=1e6,
                      max_iterations=200,max_evaluations=500,gradient_tolerance=1e-6):
    if not isinstance(likelihood,SharedEntityLikelihood):raise TypeError('Explicit shared-entity likelihood required')
    if method not in ['ml','reml']:raise ValueError('Likelihood must be ml or reml')
    if not np.isfinite(maximum_scaled_variance) or maximum_scaled_variance<=0:raise ValueError('Positive finite variance cap required')
    if not isinstance(max_iterations,int) or max_iterations<1 or not isinstance(max_evaluations,int) or max_evaluations<1:
        raise ValueError('Positive iteration and evaluation budgets required')
    if not np.isfinite(gradient_tolerance) or gradient_tolerance<=0:raise ValueError('Positive gradient tolerance required')
    qualification=ComponentKernelProducts(likelihood.labels,likelihood.incidence,likelihood.diagonal).tree(likelihood.factor).audit(likelihood.design)
    if any(qualification[k]['disposition']!='numerically_independent_covariance_bases' for k in ['raw_diagnostics','reml_diagnostics']):
        raise ValueError('Covariance bases require identifiability review before variance fitting')
    norms=np.asarray([float(z.multiply(z).sum())/likelihood.n for z in likelihood.incidence.values()]
        +[float(np.sum(likelihood.factor*likelihood.factor))/likelihood.n])/float(np.mean(likelihood.diagonal))
    if np.any(norms<=0) or not np.isfinite(norms).all():raise ValueError('Positive finite kernel norms required')
    upper=float(np.log1p(maximum_scaled_variance));tolerance=64*np.finfo(float).eps*max(1.,upper)
    rng=np.random.default_rng(20261002);middle=min(1.,maximum_scaled_variance/2)
    starts=[np.zeros(len(norms)),np.full(len(norms),np.log1p(middle)),
        np.log1p(middle*rng.uniform(.1,1.8,len(norms)))]
    histories=[];candidates=[]
    for number,start in enumerate(starts):
        evaluations=0;last=None
        def objective(point):
            nonlocal evaluations,last
            if evaluations>=max_evaluations:raise RuntimeError('Explicit likelihood evaluation budget exhausted')
            evaluations+=1
            if not np.isfinite(point).all() or np.any(point<0) or np.any(point>upper):raise ArithmeticError('Optimizer proposed an invalid bounded coordinate')
            ratios=np.expm1(point)/norms;value=likelihood.evaluate(ratios,method)
            gradient=value['gradient']*np.exp(point)/norms
            last=(point.copy(),value,gradient.copy())
            return value['negative_profiled_likelihood'],gradient
        try:
            optimized=minimize(objective,start,method='L-BFGS-B',jac=True,bounds=[(0.,upper)]*len(norms),
                options=dict(maxiter=max_iterations,maxfun=max_evaluations,gtol=gradient_tolerance,ftol=1e-12,maxls=30))
            # Final replay is explicit and does not consume the search budget.
            point=np.asarray(optimized.x);value=likelihood.evaluate(np.expm1(point)/norms,method)
            gradient=value['gradient']*np.exp(point)/norms
            projected=gradient.copy();low=point<=tolerance;high=point>=upper-tolerance
            projected[low]=np.minimum(projected[low],0);projected[high]=np.maximum(projected[high],0)
            kkt=float(np.max(abs(projected),initial=0))
            summary=dict(start=number,initial_coordinates=start.tolist(),final_coordinates=point.tolist(),
                objective=float(value['negative_profiled_likelihood']),optimizer_success=bool(optimized.success),
                optimizer_message=str(optimized.message),search_evaluations=evaluations,final_replay_evaluations=1,
                iterations=int(optimized.nit),maximum_projected_gradient=kkt,
                lower_boundary_indices=np.flatnonzero(low).tolist(),upper_boundary_indices=np.flatnonzero(high).tolist(),
                disposition='candidate_pending_independent_audit' if optimized.success and kkt<=gradient_tolerance and not np.any(high)
                    else 'optimizer_or_variance_boundary_requires_review')
            candidates.append((summary,value))
        except (ValueError,ArithmeticError,RuntimeError,np.linalg.LinAlgError) as error:
            summary=dict(start=number,initial_coordinates=start.tolist(),search_evaluations=evaluations,
                disposition='failed_search_requires_review',error_type=type(error).__name__,error_message=str(error))
            if last is not None:summary['last_finite_objective']=float(last[1]['negative_profiled_likelihood'])
        histories.append(summary)
    if not candidates:
        return dict(status='all_shared_entity_optimizer_searches_require_review',method=method,starts=histories,
            parameter_names=likelihood.parameter_names,kernel_normalization=norms.tolist(),
            maximum_scaled_variance=maximum_scaled_variance,scientific_eligibility=False)
    best,value=min(candidates,key=lambda v:v[0]['objective'])
    spread=max(c[0]['objective'] for c in candidates)-best['objective']
    objective_tolerance=1e-7*max(1.,abs(best['objective']))
    all_qualified=all(h['disposition']=='candidate_pending_independent_audit' for h in histories)
    status='optimized_shared_entity_candidate_pending_independent_audit' if all_qualified and spread<=objective_tolerance else 'shared_entity_optimizer_candidate_requires_review'
    return dict(status=status,method=method,starts=histories,selected_start=best['start'],objective_spread=float(spread),
        multistart_objective_tolerance=float(objective_tolerance),kernel_normalization=norms.tolist(),
        maximum_scaled_variance=maximum_scaled_variance,parameter_names=likelihood.parameter_names,
        variance_ratios=value['variance_ratios'].tolist(),variance_components=value['variance_components'].tolist(),
        profiled_scale=value['profiled_scale'],negative_profiled_likelihood=value['negative_profiled_likelihood'],
        beta=value['beta'].tolist(),conditional_beta_covariance=value['conditional_beta_covariance'].tolist(),
        coefficient_covariance_is_conditional=True,source_covariance_qualification=qualification,
        scientific_eligibility=False,scope='Numerical optimizer candidate only. Nonnegative expm1 coordinates retain '
            'exact zero boundaries. Every search and computational cap is recorded; failed/disagreeing starts '
            'or upper-bound solutions remain review states. Independent replay, curvature/numerical audit and '
            'calibrated inferential uncertainty remain required; no biological effect is accepted.')
