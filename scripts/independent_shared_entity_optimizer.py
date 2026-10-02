"""Independent bounded SLSQP searches using spectral likelihood and scores.

Agreement supports numerical candidate qualification. It is not a global
optimality proof, a likelihood test or calibrated biological inference.
"""
import numpy as np
from scipy.optimize import minimize
from independent_shared_entity_likelihood import IndependentSharedEntityLikelihood,audit_candidate


def compare_independent_optimizer(likelihood,candidate,max_iterations=200,max_evaluations=500,
                                 gradient_tolerance=1e-6,objective_tolerance=1e-7):
    if not isinstance(likelihood,IndependentSharedEntityLikelihood):raise TypeError('Independent spectral likelihood required')
    if not isinstance(max_iterations,int) or max_iterations<1 or not isinstance(max_evaluations,int) or max_evaluations<1:
        raise ValueError('Positive independent search budgets required')
    if any(not np.isfinite(v) or v<=0 for v in [gradient_tolerance,objective_tolerance]):
        raise ValueError('Positive finite independent tolerances required')
    replay=audit_candidate(likelihood,candidate,objective_atol=objective_tolerance,gradient_atol=gradient_tolerance)
    if replay['status']!='independent_shared_entity_replay_passed_pending_curvature_and_calibration':
        return dict(status='independent_shared_entity_search_requires_review',candidate_replay=replay,
            starts=[],scientific_eligibility=False,reason='Supplied candidate failed independent replay')
    norms=np.asarray(candidate['kernel_normalization']);cap=float(candidate['maximum_scaled_variance'])
    upper=float(np.log1p(cap));boundary_tolerance=64*np.finfo(float).eps*max(1.,upper)
    middle=min(1.,cap/2);rng=np.random.default_rng(191943)
    starts=[np.zeros(len(norms)),np.full(len(norms),np.log1p(middle)),
        np.log1p(middle*rng.uniform(.2,1.7,len(norms)))]
    histories=[]
    for number,start in enumerate(starts):
        evaluations=0;last=None
        def objective(point):
            nonlocal evaluations,last
            if evaluations>=max_evaluations:raise RuntimeError('Independent evaluation budget exhausted')
            evaluations+=1
            if not np.isfinite(point).all() or np.any(point<0) or np.any(point>upper):
                raise ArithmeticError('Independent optimizer proposed invalid coordinates')
            value=likelihood.evaluate(np.expm1(point)/norms,candidate['method'])
            last=value['negative_profiled_likelihood']
            return last,value['gradient']*np.exp(point)/norms
        try:
            fitted=minimize(objective,start,method='SLSQP',jac=True,bounds=[(0.,upper)]*len(norms),
                options=dict(maxiter=max_iterations,ftol=1e-12))
            point=np.asarray(fitted.x);value=likelihood.evaluate(np.expm1(point)/norms,candidate['method'])
            gradient=value['gradient']*np.exp(point)/norms
            low=point<=boundary_tolerance;high=point>=upper-boundary_tolerance;projected=gradient.copy()
            projected[low]=np.minimum(projected[low],0);projected[high]=np.maximum(projected[high],0)
            kkt=float(np.max(abs(projected),initial=0))
            error=abs(value['negative_profiled_likelihood']-candidate['negative_profiled_likelihood'])
            accepted=bool(fitted.success) and kkt<=gradient_tolerance and not np.any(high) and error<=objective_tolerance
            summary=dict(start=number,initial_coordinates=start.tolist(),final_coordinates=point.tolist(),
                objective=float(value['negative_profiled_likelihood']),supplied_candidate_objective_absolute_difference=float(error),
                optimizer_success=bool(fitted.success),optimizer_message=str(fitted.message),iterations=int(fitted.nit),
                search_evaluations=evaluations,final_replay_evaluations=1,maximum_projected_gradient=kkt,
                lower_boundary_indices=np.flatnonzero(low).tolist(),upper_boundary_indices=np.flatnonzero(high).tolist(),
                disposition='independent_search_agrees_pending_calibration' if accepted else 'independent_search_requires_review')
        except (ValueError,ArithmeticError,RuntimeError,np.linalg.LinAlgError) as error:
            summary=dict(start=number,initial_coordinates=start.tolist(),search_evaluations=evaluations,
                disposition='failed_independent_search_requires_review',error_type=type(error).__name__,error_message=str(error))
            if last is not None:summary['last_finite_objective']=float(last)
        histories.append(summary)
    passed=all(r['disposition']=='independent_search_agrees_pending_calibration' for r in histories)
    return dict(status='independent_multistart_searches_agree_pending_inferential_calibration' if passed
        else 'independent_shared_entity_search_requires_review',candidate_replay=replay,starts=histories,
        method='SLSQP',objective_absolute_tolerance=objective_tolerance,gradient_absolute_tolerance=gradient_tolerance,
        maximum_search_evaluations_per_start=max_evaluations,maximum_iterations_per_start=max_iterations,
        scientific_eligibility=False,scope='Three deterministic spectral SLSQP searches independent of '
            'production Cholesky/L-BFGS-B. All starts must meet boundary KKT checks and agree with the '
            'supplied candidate at the explicit absolute likelihood tolerance. Failed, disagreeing, '
            'cap-limited or budget-limited searches remain review; no global optimum proof or inference.')
