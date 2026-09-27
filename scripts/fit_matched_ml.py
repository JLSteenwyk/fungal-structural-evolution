#!/usr/bin/env python3
"""Optimize ordinary ML with analytic scores; preserve all boundary/start results."""
import itertools
import numpy as np
from scipy.optimize import minimize
from cached_matched_ml import CachedMatchedML
from matched_ml_gradient import evaluate_gradient
from collections import Counter
import time


def fit_variance_ratios(background, family, factor, design, response,
                        maximum_ratio=10000., maxiter=300):
    """Profile overall scale with residual ratio fixed to one.

    Analytic-score cancellation fails explicitly; no approximate gradient fallback.
    Enumerate all eight zero-component faces, with three starts per nonempty
    face. Returned diagnostics are necessary checks, not a guarantee of global
    optimization or valid biological inference. Never silently expand bounds.
    """
    if not np.isfinite(maximum_ratio) or maximum_ratio <= 20 or maxiter < 1:
        raise ValueError('Need finite maximum ratio >20 and positive maxiter')
    started = time.perf_counter()
    cache = CachedMatchedML(background, family, factor, design, response)
    methods = Counter()
    upper = float(np.log1p(maximum_ratio))
    def evaluate(theta):
        ratios = np.expm1(theta)
        result = cache.evaluate(*ratios)
        methods[result["evaluation_method"]] += 1
        return result
    candidates = []
    for active_tuple in itertools.product([False, True], repeat=3):
        active = np.array(active_tuple)
        if not active.any():
            value = evaluate(np.zeros(3))['negative_profiled_ml']
            candidates.append(dict(active=list(active_tuple), start_ratio=0., theta=[0.,0.,0.],
                                   objective=value, success=True, message='Exact all-zero face', evaluations=1))
            continue
        for start in [.05,1.,20.]:
            def objective(parameters):
                theta = np.zeros(3)
                theta[active] = parameters
                value = evaluate_gradient(cache, np.expm1(theta))
                methods['analytic_objective_gradient'] += 1
                return value['negative_profiled_ml'], value['log1p_ratio_gradient'][active]
            result = minimize(objective, np.full(int(active.sum()),np.log1p(start)), method='L-BFGS-B', jac=True,
                              bounds=[(0.,upper)]*int(active.sum()),
                              options=dict(maxiter=maxiter,ftol=1e-12,gtol=1e-7,maxls=40))
            theta = np.zeros(3); theta[active] = result.x
            verified = evaluate(theta)['negative_profiled_ml']
            np.testing.assert_allclose(verified, result.fun, rtol=1e-10, atol=1e-8)
            candidates.append(dict(active=list(active_tuple), start_ratio=start, theta=theta.tolist(),
                                   objective=float(result.fun), success=bool(result.success),
                                   message=str(result.message), evaluations=int(result.nfev)))
    # Preserve failed attempts; a failed but lower candidate must not be hidden.
    best = min(candidates,key=lambda item:item['objective'])
    theta = np.array(best['theta'])
    fitted = evaluate(theta)
    np.testing.assert_allclose(fitted['negative_profiled_ml'],best['objective'],rtol=1e-10,atol=1e-9)
    diagnostic = evaluate_gradient(cache, np.expm1(theta))
    derivatives = diagnostic['log1p_ratio_gradient'].tolist()
    np.testing.assert_allclose(diagnostic['negative_profiled_ml'], fitted['negative_profiled_ml'], rtol=1e-10, atol=1e-8)
    projected = []
    for k, derivative in enumerate(derivatives):
        if theta[k] <= 1e-7:
            projected.append(min(0., derivative))
        elif theta[k] >= upper-1e-7:
            projected.append(max(0., derivative))
        else:
            projected.append(derivative)
    boundary = theta <= 1e-7
    upper_contact = theta >= upper-1e-6
    gradient_pass = max(map(abs,projected)) <= 1e-3
    full_face = [c['objective'] for c in candidates if all(c['active'])]
    checks = dict(best_optimizer_success=best['success'], projected_gradient_pass=gradient_pass,
                  upper_bound_contact=bool(upper_contact.any()),
                  all_full_face_starts_agree=bool(np.ptp(full_face)<=1e-5))
    status = 'ml_candidate_passed_numerical_optimization_checks' if (checks['best_optimizer_success'] and gradient_pass
              and not checks['upper_bound_contact'] and checks['all_full_face_starts_agree']) else 'ml_candidate_requires_optimization_review'
    return dict(status=status, evaluation_methods=dict(methods), wall_seconds=time.perf_counter()-started, ratios=np.expm1(theta).tolist(), log1p_ratios=theta.tolist(),
                component_order=['background','family_component','species'], maximum_ratio=maximum_ratio,
                species_kernel_is_zero=bool(not np.any(cache.factor)),
                component_identifiability_assessed=False,
                zero_boundary=boundary.tolist(), upper_boundary=upper_contact.tolist(),
                analytic_log1p_gradient=derivatives,projected_gradient=projected,checks=checks,
                candidates=candidates, **{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in fitted.items()})
