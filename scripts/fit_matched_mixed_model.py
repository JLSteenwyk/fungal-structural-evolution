#!/usr/bin/env python3
"""Optimize a specified working mixed model; keep every boundary/start result."""
import itertools
import numpy as np
from scipy.optimize import minimize
from matched_mixed_covariance import MatchedCovariance, profiled_reml


def fit_variance_ratios(background, family, factor, design, response,
                        maximum_ratio=10000., maxiter=300):
    """Profile overall scale with residual ratio fixed to one.

    Enumerate all eight zero-component faces, with three starts per nonempty
    face. Returned diagnostics are necessary checks, not a guarantee of global
    optimization or valid biological inference. Never silently expand bounds.
    """
    if not np.isfinite(maximum_ratio) or maximum_ratio <= 20 or maxiter < 1:
        raise ValueError('Need finite maximum ratio >20 and positive maxiter')
    upper = float(np.log1p(maximum_ratio))
    def evaluate(theta):
        ratios = np.expm1(theta)
        covariance = MatchedCovariance(background, family, factor, 1., *ratios)
        return profiled_reml(covariance, design, response)
    candidates = []
    for active_tuple in itertools.product([False, True], repeat=3):
        active = np.array(active_tuple)
        if not active.any():
            value = evaluate(np.zeros(3))['negative_profiled_reml']
            candidates.append(dict(active=list(active_tuple), start_ratio=0., theta=[0.,0.,0.],
                                   objective=value, success=True, message='Exact all-zero face', evaluations=1))
            continue
        for start in [.05,1.,20.]:
            def objective(parameters):
                theta = np.zeros(3)
                theta[active] = parameters
                return evaluate(theta)['negative_profiled_reml']
            result = minimize(objective, np.full(int(active.sum()),np.log1p(start)), method='L-BFGS-B',
                              bounds=[(0.,upper)]*int(active.sum()),
                              options=dict(maxiter=maxiter,ftol=1e-12,gtol=1e-7,maxls=40))
            theta = np.zeros(3); theta[active] = result.x
            candidates.append(dict(active=list(active_tuple), start_ratio=start, theta=theta.tolist(),
                                   objective=float(result.fun), success=bool(result.success),
                                   message=str(result.message), evaluations=int(result.nfev)))
    # Preserve failed attempts; a failed but lower candidate must not be hidden.
    best = min(candidates,key=lambda item:item['objective'])
    theta = np.array(best['theta'])
    fitted = evaluate(theta)
    np.testing.assert_allclose(fitted['negative_profiled_reml'],best['objective'],rtol=1e-10,atol=1e-9)
    derivatives = []
    projected = []
    for k in range(3):
        step = 1e-4
        lo,hi = theta.copy(),theta.copy()
        lo[k] = max(0.,theta[k]-step); hi[k] = min(upper,theta[k]+step)
        derivative = (evaluate(hi)['negative_profiled_reml']-evaluate(lo)['negative_profiled_reml'])/(hi[k]-lo[k])
        derivatives.append(float(derivative))
        if theta[k] <= 1e-7:
            projected.append(min(0.,float(derivative)))
        elif theta[k] >= upper-1e-7:
            projected.append(max(0.,float(derivative)))
        else:
            projected.append(float(derivative))
    boundary = theta <= 1e-7
    upper_contact = theta >= upper-1e-6
    gradient_pass = max(map(abs,projected)) <= 1e-3
    full_face = [c['objective'] for c in candidates if all(c['active'])]
    checks = dict(best_optimizer_success=best['success'], projected_gradient_pass=gradient_pass,
                  upper_bound_contact=bool(upper_contact.any()),
                  all_full_face_starts_agree=bool(np.ptp(full_face)<=1e-5))
    status = 'candidate_passed_numerical_optimization_checks' if (checks['best_optimizer_success'] and gradient_pass
              and not checks['upper_bound_contact'] and checks['all_full_face_starts_agree']) else 'candidate_requires_optimization_review'
    return dict(status=status, ratios=np.expm1(theta).tolist(), log1p_ratios=theta.tolist(),
                component_order=['background','family_component','species'], maximum_ratio=maximum_ratio,
                zero_boundary=boundary.tolist(), upper_boundary=upper_contact.tolist(),
                finite_difference_gradient=derivatives,projected_gradient=projected,checks=checks,
                candidates=candidates, **{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in fitted.items()})
