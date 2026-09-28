#!/usr/bin/env python3
"""Separate zero-species-shortcut refinement; identical candidate grid and checks."""
import itertools
import time
import numpy as np
from scipy.optimize import minimize
from cached_matched_likelihood import CachedMatchedLikelihood
from matched_reml_gradient_fast import evaluate_gradient
from matched_mixed_covariance import MatchedCovariance, profiled_reml


def refine(background, family, factor, design, response, original_theta,
           maximum_ratio=10000., maxiter=1000):
    started = time.perf_counter()
    if not np.isfinite(maximum_ratio) or maximum_ratio <= 20 or maxiter < 1:
        raise ValueError('Invalid bounds or iteration limit')
    upper = float(np.log1p(maximum_ratio))
    initial = np.asarray(original_theta, dtype=float)
    if initial.shape != (3,) or not np.isfinite(initial).all() or np.any(initial < 0) or np.any(initial > upper):
        raise ValueError('Original parameters outside fixed bounds')
    cache = CachedMatchedLikelihood(background, family, factor, design, response)
    candidates = []

    def record(theta, active, start, success, message, evaluations):
        result = evaluate_gradient(cache, np.expm1(theta))
        candidates.append(dict(theta=np.asarray(theta).tolist(), active=list(map(bool, active)),
                               start=start, objective=result['negative_profiled_reml'],
                               success=bool(success), message=str(message), evaluations=int(evaluations)))

    record(initial, [True]*3, 'original_unoptimized_reference', False, 'Retained reference, not an optimizer termination', 1)
    for face in itertools.product([False, True], repeat=3):
        active = np.array(face)
        if not active.any():
            record(np.zeros(3), active, 'exact_zero', True, 'Exact all-zero face', 1)
            continue
        starts = [(str(value), np.full(active.sum(), np.log1p(value))) for value in [.05, 1., 20.]]
        if active.all():
            starts.append(('original_parameters', initial.copy()))
        for label, start in starts:
            def objective(parameters):
                theta = np.zeros(3)
                theta[active] = parameters
                result = evaluate_gradient(cache, np.expm1(theta))
                return result['negative_profiled_reml'], result['log1p_ratio_gradient'][active]
            fit = minimize(objective, start, jac=True, method='L-BFGS-B',
                           bounds=[(0., upper)]*int(active.sum()),
                           options=dict(maxiter=maxiter, ftol=1e-14, gtol=1e-8, maxls=80))
            theta = np.zeros(3)
            theta[active] = fit.x
            record(theta, active, label, fit.success, fit.message, fit.nfev)
    # Successful exact ties take precedence, but no lower failed candidate is hidden.
    best = min(candidates, key=lambda c: (c['objective'], not c['success']))
    theta = np.array(best['theta'])
    analytic = evaluate_gradient(cache, np.expm1(theta))
    gradient = analytic['log1p_ratio_gradient']
    projected = np.where(theta <= 1e-7, np.minimum(gradient, 0),
                         np.where(theta >= upper-1e-7, np.maximum(gradient, 0), gradient))
    maximum_error = 0.
    for candidate in candidates:
        direct = profiled_reml(MatchedCovariance(background, family, factor, 1., *np.expm1(candidate['theta'])), design, response)
        reference = direct['negative_profiled_reml']
        np.testing.assert_allclose(reference, candidate['objective'], rtol=1e-9, atol=1e-7)
        maximum_error = max(maximum_error, abs(reference-candidate['objective']))
    fitted = profiled_reml(MatchedCovariance(background, family, factor, 1., *np.expm1(theta)), design, response)
    full = [c['objective'] for c in candidates if all(c['active']) and c['start'] != 'original_unoptimized_reference']
    checks = dict(best_optimizer_success=best['success'], projected_gradient_pass=bool(np.max(abs(projected)) <= 1e-3),
                  upper_bound_contact=bool(np.any(theta >= upper-1e-6)), all_full_face_starts_agree=bool(np.ptp(full) <= 1e-5),
                  original_objective_not_worsened=bool(best['objective'] <= candidates[0]['objective']))
    passed = all(checks[k] for k in checks if k != 'upper_bound_contact') and not checks['upper_bound_contact']
    return dict(status='refinement_passed_numerical_checks' if passed else 'refinement_requires_review',
                checks=checks, ratios=np.expm1(theta).tolist(), log1p_ratios=theta.tolist(), maximum_ratio=maximum_ratio,
                analytic_gradient=gradient.tolist(), projected_gradient=projected.tolist(), candidates=candidates,
                original_objective=candidates[0]['objective'], objective_improvement=candidates[0]['objective']-best['objective'],
                maximum_direct_candidate_objective_error=maximum_error, wall_seconds=time.perf_counter()-started,
                **{k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in fitted.items()},
                scope='Numerical refinement only; not proof of global optimality, identifiable components, calibrated inference or biological adequacy.')
