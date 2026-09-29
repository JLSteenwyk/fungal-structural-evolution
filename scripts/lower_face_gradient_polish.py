#!/usr/bin/env python3
"""Local proposals preserving exact zero variance components and their KKT signs."""
import numpy as np


def propose(objective, gradient, theta, upper):
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (3,) or not np.isfinite(theta).all() or not np.isfinite(upper) or upper <= 0:
        raise ValueError('Expected three finite parameters and positive upper bound')
    if np.any(theta < 0) or np.any(theta > upper):
        raise ValueError('Parameters outside bounds')
    fixed, free = np.flatnonzero(theta == 0), np.flatnonzero(theta > 0)
    if not len(fixed) or not len(free):
        return dict(status='not_applicable_without_mixed_lower_face', proposals=[])
    if np.any(theta[free] <= 1e-5) or np.any(theta[free] >= upper - 1e-5):
        return dict(status='not_applicable_near_other_boundary', proposals=[])
    def grad(point):
        value = np.asarray(gradient(point), dtype=float)
        if value.shape != (3,) or not np.isfinite(value).all():
            raise ValueError('Invalid gradient')
        return value
    def obj(point):
        value = float(objective(point))
        if not np.isfinite(value):
            raise ValueError('Invalid objective')
        return value
    base, g = obj(theta), grad(theta)
    if np.any(g[fixed] < 0):
        return dict(status='boundary_sign_requires_review', proposals=[], original_gradient=g.tolist())
    if np.max(abs(g[free])) <= 1e-3:
        return dict(status='no_free_gradient_correction_needed', proposals=[], original_gradient=g.tolist())
    trials = []
    for h in (1e-5, 1e-6):
        raw = np.column_stack([(grad(theta + np.eye(3)[j] * h)[free] - grad(theta - np.eye(3)[j] * h)[free]) / (2 * h) for j in free])
        hessian = (raw + raw.T) / 2
        eigenvalues = np.linalg.eigvalsh(hessian)
        trial = dict(hessian_step=h, hessian=hessian.tolist(), eigenvalues=eigenvalues.tolist(), passed=False)
        trials.append(trial)
        if np.any(eigenvalues <= 0) or not np.isfinite(eigenvalues).all():
            trial['reason'] = 'nonpositive_or_nonfinite_curvature'
            continue
        point = theta.copy()
        point[free] -= np.linalg.solve(hessian, g[free])
        trial['theta'] = point.tolist()
        if not np.isfinite(point).all() or np.any(point[free] <= 1e-6) or np.any(point[free] >= upper - 1e-6):
            trial['reason'] = 'proposal_not_interior_on_free_face'
            continue
        value, derivative = obj(point), grad(point)
        differences = []
        for step in (1e-6, 1e-7):
            estimates = np.zeros(3)
            for axis in range(3):
                offset = np.eye(3)[axis] * step
                if axis in fixed:
                    estimates[axis] = (-3 * value + 4 * obj(point + offset) - obj(point + 2 * offset)) / (2 * step)
                else:
                    estimates[axis] = (obj(point + offset) - obj(point - offset)) / (2 * step)
            differences.append(dict(step=step, gradient=estimates.tolist(), free_gradient_pass=bool(np.max(abs(estimates[free])) <= 1e-3), boundary_sign_pass=bool(np.all(estimates[fixed] >= 0))))
        checks = dict(objective_not_worsened=bool(value <= base), analytic_free_gradient_pass=bool(np.max(abs(derivative[free])) <= 1e-3), analytic_boundary_sign_pass=bool(np.all(derivative[fixed] >= 0)), direct_free_gradients_pass=all(d['free_gradient_pass'] for d in differences), direct_boundary_signs_pass=all(d['boundary_sign_pass'] for d in differences))
        trial.update(objective=value, gradient=derivative.tolist(), finite_differences=differences, checks=checks, passed=all(checks.values()))
    return dict(status='two_lower_face_proposals_passed' if all(t['passed'] for t in trials) else 'lower_face_proposals_require_review', original_theta=theta.tolist(), original_objective=base, original_gradient=g.tolist(), fixed_axes=fixed.tolist(), free_axes=free.tolist(), upper=float(upper), proposals=trials, scope='Local proposals only; exact lower face preserved. Both derivative scales and both Hessian steps required. No upper-bound handling, fit replacement, other-flag resolution, global optimum or calibrated inference.')
