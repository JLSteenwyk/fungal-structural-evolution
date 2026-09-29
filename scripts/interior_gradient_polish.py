#!/usr/bin/env python3
"""Conservative local proposals; never replace a fit or relax its checks."""
import numpy as np


def propose(objective, gradient, theta, upper):
    """Require two positive-curvature proposals and direct derivative agreement.

    Bounds are log1p variance-ratio bounds. This routine applies only to
    interior points. Its result is a local numerical diagnostic, not optimizer
    convergence evidence for other starts or a global-optimum certificate.
    """
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (3,) or not np.isfinite(theta).all() or not np.isfinite(upper) or upper <= 0:
        raise ValueError('Expected three finite parameters and a positive bound')
    if (theta < 0).any() or (theta > upper).any():
        raise ValueError('Parameters outside fixed bounds')
    margin = 1e-5
    if (theta <= margin).any() or (theta >= upper - margin).any():
        return dict(status='not_applicable_near_boundary', proposals=[])
    base = float(objective(theta.copy()))
    g = np.asarray(gradient(theta.copy()), dtype=float)
    if not np.isfinite(base) or g.shape != (3,) or not np.isfinite(g).all():
        raise ValueError('Nonfinite objective or invalid gradient')
    if np.max(abs(g)) <= 1e-3:
        return dict(status='no_gradient_correction_needed', proposals=[], original_objective=base, original_gradient=g.tolist())
    trials = []
    for step in [1e-5, 1e-6]:
        columns = []
        for axis in range(3):
            offset = np.zeros(3)
            offset[axis] = step
            high = np.asarray(gradient(theta + offset), dtype=float)
            low = np.asarray(gradient(theta - offset), dtype=float)
            if high.shape != (3,) or low.shape != (3,) or not np.isfinite([high, low]).all():
                raise ValueError('Invalid finite-difference gradient')
            columns.append((high - low) / (2 * step))
        raw = np.column_stack(columns)
        hessian = (raw + raw.T) / 2
        eigenvalues = np.linalg.eigvalsh(hessian)
        trial = dict(hessian_step=step, hessian=hessian.tolist(), eigenvalues=eigenvalues.tolist(), passed=False)
        trials.append(trial)
        if not np.isfinite(eigenvalues).all() or not (eigenvalues > 0).all():
            trial['reason'] = 'nonpositive_or_nonfinite_curvature'
            continue
        candidate = theta - np.linalg.solve(hessian, g)
        trial['theta'] = candidate.tolist()
        if not np.isfinite(candidate).all() or (candidate <= 1e-6).any() or (candidate >= upper - 1e-6).any():
            trial['reason'] = 'proposal_not_interior_for_direct_difference_checks'
            continue
        value = float(objective(candidate.copy()))
        derivatives = np.asarray(gradient(candidate.copy()), dtype=float)
        if derivatives.shape != (3,) or not np.isfinite(value) or not np.isfinite(derivatives).all():
            raise ValueError('Invalid proposal objective or gradient')
        differences = []
        for h in [1e-6, 1e-7]:
            estimate = []
            for axis in range(3):
                offset = np.zeros(3)
                offset[axis] = h
                estimate.append((objective(candidate + offset) - objective(candidate - offset)) / (2 * h))
            if not np.isfinite(estimate).all():
                raise ValueError('Nonfinite direct finite difference')
            differences.append(dict(step=h, gradient=list(map(float, estimate))))
        checks = dict(objective_not_worsened=bool(value <= base),
                      analytic_gradient_pass=bool(np.max(abs(derivatives)) <= 1e-3),
                      both_direct_gradient_checks_pass=all(max(map(abs, d['gradient'])) <= 1e-3 for d in differences))
        trial.update(objective=value, gradient=derivatives.tolist(), finite_differences=differences,
                     checks=checks, passed=all(checks.values()))
    passed = all(t['passed'] for t in trials)
    return dict(status='two_local_proposals_passed' if passed else 'local_polish_requires_review',
                original_theta=theta.tolist(), original_objective=base, original_gradient=g.tolist(),
                upper=float(upper), proposals=trials,
                scope='Local proposals only. Existing1e-3 criterion unchanged; both Hessian-step proposals must pass. Original all-face/start diagnostics and full output audit remain required. No fit replacement.')
