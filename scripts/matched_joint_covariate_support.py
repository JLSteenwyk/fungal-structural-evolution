"""Numerical convex-hull support for the uncentered zero-covariate reference."""
import numpy as np
from scipy.optimize import linprog


def assess_support(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or not len(x) or not np.isfinite(x).all():
        raise ValueError('Invalid covariate matrix')
    scales = np.max(abs(x), axis=0)
    scales[scales == 0] = 1.
    z = x / scales
    n, d = z.shape
    # Minimize t with |sum_i w_i z_i| <= t, sum_i w_i=1, w_i>=0.
    objective = np.r_[np.zeros(n), 1.]
    inequalities = np.column_stack([np.vstack([z.T, -z.T]), -np.ones(2*d)])
    equality = np.r_[np.ones(n), 0.][None, :]
    fit = linprog(objective, A_ub=inequalities, b_ub=np.zeros(2*d),
                  A_eq=equality, b_eq=[1.], bounds=(0, None), method='highs-ds',
                  options={'primal_feasibility_tolerance': 1e-9,
                           'dual_feasibility_tolerance': 1e-9, 'time_limit': 30.})
    result = dict(records=n, scales=scales.tolist(), solver_status=int(fit.status),
                  solver_message=fit.message, classification='unresolved_solver')
    if not fit.success:
        return result
    weights = fit.x[:n]
    barycenter = weights @ z
    distance = float(max(abs(barycenter)))
    direction = fit.ineqlin.marginals[d:] - fit.ineqlin.marginals[:d]
    norm = float(sum(abs(direction)))
    lower = float(min(z @ direction))
    valid = bool(np.min(weights) >= 0 and abs(sum(weights)-1) <= 1e-8
                 and abs(distance-fit.fun) <= 1e-8 and norm <= 1+1e-8
                 and lower <= distance+1e-8)
    status = 'unresolved_certificate'
    if valid:
        if distance <= 1e-8:
            status = 'zero_supported_to_numeric_tolerance'
        elif lower > 1e-7:
            status = 'zero_outside_joint_convex_hull'
        else:
            status = 'unresolved_near_boundary'
    indices = np.flatnonzero(weights != 0)
    result.update(classification=status, primal_distance=distance, solver_objective=float(fit.fun),
        support_indices=indices.tolist(), support_weights=weights[indices].tolist(),
        barycenter=barycenter.tolist(), separating_direction=direction.tolist(),
        direction_l1_norm=norm, separating_lower_bound=lower, certificate_valid=valid)
    return result
