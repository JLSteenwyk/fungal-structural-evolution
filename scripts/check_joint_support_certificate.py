"""Validate saved support/separation certificates without calling an optimizer."""
import math
import numpy as np


def check_certificate(x, result):
    x = np.asarray(x, dtype=float)
    n, d = x.shape
    assert n == result['records'] and np.isfinite(x).all()
    scale = np.array([max(abs(float(v)) for v in x[:, j]) or 1. for j in range(d)])
    np.testing.assert_array_equal(scale, result['scales'])
    if result['solver_status'] != 0:
        assert result['classification'] == 'unresolved_solver' and result['solver_message']
        return
    indices = result['support_indices']
    weights = result['support_weights']
    assert len(indices) == len(weights) and len(set(indices)) == len(indices)
    assert all(type(i) is int and 0 <= i < n for i in indices)
    assert all(math.isfinite(w) for w in weights)
    z = x / scale
    barycenter = [math.fsum(w * float(z[i, j]) for i, w in zip(indices, weights)) for j in range(d)]
    np.testing.assert_allclose(barycenter, result['barycenter'], rtol=1e-9, atol=1e-12)
    distance = max(map(abs, barycenter))
    direction = np.asarray(result['separating_direction'])
    assert direction.shape == (d,) and np.isfinite(direction).all()
    norm = math.fsum(map(abs, direction))
    # Explicit column accumulation is separate from the solver's matrix product.
    projections = sum((z[:, j] * direction[j] for j in range(d)), np.zeros(n))
    lower = float(min(projections))
    for actual, expected in [(distance,result['primal_distance']), (norm,result['direction_l1_norm']),
                             (lower,result['separating_lower_bound'])]:
        assert math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-12)
    assert math.isfinite(result['solver_objective'])
    valid = (all(w >= 0 for w in weights) and abs(math.fsum(weights)-1) <= 1e-8
             and abs(distance-result['solver_objective']) <= 1e-8
             and norm <= 1+1e-8 and lower <= distance+1e-8)
    assert valid == result['certificate_valid']
    if not valid:
        expected = 'unresolved_certificate'
    elif distance <= 1e-8:
        expected = 'zero_supported_to_numeric_tolerance'
    elif lower > 1e-7:
        expected = 'zero_outside_joint_convex_hull'
    else:
        expected = 'unresolved_near_boundary'
    assert result['classification'] == expected
