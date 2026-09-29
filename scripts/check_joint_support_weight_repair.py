"""Known support and outside-hull examples for nonnegative witness repair."""
from copy import deepcopy
import numpy as np
from repair_joint_support_weights import repair_weights
from check_joint_support_certificate import check_certificate


def source(x, weights):
    barycenter = np.asarray(weights) @ x
    return dict(records=len(x), scales=[1.], solver_status=0,
                solver_message='Synthetic successful solver certificate',
                classification='unresolved_certificate', support_indices=[0, 1],
                support_weights=weights, barycenter=barycenter.tolist(),
                primal_distance=float(abs(barycenter[0])), solver_objective=0.,
                separating_direction=[0.], direction_l1_norm=0.,
                separating_lower_bound=0., certificate_valid=False)


x = np.array([[1.], [0.]])
original = source(x, [-1e-10, 1+1e-10])
saved = deepcopy(original)
result = repair_weights(x, original)
assert result['accepted'] and result['candidate']['support_indices'] == [1]
assert result['candidate']['primal_distance'] == 0
assert original == saved == result['original']

# Both points are strictly positive: clipping must not invent zero support.
epsilon = 1e-4
x = np.array([[1.], [epsilon]])
original = source(x, [-epsilon/(1-epsilon), 1/(1-epsilon)])
result = repair_weights(x, original)
assert not result['accepted']
assert result['candidate']['classification'] == 'unresolved_certificate'
assert result['candidate']['primal_distance'] == epsilon

# The input witness is checked before any repair is attempted.
bad = deepcopy(original)
bad['barycenter'] = [0.1]
try:
    repair_weights(x, bad)
except AssertionError:
    pass
else:
    raise AssertionError('Corrupt original witness accepted')

# Serialized candidate corruption must be detected by the independent checker.
bad = deepcopy(result['candidate'])
bad['classification'] = 'zero_supported_to_numeric_tolerance'
try:
    check_certificate(x, bad)
except AssertionError:
    pass
else:
    raise AssertionError('False support classification accepted')
print('Nonnegative witness repair passed: exact-support case, outside-hull rejection, original preservation, source and candidate corruption detection.')
