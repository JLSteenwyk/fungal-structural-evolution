"""Construct a separate convex-hull witness; never change a saved certificate."""
from copy import deepcopy
import math
import numpy as np
from check_joint_support_certificate import check_certificate


def repair_weights(x, original):
    """Clip negative weights and renormalize, then check the actual barycenter.

    This is a new feasibility witness, not a new optimizer result. The original
    solver objective and dual certificate remain unchanged and must agree at
    the existing tolerances. Callers must bind x to audited source identities.
    An unsuccessful repair remains explicitly unresolved.
    """
    check_certificate(x, original)
    if original['solver_status'] != 0 or original['classification'] != 'unresolved_certificate':
        raise ValueError('Expected an unresolved certificate from a successful solver')
    old_weights = np.asarray(original['support_weights'], dtype=float)
    if not np.any(old_weights < 0):
        raise ValueError('No negative weights to repair')
    clipped = np.maximum(old_weights, 0.)
    mass = math.fsum(clipped)
    if not math.isfinite(mass) or mass <= 0:
        raise ValueError('No positive weight mass')
    weights = clipped / mass
    indices = np.asarray(original['support_indices'], dtype=int)
    keep = weights > 0
    weights, indices = weights[keep], indices[keep]
    z = np.asarray(x, dtype=float) / np.asarray(original['scales'])
    barycenter = weights @ z[indices]
    distance = float(np.max(np.abs(barycenter)))
    candidate = deepcopy(original)
    valid = bool(np.all(weights >= 0) and abs(math.fsum(weights)-1) <= 1e-8
                 and abs(distance-original['solver_objective']) <= 1e-8
                 and original['direction_l1_norm'] <= 1+1e-8
                 and original['separating_lower_bound'] <= distance+1e-8)
    classification = 'unresolved_certificate'
    if valid:
        if distance <= 1e-8:
            classification = 'zero_supported_to_numeric_tolerance'
        elif original['separating_lower_bound'] > 1e-7:
            classification = 'zero_outside_joint_convex_hull'
        else:
            classification = 'unresolved_near_boundary'
    candidate.update(support_indices=indices.tolist(), support_weights=weights.tolist(),
                     barycenter=barycenter.tolist(), primal_distance=distance,
                     certificate_valid=valid, classification=classification,
                     witness_method='clip_negative_weights_then_normalize_positive_mass')
    # Independent fsum reconstruction; this does not call an optimizer.
    check_certificate(x, candidate)
    return dict(original=deepcopy(original), candidate=candidate,
                accepted=classification == 'zero_supported_to_numeric_tolerance',
                removed_negative_mass=math.fsum(-v for v in old_weights if v < 0),
                scope='Numerical zero-reference hull witness only; original solver fields retained. No interior-overlap, model-adequacy or biological inference.')
