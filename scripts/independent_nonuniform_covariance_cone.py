"""Independently reconstruct each original kernel's positive-diagonal image."""
from fractions import Fraction
import numpy as np

from covariance_exact_folds_v2 import ORIGINAL_NAMES, RELATIONS
from full_expanded_model_design_sources import digest
from reduced_covariance_basis import validate_certificate


def mapping(parent):
    validate_certificate(parent)
    c = parent['certificate']; relations = c['relations']
    assert relations[RELATIONS[1]]['exact'] is relations[RELATIONS[2]]['exact'] is True
    pair = {'target_node': Fraction(1), 'background_node': Fraction(1)} if relations[RELATIONS[0]]['exact'] else {'model_pair': Fraction(1)}
    images = {
        'residual': {'residual': Fraction(1)},
        'target_node': {'target_node': Fraction(1)},
        'background_node': {'background_node': Fraction(1)},
        'model_pair': pair,
        'gene': {'target_node': Fraction(1, 2), 'background_node': Fraction(1, 2)},
        'model': {name: amount / 2 for name, amount in pair.items()},
        'family': {} if c['loading_mode'] == 'signed' else {'family_intercept': Fraction(4)},
        'family_intercept': {'family_intercept': Fraction(1)},
        'species': {'species': Fraction(1)},
    }
    used = set().union(*(set(image) for image in images.values()))
    names = [name for name in ORIGINAL_NAMES if name in used]
    forward = [[images[old].get(name, Fraction(0)) for old in ORIGINAL_NAMES] for name in names]
    inverse = [[Fraction(int(old == name)) for name in names] for old in ORIGINAL_NAMES]
    for i in range(len(names)):
        for j in range(len(names)):
            assert sum(forward[i][k] * inverse[k][j] for k in range(9)) == Fraction(int(i == j))
    assert all(v >= 0 for row in forward + inverse for v in row)
    return dict(original_names=ORIGINAL_NAMES, retained_names=names,
        forward=[[float(v) for v in row] for row in forward],
        nonnegative_right_inverse=[[float(v) for v in row] for row in inverse],
        residual_diagonal_scope='any_finite_strictly_positive_diagonal',
        residual_and_target_identity_kept_distinct=True,
        covariance_cone_preserved_given_exact_certificates=True,
        residual_diagonal_prepared=False, raw_reml_basis_qualification_complete=False,
        nonuniform_weighting_accepted=False, component_variance_attribution_accepted=False)


def readback(saved, parent, contract):
    c = parent['certificate']
    assert saved == dict(schema='exact-positive-diagonal-covariance-cone-v1',
        source_contract=contract, cohort_id=parent['cohort_id'], loading_mode=c['loading_mode'],
        records=c['records'], cohort_rows_sha256=parent['cohort_rows_sha256'],
        ordered_case_ids_sha256=parent['ordered_case_ids_sha256'],
        parent_exact_certificate_sha256=digest(parent), variance_map=mapping(parent),
        parent_pair_identity_exact=c['relations'][RELATIONS[0]]['exact'], scientific_eligibility=False)
    f = np.asarray(saved['variance_map']['forward'])
    assert f[0].tolist() == [1.] + [0.] * 8
    assert f[1, 1] == 1. and f[1, 0] == 0.
