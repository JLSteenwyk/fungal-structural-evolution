"""Exact cone maps with a separate positive residual diagonal and target I.

These algebraic maps do not qualify a residual weighting policy or numerical
identifiability. Uniform qualification cannot be transferred to this model.
"""
import numpy as np

from covariance_exact_folds_v2 import ORIGINAL_NAMES, RELATIONS
from full_expanded_model_design_sources import digest
from reduced_covariance_basis import validate_certificate

SCHEMA = 'exact-positive-diagonal-covariance-cone-v1'


def validate_parent(parent):
    validate_certificate(parent)
    c = parent['certificate']
    assert c['relations'][RELATIONS[1]]['exact'] is True
    assert c['relations'][RELATIONS[2]]['exact'] is True
    return c


def mapping(parent):
    c = validate_parent(parent)
    pair_fold = c['relations'][RELATIONS[0]]['exact']
    names = ['residual', 'target_node', 'background_node']
    if not pair_fold:
        names.append('model_pair')
    names += ['family_intercept', 'species']
    coefficients = {name: np.eye(9)[ORIGINAL_NAMES.index(name)].copy() for name in names}
    coefficients['target_node'][ORIGINAL_NAMES.index('gene')] = .5
    coefficients['background_node'][ORIGINAL_NAMES.index('gene')] = .5
    if pair_fold:
        for name in ['target_node', 'background_node']:
            coefficients[name][ORIGINAL_NAMES.index('model_pair')] = 1.
            coefficients[name][ORIGINAL_NAMES.index('model')] = .5
    else:
        coefficients['model_pair'][ORIGINAL_NAMES.index('model')] = .5
    coefficients['family_intercept'][ORIGINAL_NAMES.index('family')] = (
        0. if c['loading_mode'] == 'signed' else 4.)
    forward = np.asarray([coefficients[name] for name in names])
    inverse = np.zeros((9, len(names)))
    for j, name in enumerate(names):
        inverse[ORIGINAL_NAMES.index(name), j] = 1.
    assert np.array_equal(forward @ inverse, np.eye(len(names)))
    assert np.all(forward >= 0) and np.all(inverse >= 0)
    assert np.array_equal(2 * forward, np.rint(2 * forward))
    return dict(original_names=ORIGINAL_NAMES, retained_names=names,
        forward=forward.tolist(), nonnegative_right_inverse=inverse.tolist(),
        residual_diagonal_scope='any_finite_strictly_positive_diagonal',
        residual_and_target_identity_kept_distinct=True,
        covariance_cone_preserved_given_exact_certificates=True,
        residual_diagonal_prepared=False, raw_reml_basis_qualification_complete=False,
        nonuniform_weighting_accepted=False, component_variance_attribution_accepted=False)


def record(parent, contract):
    c = validate_parent(parent)
    return dict(schema=SCHEMA, source_contract=contract, cohort_id=parent['cohort_id'],
        loading_mode=c['loading_mode'], records=c['records'],
        cohort_rows_sha256=parent['cohort_rows_sha256'],
        ordered_case_ids_sha256=parent['ordered_case_ids_sha256'],
        parent_exact_certificate_sha256=digest(parent), variance_map=mapping(parent),
        parent_pair_identity_exact=c['relations'][RELATIONS[0]]['exact'],
        scientific_eligibility=False)


def validate_diagonal(diagonal, records):
    diagonal = np.asarray(diagonal, dtype=float)
    assert diagonal.shape == (records,) and np.isfinite(diagonal).all()
    assert np.all(diagonal > 0)
    return diagonal
