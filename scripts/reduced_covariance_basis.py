"""Qualify exact retained kernels using principal submatrices of closed Grams.

The original numerical envelopes are copied, never narrowed. Exact named
operator identities, not numerical ranks, determine the retained kernels.
"""
from copy import deepcopy
import numpy as np
from scipy import linalg

from covariance_basis_audit import _diagnostics
from covariance_exact_folds_v2 import variance_map, independent_variance_map
from full_expanded_model_design_sources import digest

SCHEMA = 'full-exact-retained-uniform-covariance-qualification-v1'
OLD_NAMES = ['residual', 'background_node', 'model_pair', 'gene', 'model', 'family_intercept', 'species']
MATRIX_FIELDS = ['raw_gram', 'projected_gram', 'raw_roundoff_envelope', 'projected_roundoff_envelope']


def validate_certificate(row):
    c = row['certificate']; mode = row['certificate']['loading_mode']
    assert row['scientific_eligibility'] is row['raw_reml_basis_qualification_complete'] is False
    assert c['target_identity_exact'] is c['family_fold_exact'] is True
    assert type(c['records']) is int and c['records'] > 0
    expected = variance_map(c['relations'], mode)
    assert expected == independent_variance_map(c['relations'], mode) == row['variance_map']
    names = expected['retained_names']
    assert names == [k for k in OLD_NAMES if k in names]
    assert names[0] == 'residual' and names[-2:] == ['family_intercept', 'species']
    for relation in c['relations'].values():
        assert type(relation['exact']) is bool
        assert type(relation['nonzero_entries']) is int and relation['nonzero_entries'] >= 0
        assert relation['exact'] == (relation['nonzero_entries'] == 0) == (relation['witness'] is None)
    return names


def disposition(value):
    return ('qualified_exact_retained_uniform_covariance_basis' if all(
        value[k]['disposition'] == 'numerically_independent_covariance_bases'
        for k in ['raw_diagnostics', 'reml_diagnostics']) else 'exact_retained_covariance_basis_requires_review')


def identity(contract, old_id):
    return digest([SCHEMA, contract, old_id])


def reduce_record(original, certificate, contract):
    names = validate_certificate(certificate)
    assert original['cohort_id'] == certificate['cohort_id']
    assert original['loading_mode'] == certificate['certificate']['loading_mode']
    assert original['records'] == certificate['certificate']['records']
    assert original['residual_diagonal'] == 'uniform_one'
    result = deepcopy(original)
    result.update(audit_id=identity(contract, original['audit_id']), source_audit_id=original['audit_id'],
        source_covariance_disposition=original['disposition'], retained_source_contract=contract,
        exact_certificate_sha256=digest(certificate), retained_kernel_names=names,
        removed_kernel_names=[k for k in OLD_NAMES if k not in names],
        scientific_eligibility=False, component_variance_attribution_accepted=False,
        covariance_cone_preserved_given_exact_certificates=True, nonuniform_weighting_accepted=False)
    if original['numerical_audit'] is None:
        assert original['disposition'] != 'qualified_uniform_working_covariance_basis'
        return result
    value = original['numerical_audit']
    assert original['source_design_disposition'] == 'full_rank_design'
    assert value['kernel_names'] == OLD_NAMES and value['records'] == original['records']
    assert value['covariance_model_selected'] is False
    selected = [OLD_NAMES.index(k) for k in names]
    reduced = deepcopy(value); reduced['kernel_names'] = names
    for field in MATRIX_FIELDS:
        matrix = np.asarray(value[field], dtype=float)
        assert matrix.shape == (7, 7) and np.isfinite(matrix).all()
        assert np.array_equal(matrix, matrix.T)
        if field.endswith('envelope'): assert np.all(matrix >= 0)
        reduced[field] = matrix[np.ix_(selected, selected)].tolist()
    reduced['exactly_zero_incidence_names'] = [k for k in value['exactly_zero_incidence_names'] if k in names]
    reduced['raw_diagnostics'] = _diagnostics(np.asarray(reduced['raw_gram']), np.asarray(reduced['raw_roundoff_envelope']), names)
    reduced['reml_diagnostics'] = _diagnostics(np.asarray(reduced['projected_gram']), np.asarray(reduced['projected_roundoff_envelope']), names)
    reduced['scope'] = ('Principal submatrices of independently closed full seven-kernel raw/REML Grams; '
        'original numerical error envelopes copied exactly. Retention authorized by closed named integer '
        'operator certificates and nonnegative cone maps. Numerical identifiability only; no variance fit, '
        'separate component attribution, nonuniform weighting, biological acceptance or calibration.')
    result['numerical_audit'] = reduced
    result['disposition'] = disposition(reduced)
    return result


def independent_diagnostics(gram, error, names):
    """Independent SciPy gesvd classification with the same review policy."""
    gram = np.asarray(gram, dtype=float); error = np.asarray(error, dtype=float)
    diagonal = np.diag(gram)
    assert np.isfinite(gram).all() and np.isfinite(error).all() and np.all(error >= 0)
    if np.any(diagonal < -np.diag(error)): raise ArithmeticError('Negative Gram norm beyond envelope')
    active = [i for i in range(len(names)) if diagonal[i] > error[i, i]]
    unresolved = [i for i in range(len(names)) if i not in active]
    if active:
        norm = np.sqrt(np.outer(diagonal[active], diagonal[active]))
        normalized = np.asarray([[gram[i, j] for j in active] for i in active]) / norm
        bound = np.asarray([[error[i, j] for j in active] for i in active]) / norm
        singular = linalg.svd(normalized, compute_uv=False, lapack_driver='gesvd')
        roundoff = float(linalg.svd(bound, compute_uv=False, lapack_driver='gesvd')[0])
        if linalg.eigvalsh((normalized + normalized.T) / 2)[0] < -roundoff:
            raise ArithmeticError('Non-PSD Gram beyond envelope')
        tolerance = max(len(active) * np.finfo(float).eps * singular[0], roundoff)
        ranks = [sum(float(s) > tolerance * multiplier for s in singular) for multiplier in [.1, 1., 10.]]
    else:
        singular = np.empty(0); roundoff = tolerance = 0.; ranks = [0, 0, 0]
    status = ('unresolved_kernel_norm_requires_review' if unresolved else
        'rank_boundary_requires_review' if ranks[0] != ranks[2] else
        'dependent_covariance_bases_require_review' if ranks[1] < len(active) else
        'numerically_independent_covariance_bases')
    return dict(disposition=status, resolved_kernel_indices=active,
        unresolved_kernel_names=[names[i] for i in unresolved], normalized_singular_values=singular.tolist(),
        propagated_normalized_roundoff_bound=roundoff, rank_tolerance=float(tolerance),
        rank_at_one_tenth_tolerance=ranks[0], rank=ranks[1], rank_at_ten_times_tolerance=ranks[2])


def readback_record(saved, original, certificate, contract):
    """Independently select every matrix entry and recompute all classifications."""
    names = validate_certificate(certificate)
    assert saved['source_audit_id'] == original['audit_id']
    assert saved['audit_id'] == digest([SCHEMA, contract, original['audit_id']])
    assert saved['retained_source_contract'] == contract
    assert saved['source_covariance_disposition'] == original['disposition']
    assert saved['exact_certificate_sha256'] == digest(certificate)
    assert saved['retained_kernel_names'] == names
    assert saved['removed_kernel_names'] == [k for k in OLD_NAMES if k not in names]
    assert saved['scientific_eligibility'] is saved['component_variance_attribution_accepted'] is saved['nonuniform_weighting_accepted'] is False
    assert saved['covariance_cone_preserved_given_exact_certificates'] is True
    assert original['cohort_id'] == certificate['cohort_id'] == saved['cohort_id']
    assert original['loading_mode'] == certificate['certificate']['loading_mode'] == saved['loading_mode']
    assert original['records'] == certificate['certificate']['records'] == saved['records']
    changed = {'audit_id', 'disposition', 'numerical_audit'}
    for key in set(original) - changed: assert saved[key] == original[key], key
    if original['numerical_audit'] is None:
        assert saved['numerical_audit'] is None and saved['disposition'] == original['disposition']
        return
    old = original['numerical_audit']; value = saved['numerical_audit']
    assert old['kernel_names'] == OLD_NAMES and value['kernel_names'] == names
    assert old['covariance_model_selected'] is value['covariance_model_selected'] is False
    indices = [old['kernel_names'].index(k) for k in names]
    for key in set(old) - set(MATRIX_FIELDS) - {'kernel_names', 'raw_diagnostics', 'reml_diagnostics', 'exactly_zero_incidence_names', 'scope'}:
        assert value[key] == old[key], key
    assert value['exactly_zero_incidence_names'] == [k for k in old['exactly_zero_incidence_names'] if k in names]
    for key in MATRIX_FIELDS:
        expected = [[old[key][i][j] for j in indices] for i in indices]
        assert value[key] == expected, key
    for field, gram, error in [('raw_diagnostics', 'raw_gram', 'raw_roundoff_envelope'),
                               ('reml_diagnostics', 'projected_gram', 'projected_roundoff_envelope')]:
        reference = independent_diagnostics(value[gram], value[error], names)
        for key in reference:
            if key in ['normalized_singular_values', 'propagated_normalized_roundoff_bound', 'rank_tolerance']:
                np.testing.assert_allclose(value[field][key], reference[key], rtol=1e-12, atol=32*np.finfo(float).eps)
            else: assert value[field][key] == reference[key], (field, key)
    assert saved['disposition'] == disposition(value)
