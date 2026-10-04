"""Fit a source-qualified named covariance basis with its actual residual D.

The whole-grid adapter must require closed source and numerical readbacks.
This bridge verifies a single identity and adds fresh production guards; it
does not confer source closure or inferential acceptance.
"""
import numpy as np

from covariance_basis_context import ComponentKernelProducts
from fit_shared_entity_likelihood import fit_shared_entity
from full_expanded_model_design_sources import array_digest, digest
from full_weighted_covariance_qualification import disposition
from full_weighted_covariance_sources_v2 import POLICIES
from readback_full_covariance_qualification import numeric as audit_qualification
from shared_entity_likelihood import SharedEntityLikelihood

READY = 'numerically_qualified_exact_four_control_covariance_basis'


def validate_source(record, audit, route, diagonal):
    assert audit is not None and route is not None
    d = np.asarray(diagonal, dtype=float)
    assert d.shape == (audit['records'],) and np.isfinite(d).all() and np.all(d > 0)
    assert record['covariance_audit_sha256'] == digest(audit)
    assert record['covariance_audit_id'] == audit['audit_id']
    assert record['source_covariance_disposition'] == audit['disposition']
    assert record['control_policy'] == audit['control_policy'] and record['control_policy'] in POLICIES
    for name in ['loading_mode', 'tree', 'records', 'diagonal_sha256',
                 'control_record_sha256', 'exact_certificate_sha256', 'retained_kernel_names']:
        assert record[name] == audit[name], name
    assert record['method'] in ['ml', 'reml']
    for value in [record, audit]:
        for name in ['scientific_eligibility', 'nonuniform_weighting_accepted',
                     'component_variance_attribution_accepted']:
            assert value[name] is False
    assert audit['covariance_model_selected'] is False
    assert array_digest(d, '<f8') == audit['diagonal_sha256']
    assert route['certificate_sha256'] == audit['exact_certificate_sha256']
    assert route['names'] == audit['retained_kernel_names']
    assert route['route'] == audit['basis_route']
    assert route['exact_uniform_one'] is audit['residual_diagonal_is_exact_uniform_one']
    assert route['exact_uniform_one'] is bool(np.array_equal(d, np.ones(len(d))))
    expected_names = ['residual'] + ([] if route['exact_uniform_one'] else ['target_node']) + ['background_node']
    if 'model_pair' in route['names']: expected_names.append('model_pair')
    expected_names += ['family_intercept', 'species']
    assert route['names'] == expected_names
    assert route['route'] == ('closed_exact_uniform_named_fold' if route['exact_uniform_one'] else 'closed_positive_diagonal_cone')
    if audit['numerical_audit'] is not None:
        assert audit['numerical_audit']['kernel_names'] == route['names']
        assert audit['disposition'] == disposition(audit['numerical_audit'])
    expected = (record['source_fit_disposition'] if record['source_fit_disposition'] != 'ready_for_working_covariance_fit'
                else audit['disposition'])
    assert record['source_combined_disposition'] == expected
    return d


def backend_guard(source, rows, matrix, operators, audit, diagonal):
    """Fresh actual-D guard in the optimizer's own component implementation."""
    assert audit['disposition'] == READY
    names = audit['retained_kernel_names']; assert list(operators) == names[1:-1]
    inherited = audit['numerical_audit']
    assert all(inherited[k]['disposition'] == 'numerically_independent_covariance_bases'
               for k in ['raw_diagnostics', 'reml_diagnostics'])
    fresh = ComponentKernelProducts(source['labels'][rows], operators, diagonal).tree(
        source['factors'][audit['tree']][rows]).audit(matrix)
    for name in ['kernel_names', 'records', 'fixed_effect_columns', 'residual_dimension',
                 'family_components', 'exactly_zero_incidence_names']:
        assert fresh[name] == inherited[name], name
    reference = dict(raw=np.asarray(fresh['raw_gram']), projected=np.asarray(fresh['projected_gram']),
        raw_error=np.asarray(fresh['raw_roundoff_envelope']), projected_error=np.asarray(fresh['projected_roundoff_envelope']),
        normalized_design_condition_number=fresh['normalized_design_condition_number'])
    review, _ = audit_qualification(inherited, reference, names, len(rows))
    assert review is False
    if any(fresh[k]['disposition'] != 'numerically_independent_covariance_bases'
           for k in ['raw_diagnostics', 'reml_diagnostics']):
        raise ValueError('Fresh actual-D backend requires additional covariance review')
    return fresh


def candidate(source, plan, rows, record, matrix, response, operators, audit, route, diagonal):
    d = validate_source(record, audit, route, diagonal)
    assert np.asarray(response).shape == (len(rows),) and np.isfinite(response).all()
    assert record['response_sha256'] == array_digest(response, '<f8')
    if record['source_combined_disposition'] != READY:
        return dict(**record, disposition=record['source_combined_disposition'], fit=None, numerical_attempted=False)
    try:
        guard = backend_guard(source, rows, matrix, operators, audit, d)
        likelihood = SharedEntityLikelihood(source['labels'][rows], operators,
            source['factors'][record['tree']][rows], d, matrix, response)
        fitted = fit_shared_entity(likelihood, record['method'], **plan['optimizer'])
        if 'source_covariance_qualification' in fitted:
            assert fitted['source_covariance_qualification'] == guard
        return dict(**record, disposition=fitted['status'], fit=fitted, numerical_attempted=True)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        return dict(**record, disposition='shared_entity_fit_error_requires_review', fit=None,
            numerical_attempted=True, error_type=type(error).__name__, error_message=str(error))
