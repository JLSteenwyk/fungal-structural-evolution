"""Fit exact retained kernels while preserving the closed inherited guard.

Whole-grid consumers must independently require closed original, exact and
retained-basis source proofs. This numerical bridge does not provide closure.
"""
import numpy as np

from covariance_basis_context import ComponentKernelProducts
from fit_shared_entity_likelihood import fit_shared_entity
from full_expanded_model_design_sources import digest
from reduced_covariance_basis import readback_record
from shared_entity_likelihood import SharedEntityLikelihood

READY = 'qualified_exact_retained_uniform_covariance_basis'


def validate_source(record, source_audit, original_audit, certificate):
    assert all(value is not None for value in [source_audit, original_audit, certificate])
    assert record['covariance_audit_sha256'] == digest(source_audit)
    assert record['covariance_audit_id'] == source_audit['audit_id']
    assert record['source_covariance_disposition'] == source_audit['disposition']
    assert record['records'] == source_audit['records']
    assert record['loading_mode'] == source_audit['loading_mode']
    assert record['tree'] == source_audit['tree']
    assert record['scientific_eligibility'] is False
    readback_record(source_audit, original_audit, certificate, source_audit['retained_source_contract'])
    names = source_audit['retained_kernel_names']
    assert names in [['residual', 'background_node', 'family_intercept', 'species'],
                     ['residual', 'background_node', 'model_pair', 'family_intercept', 'species']]
    return names


def backend_guard(source, rows, matrix, operators, source_audit):
    """Additional fresh guard; it never replaces inherited error envelopes."""
    names = source_audit['retained_kernel_names']
    assert list(operators) == names[1:-1]
    assert source_audit['disposition'] == READY
    inherited = source_audit['numerical_audit']
    assert inherited['kernel_names'] == names
    assert all(inherited[k]['disposition'] == 'numerically_independent_covariance_bases'
               for k in ['raw_diagnostics', 'reml_diagnostics'])
    fresh = ComponentKernelProducts(source['labels'][rows], operators, np.ones(len(rows))).tree(
        source['factors'][source_audit['tree']][rows]).audit(matrix)
    assert fresh['kernel_names'] == names
    for key in ['records', 'fixed_effect_columns', 'residual_dimension', 'family_components',
                'exactly_zero_incidence_names']:
        assert fresh[key] == inherited[key], key
    condition = fresh['normalized_design_condition_number']
    assert abs(1 / condition - 1 / inherited['normalized_design_condition_number']) <= (
        4 * np.finfo(float).eps * max(len(rows), matrix.shape[1]))
    for key, error in [('raw_gram', 'raw_roundoff_envelope'),
                       ('projected_gram', 'projected_roundoff_envelope')]:
        a = np.asarray(fresh[key]); b = np.asarray(inherited[key])
        np.testing.assert_allclose(a, b, rtol=3e-9, atol=2e-8)
        # Both previously specified contraction bounds must cover differences.
        # Neither bound is overwritten, shrunk, inferred by subtraction or used
        # to promote a source review. The fresh bound can differ with fewer
        # incidence columns; both numerical independence checks must pass.
        total = np.asarray(fresh[error]) + np.asarray(inherited[error])
        rounding = 8 * np.finfo(float).eps * np.maximum(abs(a), abs(b))
        if np.any(abs(a - b) > total + rounding):
            raise ArithmeticError('Retained Gram disagreement exceeds preserved contraction bounds')
    if any(fresh[k]['disposition'] != 'numerically_independent_covariance_bases'
           for k in ['raw_diagnostics', 'reml_diagnostics']):
        raise ValueError('Fresh retained backend requires identifiability review in addition to inherited guard')
    return fresh


def candidate(source, plan, rows, record, matrix, response, operators,
              source_audit, original_audit, certificate):
    validate_source(record, source_audit, original_audit, certificate)
    if record['source_combined_disposition'] != READY:
        assert record['source_combined_disposition'] == (
            record['source_fit_disposition'] if record['source_fit_disposition'] != 'ready_for_working_covariance_fit'
            else source_audit['disposition'])
        return dict(**record, disposition=record['source_combined_disposition'], fit=None, numerical_attempted=False)
    assert record['source_fit_disposition'] == 'ready_for_working_covariance_fit'
    try:
        guard = backend_guard(source, rows, matrix, operators, source_audit)
        likelihood = SharedEntityLikelihood(source['labels'][rows], operators,
            source['factors'][record['tree']][rows], np.ones(len(rows)), matrix, response)
        fitted = fit_shared_entity(likelihood, record['method'], **plan['optimizer'])
        if 'source_covariance_qualification' in fitted:
            assert fitted['source_covariance_qualification'] == guard
        return dict(**record, disposition=fitted['status'], fit=fitted, numerical_attempted=True)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        return dict(**record, disposition='shared_entity_fit_error_requires_review', fit=None,
                    numerical_attempted=True, error_type=type(error).__name__, error_message=str(error))
