"""Exact four-control probe exports, numeric replay and review-input custody."""
import hashlib
import json
from pathlib import Path
import numpy as np

from full_expanded_model_design_sources import array_digest, digest
from full_weighted_fit_exports import atomic
from full_weighted_shared_entity_fit_sources import operators_for
from full_weighted_shared_entity_timing import probe

COMPARISONS = ['coefficient_and_conditional_covariance_comparison_passed',
    'profiled_scale_comparison_passed', 'variance_components_comparison_passed']
SETUP_TIMES = ['source_validation_seconds', 'primary_constructor_seconds',
    'production_qualification_guard_seconds', 'reader_qualification_seconds', 'independent_constructor_seconds']
BASE = {'group_id', 'representative', 'eligible_candidates', 'selection_active_columns',
    'selection_condition', 'variance_points', 'scientific_eligibility', 'status'}
POINT = {'scaled_variance', 'status', 'primary_evaluation_seconds', 'independent_evaluation_seconds',
    'objective_absolute_error', 'maximum_coordinate_gradient_error', *COMPARISONS,
    'independent_inverse_relative_residual', 'explicit_whitening_fallback_batches', 'independent_cached_kernel_bytes'}
ERRORS = ['ValueError', 'ArithmeticError', 'LinAlgError']


def finite(value, positive=False):
    assert type(value) in [int, float] and np.isfinite(value) and (value > 0 if positive else value >= 0)


def validate_probes(records, selected, counts, plan, fit):
    assert [r['group_id'] for r in records] == sorted(selected)
    for record in records:
        expected = selected[record['group_id']]
        assert record['representative'] == expected['identity']
        assert record['eligible_candidates'] == counts[record['group_id']]
        assert record['selection_active_columns'] == expected['matrix'].shape[1]
        assert record['selection_condition'] == expected['rank'][1]
        assert record['variance_points'] == plan['scaled_variance_points'] and record['scientific_eligibility'] is False
        for name, dtype, values in [('response', '<f8', expected['response']), ('diagonal', '<f8', expected['diagonal'])]:
            assert array_digest(values, dtype) == expected['identity'][name + '_sha256']
        if record['status'] == 'timing_group_construction_requires_review':
            assert set(record) <= BASE | set(SETUP_TIMES) | {'error_type', 'error_message'}
            assert BASE | {'error_type', 'error_message'} <= set(record)
            assert record['error_type'] in ERRORS and isinstance(record['error_message'], str)
            for name in set(record) & set(SETUP_TIMES): finite(record[name])
            continue
        assert record['status'] in ['timed_all_declared_points_agree_only', 'timing_group_requires_review']
        assert set(record) == BASE | set(SETUP_TIMES) | {'parameter_names', 'kernel_normalization', 'points'}
        for name in SETUP_TIMES: finite(record[name])
        assert record['parameter_names'] == expected['source_audit']['retained_kernel_names'][1:]
        assert len(record['kernel_normalization']) == len(record['parameter_names'])
        for value in record['kernel_normalization']: finite(value, positive=True)
        assert [p['scaled_variance'] for p in record['points']] == plan['scaled_variance_points']
        for p in record['points']:
            if p['status'] == 'timing_probe_precision_requires_review':
                assert set(p) == {'scaled_variance', 'status', 'error_type', 'error_message', 'elapsed_seconds'}
                assert p['error_type'] in ERRORS and isinstance(p['error_message'], str); finite(p['elapsed_seconds'])
                continue
            assert p['status'] in ['timed_likelihood_agreement_only', 'timing_probe_numerical_agreement_requires_review']
            assert set(p) == POINT
            for name in ['primary_evaluation_seconds', 'independent_evaluation_seconds',
                         'objective_absolute_error', 'maximum_coordinate_gradient_error', 'independent_inverse_relative_residual']:
                finite(p[name])
            for name in ['explicit_whitening_fallback_batches', 'independent_cached_kernel_bytes']:
                assert type(p[name]) is int and p[name] >= 0
            assert all(type(p[k]) is bool for k in COMPARISONS)
            agreement = (p['objective_absolute_error'] <= fit['independent_audit']['replay']['objective_atol'] and
                p['maximum_coordinate_gradient_error'] <= fit['independent_audit']['replay']['gradient_atol'] and all(p[k] for k in COMPARISONS))
            assert (p['status'] == 'timed_likelihood_agreement_only') is agreement
        assert (record['status'] == 'timed_all_declared_points_agree_only') == all(
            p['status'] == 'timed_likelihood_agreement_only' for p in record['points'])


def semantic(value):
    """Hardware durations are observations, not reproducible numeric results."""
    if isinstance(value, dict): return {k:semantic(v) for k,v in value.items() if not k.endswith('_seconds')}
    if isinstance(value, list): return [semantic(v) for v in value]
    return value


def same_numeric(expected, actual):
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and set(expected) == set(actual)
        for k in expected: same_numeric(expected[k], actual[k])
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(expected) == len(actual)
        for a,b in zip(expected, actual): same_numeric(a,b)
    elif type(expected) is float:
        # Check serialization/replay reproducibility without altering any
        # qualification or primary-versus-independent acceptance tolerance.
        assert type(actual) in [int,float] and np.isclose(expected, actual, rtol=1e-10, atol=1e-12)
    else: assert type(expected) is type(actual) and expected == actual


def replay(source, fit, plan, rows, records, selected, counts):
    operators = {}
    for record in records:
        r = selected[record['group_id']]
        key = (r['identity']['loading_mode'], tuple(r['source_audit']['retained_kernel_names']))
        if key not in operators: operators[key] = operators_for(source, rows, r['source_audit'])
        repeated = probe(source, fit, plan, rows, operators[key], r, counts[record['group_id']], record['group_id'])
        same_numeric(semantic(record), semantic(repeated))


def review_inputs(root, record, representative, source, rows, read=False):
    """Retain exact numeric inputs for every construction/precision review."""
    if record['status'] == 'timed_all_declared_points_agree_only': return []
    folder = Path(root) / 'reviews' / record['group_id']
    arrays = dict(case_rows=rows, labels=source['labels'][rows], diagonal=representative['diagonal'],
        design=representative['matrix'], response=representative['response'],
        species_factor=source['factors'][representative['identity']['tree']][rows])
    paths=[];references={};bank=Path(root)/'review_input_arrays'
    for key,value in arrays.items():
        assert value.dtype.kind != 'O'
        specification=dict(dtype=value.dtype.str,shape=list(value.shape),
            original_bytes_sha256=hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest())
        npz=bank/(digest(specification)+'.npz')
        if not read and not npz.exists():
            bank.mkdir(exist_ok=True)
            with npz.open('xb') as f: np.savez_compressed(f, original_array=value)
        with np.load(npz,allow_pickle=False) as saved:
            assert saved.files == ['original_array']
            restored=saved['original_array']
            assert restored.dtype == value.dtype and restored.shape == value.shape
            assert np.ascontiguousarray(restored).tobytes() == np.ascontiguousarray(value).tobytes()
        from ancestral_chain_attempt import sha
        references[key]=dict(**specification,path=str(npz.relative_to(root)),sha256=sha(npz));paths.append(npz)
    description=dict(probe=record,source_audit=representative['source_audit'],route=representative['route'],
        original_numeric_inputs=references,scientific_eligibility=False)
    jp=folder/'review.json'
    if read:
        assert json.loads(jp.read_text()) == description
    else:
        folder.mkdir(parents=True, exist_ok=False); atomic(jp, description)
    return [jp,*paths]
