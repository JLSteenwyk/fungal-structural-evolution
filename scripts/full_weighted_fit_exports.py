"""Segmented complete-grid weighted fitting exports and source-link contracts."""
from collections import Counter
import csv
import gzip
import itertools
import json
import os
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from full_expanded_model_design_sources import SETTING_FIELDS, digest
from full_weighted_covariance_qualification import LINK_FIELDS as SOURCE_LINK_FIELDS
from full_weighted_shared_entity_fit_sources import SCHEMA
from weighted_shared_entity_candidate import READY

PRODUCER = 'complete_full_four_control_shared_entity_fits_pending_readback_v1'
READER = 'passed_full_four_control_shared_entity_grid_and_spectral_readback_v1'
LINK_FIELDS = SOURCE_LINK_FIELDS + ['likelihood_method', 'candidate_id', 'candidate_disposition']
INDEPENDENT_LINK_FIELDS = LINK_FIELDS + ['independent_disposition']
SUMMARY = ['logical_cases', 'cohorts', 'designs', 'fit_inputs', 'settings', 'candidate_rows',
    'setting_fit_links', 'candidate_status_counts', 'setting_status_counts', 'methods', 'trees', 'loading_modes', 'policies']
INDEPENDENT_STATUSES = {
    'source_review_or_exclusion_retained', 'original_numerical_failure_reproduced_requires_review',
    'producer_candidate_requires_review_independent_numeric_replay_passed',
    'independently_audited_working_candidate_pending_inferential_calibration',
    'independent_curvature_or_search_requires_review', 'independent_numerical_precision_requires_review'}


def atomic(path, value):
    assert not path.exists(), 'Immutable completed artifact already exists: ' + str(path)
    temporary = path.with_name(path.name + '.partial')
    with temporary.open('x') as f:
        json.dump(value, f, sort_keys=True, allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(temporary, path)


def finish_gzip(temporary, path):
    assert not path.exists()
    with temporary.open('rb') as f: os.fsync(f.fileno())
    os.replace(temporary, path)


def validate_export(expected, value):
    assert all(value[k] == v for k, v in expected.items())
    assert all(value[k] is False for k in ['scientific_eligibility', 'nonuniform_weighting_accepted',
                                         'component_variance_attribution_accepted'])
    keys = set(expected) | {'disposition', 'fit', 'numerical_attempted'}
    if value['disposition'] == 'shared_entity_fit_error_requires_review': keys |= {'error_type', 'error_message'}
    assert set(value) == keys
    if expected['source_combined_disposition'] != READY:
        assert value == dict(**expected, disposition=expected['source_combined_disposition'], fit=None, numerical_attempted=False)
    else:
        assert value['numerical_attempted'] is True
        if value['fit'] is None:
            assert value['disposition'] == 'shared_entity_fit_error_requires_review'
            assert value['error_type'] in ['ValueError', 'ArithmeticError', 'LinAlgError']
            assert isinstance(value['error_message'], str)
        else:
            fitted = value['fit']; assert fitted['scientific_eligibility'] is False
            assert fitted['status'] == value['disposition'] and fitted['method'] == expected['method']
            assert fitted['status'] in ['all_shared_entity_optimizer_searches_require_review',
                'optimized_shared_entity_candidate_pending_independent_audit', 'shared_entity_optimizer_candidate_requires_review']
            assert fitted['parameter_names'] == expected['retained_kernel_names'][1:]
    # This structural check does not replace the independent numerical reader.
    json.dumps(value, allow_nan=False)


def key(identity):
    return tuple(identity[k] for k in ['fit_input_id', 'loading_mode', 'tree', 'control_policy', 'method'])


def link_rows(setting_rows, plan, index):
    for ordinal, row in setting_rows:
        assert len(row) == len(SETTING_FIELDS); setting = dict(zip(SETTING_FIELDS, row))
        for mode, tree, policy, method in itertools.product(plan['loading_modes'], plan['trees'], plan['policies'], plan['methods']):
            identity, status = index[setting['fit_input_id'], mode, tree, policy, method][:2]
            assert identity['design_id'] == setting['design_id'] and identity['cohort_id'] == setting['cohort_id']
            assert identity['outcome'] == setting['outcome'] and identity['records'] == int(setting['records'])
            covariance = identity['source_covariance_disposition']
            combined = setting['disposition'] if setting['disposition'] != 'ready_for_working_covariance_fit' else covariance
            assert combined == identity['source_combined_disposition']
            yield [str(ordinal), digest(row), *row, mode, tree, policy, identity['covariance_audit_id'],
                covariance, combined, method, identity['candidate_id'], status]


def process_links(path, expected_rows, reader=False, independent_path=None, index=None):
    counts = Counter(); audits = Counter(); n = 0
    if reader:
        source = gzip.open(path, 'rt'); parsed = csv.reader(source, delimiter='\t')
        assert next(parsed) == LINK_FIELDS
    else:
        temporary = path.with_name(path.name + '.partial')
        source = gzip.open(temporary, 'xt'); parsed = csv.writer(source, delimiter='\t', lineterminator='\n')
        parsed.writerow(LINK_FIELDS)
    target = output = temp = None
    if independent_path is not None:
        existing = independent_path.exists()
        if existing:
            target = gzip.open(independent_path, 'rt'); output = csv.reader(target, delimiter='\t')
            assert next(output) == INDEPENDENT_LINK_FIELDS
        else:
            temp = independent_path.with_name(independent_path.name + '.partial')
            target = gzip.open(temp, 'xt'); output = csv.writer(target, delimiter='\t', lineterminator='\n')
            output.writerow(INDEPENDENT_LINK_FIELDS)
    try:
        for expected in expected_rows:
            if reader: assert next(parsed) == expected
            else: parsed.writerow(expected)
            if target is not None:
                candidate = index[expected[LINK_FIELDS.index('fit_input_id')],
                    expected[LINK_FIELDS.index('loading_mode')], expected[LINK_FIELDS.index('tree')],
                    expected[LINK_FIELDS.index('control_policy')], expected[LINK_FIELDS.index('likelihood_method')]]
                independent = candidate[2]; assert independent in INDEPENDENT_STATUSES
                row = expected + [independent]
                if existing: assert next(output) == row
                else: output.writerow(row)
                audits[independent] += 1
            counts[expected[-1]] += 1; n += 1
        if reader: assert next(parsed, None) is None
        if target is not None and existing: assert next(output, None) is None
    finally:
        source.close()
        if target is not None: target.close()
    if not reader: finish_gzip(temporary, path)
    if temp is not None: finish_gzip(temp, independent_path)
    return n, counts, audits


def failure_capture(directory, identity, value, rows, diagonal, matrix, response, factor, error):
    path = directory / identity['candidate_id']; path.mkdir(parents=True, exist_ok=False)
    atomic(path / 'failure.json', dict(expected_identity=identity, saved_candidate_representation=repr(value),
        error_type=type(error).__name__, error_message=str(error), scientific_eligibility=False))
    with (path / 'original_numeric_inputs.npz').open('xb') as f:
        np.savez_compressed(f, case_rows=rows, diagonal=diagonal, design=matrix, response=response, species_factor=factor)


def summary(source, plan, count, links, candidates, settings):
    q = source['numerical_completion']
    assert count == q['designs'] * 2 * 2 * 2 * 5 * 4 == plan['expected']['candidate_rows']
    assert links == q['settings'] * 2 * 2 * 5 * 4 == plan['expected']['setting_fit_links']
    return dict(logical_cases=len(source['ids']), cohorts=q['cohorts'], designs=q['designs'], fit_inputs=q['designs'] * 2,
        settings=q['settings'], candidate_rows=count, setting_fit_links=links, candidate_status_counts=dict(candidates),
        setting_status_counts=dict(settings), methods=plan['methods'], trees=plan['trees'],
        loading_modes=plan['loading_modes'], policies=plan['policies'])
