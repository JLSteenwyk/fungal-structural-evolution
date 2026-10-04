#!/usr/bin/env python3
"""Reconstruct all original X/y inputs and census the full four-control grid.

This is a source and mapping check. It computes no covariance Grams or fits.
The reader reconstructs every exported cohort record and all original links.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import itertools
import json
from pathlib import Path
import sqlite3
import tempfile

import numpy as np

from ancestral_chain_attempt import sha
from full_expanded_model_design_sources import AXES, DEGREES, OUTCOMES, SETTING_FIELDS, array_digest, fit_id
from full_expanded_model_input_sources import ORDERS
from full_exact_covariance_sources import runtime_caps
from full_weighted_covariance_sources_v2 import MODES, POLICIES, SCHEMA, basis, cohort_controls, design_matrix, jsonl, load
from reference_measurement_union_sources import bind, verify

PRODUCER = 'complete_full_four_control_covariance_source_census_pending_readback_v1'
READER = 'passed_full_four_control_covariance_source_census_readback_v1'
SUMMARY = ['logical_cases', 'cohorts', 'case_row_occurrences', 'designs', 'fit_inputs', 'settings',
    'design_status_counts', 'fit_input_status_counts', 'setting_status_counts', 'policy_basis_counts',
    'expected_numerical_audit_rows', 'expected_setting_audit_links', 'policies', 'trees', 'loading_modes']


def write(path, value):
    with path.open('x') as f:
        json.dump(value, f, sort_keys=True, allow_nan=False); f.write('\n')


def run(path, reader=False):
    runtime_caps()
    plan = json.loads(path.read_text()); source, bindings = load(plan, path)
    original = dict(bindings); root = Path(plan['output'])
    name = 'readback.json' if reader else 'receipt.json'; assert not (root / name).exists()
    if reader:
        producer = json.loads((root / 'receipt.json').read_text())
        assert producer['status'] == PRODUCER and producer['plan_sha256'] == sha(path)
        assert producer['source_contract'] == source['contract'] and producer['source_hashes'] == original
        assert producer['scientific_eligibility'] is producer['raw_reml_basis_qualification_complete'] is False
        assert set(producer['artifacts']) == {'stage_plan.json', 'cohort_source_census.jsonl.gz'}
        for n, h in producer['artifacts'].items(): bind(bindings, root / n, h)
        assert json.loads((root / 'stage_plan.json').read_text()) == dict(schema=SCHEMA, plan_sha256=sha(path), source_contract=source['contract'])
        exported = jsonl(root / 'cohort_source_census.jsonl.gz')
    else:
        root.mkdir(parents=True, exist_ok=False)
        write(root / 'stage_plan.json', dict(schema=SCHEMA, plan_sha256=sha(path), source_contract=source['contract']))
    designs = jsonl(source['root'] / 'unique_designs.jsonl')
    fits = jsonl(source['root'] / 'unique_fit_inputs.jsonl')
    design_counts = Counter(); fit_counts = Counter(); setting_counts = Counter(); basis_counts = Counter()
    design_n = fit_n = settings_n = rows_n = 0
    with tempfile.TemporaryDirectory(prefix='original-weighted-source-', dir=root) as directory:
        db = sqlite3.connect(str(Path(directory) / 'identities.sqlite'))
        db.execute('CREATE TABLE designs(did TEXT PRIMARY KEY,cid TEXT,guide TEXT,mask TEXT,ord TEXT,axis TEXT,degree INTEGER,n INTEGER,status TEXT)')
        db.execute('CREATE TABLE fits(fid TEXT PRIMARY KEY,did TEXT,outcome TEXT,status TEXT,UNIQUE(did,outcome))')
        out = None if reader else gzip.open(root / 'cohort_source_census.jsonl.gz', 'wt')
        try:
            for i, (c, entry) in enumerate(zip(source['cohorts'], source['controls'])):
                selected, diagonals, control = cohort_controls(source, c, entry)
                rows_n += len(selected); choices = []
                for j, policy in enumerate(POLICIES):
                    for mode in MODES:
                        route = basis(source, c, mode, diagonals[j])
                        assert route['exact_uniform_one'] == control['policy_records'][j]['diagonal_is_exact_uniform_one']
                        choices.append(dict(policy=policy, loading_mode=mode, **route,
                            diagonal_sha256=array_digest(diagonals[j], '<f8')))
                        basis_counts[policy + ':' + mode + ':q' + str(len(route['names']))] += 1
                design_ids = []; fit_ids = []; seen = set()
                for _ in range(30):
                    d = next(designs); x = design_matrix(source, c, selected, d)
                    values = source['arrays'][d['mask'], d['order_contrast']]
                    columns = d['predictor_columns']
                    active = [j for j, k in enumerate(columns) if k == 'intercept' or np.any(values[k][selected] != 0)]
                    assert d['active_column_indices'] == active and d['active_columns'] == len(active)
                    assert d['exactly_zero_column_indices'] == [j for j in range(len(columns)) if j not in active]
                    assert x.shape == (len(selected), len(active)) and np.isfinite(x).all()
                    key = d['order_contrast'], d['sequence_axis'], d['degree']; assert key not in seen; seen.add(key)
                    db.execute('INSERT INTO designs VALUES(?,?,?,?,?,?,?,?,?)', (d['design_id'], c['cohort_id'], c['guide'], c['mask'],
                        d['order_contrast'], d['sequence_axis'], d['degree'], len(selected), d['disposition']))
                    design_ids.append(d['design_id']); design_counts[d['disposition']] += 1; design_n += 1
                    for outcome in OUTCOMES:
                        f = next(fits); y = values[outcome][selected]
                        assert f['design_id'] == d['design_id'] and f['cohort_id'] == c['cohort_id']
                        assert f['records'] == len(selected) and f['outcome'] == outcome and f['trees'] == plan['trees']
                        assert f['response_sha256'] == array_digest(y, '<f8')
                        assert f['fit_input_id'] == fit_id(d['design_id'], outcome, f['response_sha256'])
                        assert np.isfinite(y).all() and f['response_min'] == float(y.min()) and f['response_max'] == float(y.max())
                        status = d['disposition'] if d['disposition'] != 'full_rank_design' else (
                            'constant_response_requires_review' if np.ptp(y) == 0 else 'ready_for_working_covariance_fit')
                        assert f['disposition'] == status
                        db.execute('INSERT INTO fits VALUES(?,?,?,?)', (f['fit_input_id'], d['design_id'], outcome, status))
                        fit_ids.append(f['fit_input_id']); fit_counts[status] += 1; fit_n += 1
                assert seen == set(itertools.product(ORDERS, AXES, DEGREES))
                record = dict(cohort_id=c['cohort_id'], records=len(selected), guide=c['guide'], mask=c['mask'],
                    cohort_rows_sha256=c['case_rows_sha256'], ordered_case_ids_sha256=c['ordered_case_ids_sha256'],
                    control_record_sha256=entry['record_sha256'], control_array_sha256=entry['array_sha256'],
                    design_ids=design_ids, fit_input_ids=fit_ids, policy_basis_choices=choices,
                    prospective_numerical_audits=30 * 2 * len(plan['trees']) * 4,
                    raw_reml_basis_qualification_complete=False, scientific_eligibility=False)
                if reader: assert next(exported) == record
                else: out.write(json.dumps(record, sort_keys=True, allow_nan=False) + '\n')
                db.commit()
                if (i + 1) % 100 == 0: print('reconstructed_original_weighted_source_cohorts', i + 1, flush=True)
        finally:
            if out is not None: out.close()
        assert next(designs, None) is next(fits, None) is None
        if reader: assert next(exported, None) is None
        with gzip.open(source['root'] / 'model_settings.tsv.gz', 'rt') as f:
            settings = csv.DictReader(f, delimiter='\t'); assert settings.fieldnames == SETTING_FIELDS
            for s in settings:
                d = db.execute('SELECT cid,guide,mask,ord,axis,degree,n,status FROM designs WHERE did=?', (s['design_id'],)).fetchone()
                assert d is not None and d[:7] == (s['cohort_id'], s['guide'], s['mask'], s['order_contrast'], s['sequence_axis'], int(s['degree']), int(s['records']))
                fit = db.execute('SELECT did,outcome,status FROM fits WHERE fid=?', (s['fit_input_id'],)).fetchone()
                assert fit == (s['design_id'], s['outcome'], s['disposition'])
                assert int(s['nominal_tree_fits']) == len(plan['trees'])
                setting_counts[s['disposition']] += 1; settings_n += 1
        assert db.execute('SELECT COUNT(*) FROM designs').fetchone()[0] == design_n
        assert db.execute('SELECT COUNT(*) FROM fits').fetchone()[0] == fit_n
        db.close()
    expected = plan['expected']
    assert (design_n, fit_n, settings_n, rows_n) == (expected['designs'], expected['fit_inputs'], expected['settings'], expected['case_row_occurrences'])
    summary = dict(logical_cases=len(source['ids']), cohorts=len(source['cohorts']), case_row_occurrences=rows_n,
        designs=design_n, fit_inputs=fit_n, settings=settings_n, design_status_counts=dict(design_counts),
        fit_input_status_counts=dict(fit_counts), setting_status_counts=dict(setting_counts), policy_basis_counts=dict(basis_counts),
        expected_numerical_audit_rows=design_n * 2 * 5 * 4,
        expected_setting_audit_links=settings_n * 2 * 5 * 4, policies=POLICIES, trees=plan['trees'], loading_modes=MODES)
    assert summary['expected_numerical_audit_rows'] == expected['numerical_audit_rows']
    assert summary['expected_setting_audit_links'] == expected['setting_audit_links']
    if reader:
        assert all(producer[k] == v for k, v in summary.items()); bind(bindings, root / 'receipt.json')
    verify(bindings)
    artifacts = {n: sha(root / n) for n in ['stage_plan.json', 'cohort_source_census.jsonl.gz']}
    result = dict(status=READER if reader else PRODUCER, checked_utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=sha(path), source_contract=source['contract'], **summary, source_hashes=bindings,
        artifacts=artifacts, numerical_audits_computed=0, working_model_fits_computed=0,
        raw_reml_basis_qualification_complete=False, scientific_eligibility=False, scope=plan['scope'])
    if reader: result['producer_receipt_sha256'] = sha(root / 'receipt.json')
    write(root / name, result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes', 'artifacts', 'scope']}), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--reader', action='store_true'); a = p.parse_args(); run(a.plan, a.reader)
