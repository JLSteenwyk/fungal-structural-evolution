#!/usr/bin/env python3
"""Export all original-cohort controls or independently check every saved row."""
import argparse
import fcntl
import json
from pathlib import Path
import shutil

import numpy as np

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import rows, runtime_caps
from full_expanded_model_design_sources import array_digest
from full_inverse_reuse_weight_sources import load
from independent_inverse_reuse_weights import database, reconstruct, check_policy_records
from inverse_reuse_weight_controls import SCHEMA, POLICIES, ARRAYS, CLAIMS, calculate, policy_records
from reference_measurement_union_sources import bind, verify

PRODUCER = 'complete_full_original_cohort_inverse_reuse_controls_pending_readback_v1'
READER = 'passed_full_original_cohort_inverse_reuse_controls_sql_fraction_readback_v1'
SUMMARY = ['logical_cases', 'cohorts', 'case_row_occurrences', 'case_control_occurrences',
           'policies', 'policy_census', 'residual_diagonal_prepared',
           'raw_reml_basis_qualification_complete', 'nonuniform_weighting_accepted']
RECIPE = dict(weight='n/(represented_groups*cohort_local_reuse_count)',
              reciprocal_diagonal='(represented_groups*cohort_local_reuse_count)/n',
              uniform='weight=diagonal=1', counting_unit='one original logical case per cohort',
              group_scope='within original cohort after its existing quality/matching gates',
              interpretation='prespecified reuse sensitivity control; reciprocal diagonal is a working residual variance assumption',
              selection_record_multiplicity_used=False, confidence_used=False,
              physical_pair_scope='original versioned background_pair_key',
              family_scope='original connected family_component; not target family alone')


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.partial')
    assert not path.exists() and not temporary.exists()
    with temporary.open('x') as f:
        json.dump(value, f, sort_keys=True, allow_nan=False); f.write('\n')
    temporary.rename(path)


def record(cohort, contract, counts, weights, diagonals, groups, selected, npz_path):
    return dict(schema=SCHEMA, source_contract=contract, cohort_id=cohort['cohort_id'],
        guide=cohort['guide'], mask=cohort['mask'], records=len(selected),
        membership_occurrences=cohort['membership_occurrences'],
        cohort_rows_sha256=cohort['case_rows_sha256'],
        ordered_case_ids_sha256=cohort['ordered_case_ids_sha256'],
        original_cohort_npz_sha256=cohort['sha256'], array_file_sha256=sha(npz_path),
        array_sha256={name: array_digest(a, '<i8' if name in ARRAYS[:2] else '<f8')
            for name, a in zip(ARRAYS, [selected, counts, weights, diagonals])},
        policy_records=policy_records(len(selected), counts, weights, diagonals, groups),
        recipe=RECIPE, residual_diagonal_prepared=True, **CLAIMS)


def validate_saved(root, cohort, contract, selected, expected, entry=None):
    cid = cohort['cohort_id']; jp = root / 'cohorts' / (cid + '.json')
    ap = root / 'arrays' / (cid + '.npz'); saved = json.loads(jp.read_text())
    if entry is not None:
        assert entry == dict(cohort_id=cid, records=len(selected),
            record_path=str(jp.relative_to(root)), record_sha256=sha(jp),
            array_path=str(ap.relative_to(root)), array_sha256=sha(ap))
    with np.load(ap, allow_pickle=False) as z:
        assert z.files == ARRAYS
        arrays = [z[k] for k in ARRAYS]
    counts, weights, diagonals, groups = expected
    for name, a, b in zip(ARRAYS, arrays, [selected, counts, weights, diagonals]):
        assert a.dtype == np.dtype('int64' if name in ARRAYS[:2] else 'float64')
        assert a.shape == b.shape and np.array_equal(a, b), (cid, name)
    assert saved == record(cohort, contract, counts, weights, diagonals, groups, selected, ap)
    check_policy_records(saved['policy_records'], len(selected), counts, weights, diagonals, groups)
    return saved


def run(path, reader=False):
    runtime_caps()
    plan = json.loads(path.read_text()); source, contract, bindings = load(plan, path)
    root = Path(plan['output']); assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    assert not (root / ('readback.json' if reader else 'receipt.json')).exists()
    artifacts = {}; manifest = []; db = database(source) if reader else None
    census = {policy: dict(exact_uniform_one_cohorts=0, nonuniform_cohorts=0,
        represented_groups_min=None, represented_groups_max=0, maximum_reuse=0,
        weight_min=None, weight_max=0.) for policy in POLICIES}
    occurrences = 0
    if reader:
        rp = root / 'receipt.json'; producer = json.loads(rp.read_text())
        assert producer['status'] == PRODUCER and producer['plan_sha256'] == sha(path)
        assert producer['source_contract'] == contract and producer['source_hashes'] == bindings
        assert producer['scientific_eligibility'] is False and producer['recipe'] == RECIPE
        mp = root / 'control_manifest.json'; assert producer['artifacts'][mp.name] == sha(mp)
        manifest = json.loads(mp.read_text())
        assert [m['cohort_id'] for m in manifest] == [c['cohort_id'] for c in source['cohorts']]
    else:
        root.mkdir(exist_ok=True)
        lock = (root / 'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        stage = root / 'stage_plan.json'
        value = dict(plan_sha256=sha(path), source_contract=contract)
        if stage.exists(): assert json.loads(stage.read_text()) == value
        else: write_json(stage, value)
        (root / 'cohorts').mkdir(exist_ok=True); (root / 'arrays').mkdir(exist_ok=True)
    for i, cohort in enumerate(source['cohorts']):
        selected = rows(source, cohort)
        ap = root / 'arrays' / (cohort['cohort_id'] + '.npz')
        jp = root / 'cohorts' / (cohort['cohort_id'] + '.json')
        expected = (reconstruct(db, selected) if reader else
                    calculate([source['keys'][k][selected] for k in POLICIES[1:]]))
        if reader:
            saved = validate_saved(root, cohort, contract, selected, expected, manifest[i])
        else:
            if jp.exists(): saved = validate_saved(root, cohort, contract, selected, expected)
            else:
                assert not ap.exists()
                tmp = ap.with_suffix('.npz.partial'); assert not tmp.exists()
                with tmp.open('xb') as f:
                    np.savez_compressed(f, **dict(zip(ARRAYS, [selected, *expected[:3]])))
                tmp.rename(ap)
                saved = record(cohort, contract, *expected, selected, ap)
                write_json(jp, saved)
            manifest.append(dict(cohort_id=cohort['cohort_id'], records=len(selected),
                record_path=str(jp.relative_to(root)), record_sha256=sha(jp),
                array_path=str(ap.relative_to(root)), array_sha256=sha(ap)))
        for p in [ap, jp]: artifacts[str(p.relative_to(root))] = sha(p)
        for r in saved['policy_records']:
            c = census[r['policy']]
            c['exact_uniform_one_cohorts' if r['diagonal_is_exact_uniform_one'] else 'nonuniform_cohorts'] += 1
            g = r['represented_groups']
            c['represented_groups_min'] = g if c['represented_groups_min'] is None else min(c['represented_groups_min'], g)
            c['represented_groups_max'] = max(c['represented_groups_max'], g)
            c['maximum_reuse'] = max(c['maximum_reuse'], r['maximum_reuse'])
            c['weight_min'] = r['weight_min'] if c['weight_min'] is None else min(c['weight_min'], r['weight_min'])
            c['weight_max'] = max(c['weight_max'], r['weight_max'])
        occurrences += len(selected)
        if (i + 1) % 100 == 0: print('verified_controls' if reader else 'exported_controls', i + 1, flush=True)
    if db is not None: db.close()
    mp = root / 'control_manifest.json'
    if not reader: write_json(mp, manifest)
    artifacts[mp.name] = sha(mp)
    summary = dict(logical_cases=len(source['ids']), cohorts=len(source['cohorts']),
        case_row_occurrences=occurrences, case_control_occurrences=4 * occurrences,
        policies=POLICIES, policy_census=census, residual_diagonal_prepared=True,
        raw_reml_basis_qualification_complete=False, nonuniform_weighting_accepted=False)
    assert occurrences == plan['expected']['case_row_occurrences']
    if reader:
        assert producer['artifacts'] == artifacts
        assert all(producer[k] == v for k, v in summary.items())
        bind(bindings, rp)
        for name, h in artifacts.items(): bind(bindings, root / name, h)
    verify(bindings)
    result = dict(status=READER if reader else PRODUCER, plan_sha256=sha(path),
        source_contract=contract, **summary, recipe=RECIPE, source_hashes=bindings,
        scientific_eligibility=False, scope=plan['scope'])
    if reader: result['producer_receipt_sha256'] = sha(rp)
    else: result['artifacts'] = artifacts
    write_json(root / ('readback.json' if reader else 'receipt.json'), result)
    print(json.dumps(dict(status=result['status'], **summary)), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True); p.add_argument('--reader', action='store_true')
    a = p.parse_args(); run(a.plan, a.reader)
