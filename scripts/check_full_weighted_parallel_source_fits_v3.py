#!/usr/bin/env python3
"""Qualify complete fitting exports from parallel sources without real-data fits.

Complete export grids explicitly mock native fit outcomes. Separate retained
64-case actual-fit replay and guarded native producer/readback uses unchanged independent numerical arithmetic.
Private copied negatives never change qualified parent files.
"""
import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import shutil
from unittest.mock import patch

import numpy as np

from ancestral_chain_attempt import sha
from check_full_weighted_shared_entity_fits_v2 import actual_serialized, synthetic_candidate, synthetic_reader
SCHEMA = 'full-four-control-guarded-parallel-source-fit-v2'
from weighted_fit_numeric_memory_guard_v1 import ArrayGuard, guarded_call
from check_weighted_parallel_fitting_memory_v3 import memory_negatives, actual_guarded_native
from full_weighted_shared_entity_fit_sources_parallel_v1 import load, cohorts, cases
from prepare_full_weighted_parallel_source_fits_v2 import run
from readback_full_weighted_parallel_source_fits_v2 import run as readback
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify


def write(path, value):
    with Path(path).open('x') as handle:json.dump(value, handle, indent=2);handle.write('\n')


def overwrite(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n')


def rejected(action):
    try:action()
    except (AssertionError, ValueError, ArithmeticError, KeyError, FileExistsError, FileNotFoundError, StopIteration):return
    raise AssertionError('Malformed private fitting export or completed restart accepted')


def private_copy(parent_plan, destination, root):
    """Rebind headers/identities only in a copied serialized fixture."""
    plan = json.loads(parent_plan.read_text())
    parent = Path(plan['output'])
    shutil.copytree(parent, root)
    shutil.rmtree(root/'independent')
    (root/'readback.json').unlink()
    (root/'independent_readback_completed.json').unlink()
    plan.update(output=str(root), scope='Private copied fitting exports; no producer executed '
        'for this rebound plan. Mocked fit outcomes originate from the qualified parent.')
    write(destination, plan)
    source, bindings = load(plan, destination)
    stage = dict(schema=SCHEMA, plan_sha256=sha(destination), source_contract=source['fit_contract'])
    overwrite(root/'stage_plan.json', stage)
    manifest = json.loads((root/'cohort_manifest.json').read_text())
    for (cohort, rows, entries), part in zip(cohorts(source, plan), manifest):
        assert cohort['cohort_id'] == part['cohort_id']
        fp, lp, cp = [root/part[key] for key in ['candidates', 'links', 'receipt']]
        old = [json.loads(line) for line in gzip.decompress(fp.read_bytes()).decode().splitlines()]
        rebound = []
        mapping = {}
        for (expected, x, y, audit, route, diagonal), value in zip(cases(source, plan, cohort, entries), old):
            assert {k for k in expected if expected[k] != value[k]} <= {'candidate_id', 'source_contract'}
            mapping[value['candidate_id']] = expected['candidate_id']
            value.update(expected)
            rebound.append(value)
        assert len(rebound) == len(old) == 4800 and len(mapping) == 4800
        fp.write_bytes(gzip.compress((''.join(json.dumps(row)+'\n' for row in rebound)).encode(), mtime=0))
        lines = gzip.decompress(lp.read_bytes()).decode().splitlines()
        field = lines[0].split('\t').index('candidate_id')
        output = [lines[0]]
        for line in lines[1:]:
            cells = line.split('\t')
            cells[field] = mapping[cells[field]]
            output.append('\t'.join(cells))
        lp.write_bytes(gzip.compress(('\n'.join(output)+'\n').encode(), mtime=0))
        checkpoint = json.loads(cp.read_text())
        checkpoint.update(stage=stage, candidate_sha256=sha(fp), links_sha256=sha(lp))
        overwrite(cp, checkpoint)
        part.update(candidates_sha256=sha(fp), links_sha256=sha(lp), receipt_sha256=sha(cp))
    overwrite(root/'cohort_manifest.json', manifest)
    rp = root/'receipt.json'
    receipt = json.loads(rp.read_text())
    receipt.update(plan_sha256=sha(destination), source_contract=source['fit_contract'], source_hashes=bindings, scope=plan['scope'])
    receipt['artifacts'] = {name:sha(root/name) for name in receipt['artifacts']}
    overwrite(rp, receipt)
    return manifest


def corrupt(root, manifest, case):
    part = manifest[0]
    fp, lp, cp = [root/part[key] for key in ['candidates', 'links', 'receipt']]
    rp = root/'receipt.json'
    receipt = json.loads(rp.read_text())
    candidates = [json.loads(line) for line in gzip.decompress(fp.read_bytes()).decode().splitlines()]
    lines = gzip.decompress(lp.read_bytes()).decode().splitlines()
    checkpoint = json.loads(cp.read_text())
    if case == 'omit_candidate':candidates.pop()
    elif case == 'duplicate_candidate':candidates[-1] = deepcopy(candidates[0])
    elif case == 'foreign_candidate':candidates[0]['candidate_id'] = 'foreign'
    elif case == 'changed_diagonal':candidates[0]['diagonal_sha256'] = 'foreign'
    elif case == 'promote_science':candidates[0]['scientific_eligibility'] = True
    elif case == 'false_fit_failure':
        candidates[0].update(disposition='shared_entity_fit_error_requires_review', numerical_attempted=True,
            error_type='ArithmeticError', error_message='invented')
    elif case == 'missing_link':lines.pop()
    elif case == 'duplicate_link':lines[-1] = lines[1]
    elif case in ['foreign_link', 'changed_ordinal', 'changed_original_setting', 'changed_control']:
        field = {'foreign_link':'candidate_id', 'changed_ordinal':'source_setting_ordinal',
                 'changed_original_setting':'scenario_id', 'changed_control':'control_policy'}[case]
        fields, values = lines[0].split('\t'), lines[1].split('\t')
        values[fields.index(field)] = 'foreign'
        lines[1] = '\t'.join(values)
    elif case in ['cached_guard_false', 'cohort_guard_false', 'native_guard_false', 'cached_guard_count_zero', 'cohort_guard_count_zero']:
        field = {'cached_guard_false':'cached_numeric_inputs_preserved', 'cohort_guard_false':'cohort_numeric_inputs_preserved', 'native_guard_false':'native_call_numeric_inputs_checked', 'cached_guard_count_zero':'cached_input_array_bindings_checked', 'cohort_guard_count_zero':'cohort_input_array_bindings_checked'}[case]
        checkpoint[field] = 0 if case.endswith('zero') else False
    elif case in ['receipt_cached_guard_false', 'receipt_native_guard_false']:
        receipt['cached_numeric_inputs_preserved_for_all_cohorts' if case == 'receipt_cached_guard_false' else 'native_call_numeric_inputs_checked'] = False
    elif case == 'missing_artifact':receipt['artifacts'].pop(part['links'])
    elif case == 'changed_source_hash':receipt['source_hashes'][next(iter(receipt['source_hashes']))] = 'foreign'
    elif case == 'extra_artifact':receipt['artifacts']['invented'] = 'foreign'
    elif case == 'false_stage_contract':receipt['source_contract'] = 'foreign'
    else:raise AssertionError(case)
    fp.write_bytes(gzip.compress((''.join(json.dumps(row)+'\n' for row in candidates)).encode(), mtime=0))
    lp.write_bytes(gzip.compress(('\n'.join(lines)+'\n').encode(), mtime=0))
    checkpoint.update(candidate_sha256=sha(fp), links_sha256=sha(lp))
    overwrite(cp, checkpoint)
    part.update(candidates_sha256=sha(fp), links_sha256=sha(lp), receipt_sha256=sha(cp))
    overwrite(root/'cohort_manifest.json', manifest)
    for path in [fp, lp, cp, root/'cohort_manifest.json']:
        name = path.relative_to(root).as_posix()
        if name in receipt['artifacts']:receipt['artifacts'][name] = sha(path)
    overwrite(rp, receipt)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args()
    assert not a.receipt.exists()
    a.output.mkdir(exist_ok=False)
    gate = Path('metadata/weighted_parallel_fit_source_adapter_software_validation_20261004_v1.json')
    proven = json.loads(gate.read_text())
    verify(proven['source_hashes'])
    assert proven['serial_parallel_candidate_inputs_compared'] == 72000
    pins = {}
    project_sources(pins, [Path(__file__)])
    bind(pins, gate)
    fixture = Path('data/software_audits/full-weighted-parallel-fit-source-adapter-20261004-v1')
    completed = []
    restart_refusals = 0
    snapshots = {}
    with patch('prepare_full_weighted_parallel_source_fits_v2.candidate', side_effect=synthetic_candidate), \
         patch('readback_full_weighted_parallel_source_fits_v2.numeric', side_effect=synthetic_reader):
        for name in ['unqualified', 'qualified-common-pair', 'qualified-pair-exception']:
            parent = fixture/(name+'-parallel.fit.plan.json')
            fit = json.loads(parent.read_text())
            root = a.output/(name+'-exports')
            fit.update(output=str(root), resources=dict(minimum_free_disk_gib=128), pins={**fit['pins'], **pins},
                scope='Complete synthetic parallel-source fitting exports; source/journal closures '
                      'synthetic and native fit producer/reader branches explicitly mocked. No biological fit.')
            pp = a.output/(name+'.fit.plan.json')
            write(pp, fit)
            if name == 'unqualified':
                try:run(pp, stop_after_cohorts=2)
                except InterruptedError:pass
                else:raise AssertionError('Controlled producer checkpoint interruption absent')
                prior = {str(path):sha(path) for path in (root/'cohorts').glob('*') if path.is_file()}
            producer = run(pp)
            if name == 'unqualified':assert all(sha(path) == digest for path, digest in prior.items())
            try:readback(pp, root/'readback.json', stop_after_cohorts=2)
            except InterruptedError:pass
            else:raise AssertionError('Controlled reader checkpoint interruption absent')
            previous = {str(path):sha(path) for path in (root/'independent').glob('*') if path.is_file()}
            reader = readback(pp, root/'readback.json')
            assert all(sha(path) == digest for path, digest in previous.items())
            assert producer['candidate_rows'] == reader['candidate_rows'] == 24000
            assert producer['setting_fit_links'] == reader['setting_fit_links'] == 48000
            rejected(lambda:run(pp))
            rejected(lambda:readback(pp, root/'other-reader.json'))
            restart_refusals += 2
            completed.append(dict(name=name, candidate_rows=24000, setting_fit_links=48000,
                producer_status_counts=producer['candidate_status_counts'], independent_status_counts=reader['independent_candidate_status_counts']))
            for path in root.rglob('*'):
                if path.is_file():snapshots[str(path)] = sha(path)
            verify(reader['source_hashes'])
            for path, digest in reader['source_hashes'].items():bind(pins, path, digest)
        # Positive rebound copy demonstrates that plan/source identities have
        # been fully rebound before independent private negatives are meaningful.
        original = a.output/'unqualified.fit.plan.json'
        pp, root = a.output/'private-control.plan.json', a.output/'private-control'
        private_copy(original, pp, root)
        readback(pp, root/'readback.json')
        rejections = []
        captures = 0
        for case in ['omit_candidate', 'duplicate_candidate', 'foreign_candidate', 'changed_diagonal',
            'promote_science', 'false_fit_failure', 'missing_link', 'duplicate_link', 'foreign_link',
            'changed_ordinal', 'changed_original_setting', 'changed_control', 'missing_artifact',
            'changed_source_hash', 'extra_artifact', 'false_stage_contract', 'cached_guard_false', 'cohort_guard_false',
            'native_guard_false', 'cached_guard_count_zero', 'cohort_guard_count_zero',
            'receipt_cached_guard_false', 'receipt_native_guard_false']:
            pp, root = a.output/('private-'+case+'.plan.json'), a.output/('private-'+case)
            manifest = private_copy(original, pp, root)
            corrupt(root, manifest, case)
            rejected(lambda:readback(pp, root/'rejected-reader.json'))
            for failure in (root/'independent'/'failures').glob('*/failure.json'):
                assert json.loads(failure.read_text())['scientific_eligibility'] is False
                with np.load(failure.parent/'original_numeric_inputs.npz', allow_pickle=False) as values:
                    assert set(values.files) == {'case_rows','diagonal','design','response','species_factor'}
                captures += 1
            rejections.append(case)
    memory_rejections = memory_negatives(a.output, original, private_copy)
    actual_native_checks = actual_guarded_native(a.output)
    actual_pins, replays = actual_serialized(a.output)
    for path, digest in actual_pins.items():bind(pins, path, digest)
    verify(snapshots)
    # Preserve every private outcome and file inventory, including rejected
    # partial arrays and JSON; no negative or parent namespace is cleaned up.
    for path in a.output.rglob('*'):
        if path.is_file():bind(pins, path)
    verify(pins)
    result = dict(status='passed_complete_guarded_parallel_source_fitting_exports_and_numeric_reader_v3',
        checked_utc=datetime.now(timezone.utc).isoformat(), complete_synthetic_grids=completed,
        total_candidate_rows=72000, total_setting_fit_links=144000, actual_serialized_numeric_replays=len(replays),
        private_rehashed_corruptions_rejected=rejections, positive_rebound_private_control_passed=True,
        in_memory_mutations_rejected=memory_rejections, actual_guarded_native_candidate_and_reader_checks=actual_native_checks,
        exact_numeric_failure_captures=captures, qualified_parent_artifacts_preserved=True,
        producer_closed_chunks_reused_without_rewrite=True, reader_closed_chunks_replayed_without_rewrite=True,
        completed_role_restarts_refused=restart_refusals,
        full_export_native_fit_branches_explicitly_mocked=True, source_and_journal_fixtures_synthetic=True,
        actual_numeric_replays_use_unchanged_frozen_independent_reader=True, source_hashes=pins,
        production_fitting_launched=False, scientific_eligibility=False, gpu=False,
        scope='Three complete five-cohort parallel-source fitting export grids:72000candidates/144000links, '
              'explicitly mocked native fit branches, independent source/checkpoint inventory. '
              'Separate64retained actual native fits replay unchanged arithmetic. Private copied '
              'negative headers/identities rebound without claiming new producer execution; '
              'positive copy precedes twenty-three rejected corruptions; separate in-memory mutations are refused. Original qualified parents '
              'are byte preserved. No actual full biological fit or uncertainty calibration.')
    write(a.receipt, result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','scope']}, indent=2))


if __name__ == '__main__':main()
