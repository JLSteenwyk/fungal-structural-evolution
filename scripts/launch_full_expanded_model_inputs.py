#!/usr/bin/env python3
"""Launch full expanded model inputs, independent readback and provenance closure."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from full_expanded_model_input_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch, create
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args(); plan = json.loads(a.plan.read_text())
    for path, digest in plan['pins'].items(): assert sha(path) == digest, path
    expected_status = {
        'cases': 'complete_verified_full_matching_logical_case_index',
        'catalog': 'complete_verified_full_expanded_directed_measurement_catalog',
        'measurements': 'complete_verified_full_expanded_matched_case_measurements',
        'covariance': 'complete_verified_full_expanded_covariance'}
    for name, status in expected_status.items():
        c = json.loads(Path(plan[name + '_completion']).read_text())
        assert c['status'] == status and c['exact_process_journals_checked'] == 2
        assert sha(c['full_hash_archive']) == c['full_hash_archive_sha256']
    f = json.loads(Path(plan['fixture_validation']).read_text())
    assert f['status'] == 'passed_full_expanded_model_input_software_contracts'
    assert len(f['rejected_rehashed_exports']) == 18
    assert f['completed_restart_refused'] and f['full_interrupt_replay_passed']
    assert f['source_and_journal_fixtures_synthetic']
    for path, digest in f['source_hashes'].items(): assert sha(path) == digest, path
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    assert psutil.virtual_memory().available >= 32 * 2**30
    root = Path(plan['output']); assert root.parent.is_dir() and not root.exists()
    prefix = 'metadata/full_expanded_model_inputs'; name = 'full-expanded-model-inputs'
    producer = launch(name, prefix + '_wait_plan_20261002.json',
        [sys.executable, 'scripts/prepare_full_expanded_model_inputs.py', '--plan', str(a.plan)], [], str(a.plan))
    reader = launch(name + '-readback', prefix + '_readback_wait_plan_20261002.json',
        [sys.executable, 'scripts/readback_full_expanded_model_inputs.py', '--plan', str(a.plan),
         '--output', str(root / 'readback.json')], [producer], str(a.plan))
    cp = prefix + '_completion_plan_20261002.json'
    create(cp, dict(source_plan=str(a.plan), producer_receipt=str(root / 'receipt.json'),
        independent_readback=str(root / 'readback.json'),
        producer_status='complete_full_expanded_model_inputs_pending_independent_readback',
        reader_status='passed_full_expanded_model_input_decimal_and_sql_readback',
        completed_status='complete_verified_full_expanded_model_inputs', summary_fields=SUMMARY_FIELDS,
        launches=[producer, reader], pins={str(x): sha(x) for x in [a.plan, producer, reader,
            'scripts/close_full_triad_sequence_stage.py', 'scripts/record_completed_process_handoffs_v2.py']},
        output=prefix + '_completed_20261002.json', scope=plan['scope']))
    closure = launch(name + '-closure', prefix + '_closure_wait_plan_20261002.json',
        [sys.executable, 'scripts/close_full_triad_sequence_stage.py', '--plan', cp], [producer, reader], cp)
    create(prefix + '_launches_20261002.json', dict(status='launched_full_expanded_model_inputs',
        source_plan=str(a.plan), source_plan_sha256=sha(a.plan), launches=[producer, reader, closure], scope=plan['scope']))


if __name__ == '__main__': main()
