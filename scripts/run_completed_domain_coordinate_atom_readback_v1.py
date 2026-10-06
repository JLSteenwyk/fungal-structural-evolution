#!/usr/bin/env python3
"""Run the unchanged full atom reader behind exact original native closure."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    parser.add_argument('--receipt', required=True, type=Path)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text())
    pins = dict(plan['pins'])
    verify(pins)
    closure_path = Path(plan['original_native_closure'])
    closure = json.loads(closure_path.read_text())
    assert closure['status'] == 'verified_original_domain_native_completion_api_terminal_unavailable'
    assert closure['original_native_exit_code'] == 0 and closure['original_native_exit_proof_available']
    assert closure['original_tool_session_id'] == 64749 and closure['original_tool_terminal_exit_code'] is None
    assert closure['original_tool_terminal_available'] is False
    assert closure['whole_wrapper_initial_and_terminal_payloads_matched']
    assert closure['manager_start_records'] == closure['manager_completion_records'] == 1
    assert closure['shards'] == plan['expected_shards'] == 977
    assert closure['counts']['models'] == plan['expected_models'] == 976357
    assert closure['counts']['intervals'] == plan['expected_intervals'] == 2454565
    assert closure['native_work_repeated'] is False
    transport = json.loads(Path(plan['original_native_closure_transport']).read_text())
    assert transport['status'] == 'verified_actual_original_review_wait_whole_wrapper_optional_cpu_accounting'
    assert transport['validation_sha256'] == sha(closure_path)
    assert transport['original_tool_session_id'] == 92757
    assert transport['original_tool_terminal_exit_code'] == 0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    assert transport['actual_manager_wait_terminal_verified']
    assert transport['optional_manager_cpu_record_fabricated'] is False
    for p, h in closure['source_hashes'].items():
        bind(pins, p, h)
    verify(pins)
    fixture = json.loads(Path(plan['reader_fixture']).read_text())
    assert fixture['status'] == 'passed_full_fixture_archive_atom_readback_and_mutation_rejection'
    assert fixture['script_sha256'] == sha('scripts/readback_domain_coordinate_archives.py')
    assert fixture['coordinate_mutation_rejected']
    # The qualified all-atom implementation, policy scope and rounding
    # tolerances are unchanged. Its full source/shard/CIF checks run once.
    subprocess.run([sys.executable, 'scripts/readback_domain_coordinate_archives.py',
                    '--plan', str(args.plan)], check=True)
    output = Path(plan['output'])
    raw_path = output / 'receipt.json'
    raw = json.loads(raw_path.read_text())
    assert raw['status'] == 'passed_all_exported_domain_atoms_and_full_disposition_scope'
    assert raw['plan_sha256'] == sha(args.plan)
    assert raw['shards'] == 977 and raw['counts'] == closure['counts']
    producer = json.loads(Path(plan['producer_plan']).read_text())
    assert raw['producer_receipt_sha256'] == sha(Path(producer['output']) / 'receipt.json')
    for p in [args.plan, raw_path, output / 'state.json', closure_path, Path(__file__),
              Path('scripts/readback_domain_coordinate_archives.py'), Path(plan['reader_fixture'])]:
        bind(pins, p)
    verify(pins)
    result = dict(status='complete_full_refreshed_domain_all_atom_readback_pending_original_reader_closure',
        checked_utc=datetime.now(timezone.utc).isoformat(), counts=raw['counts'], shards=raw['shards'],
        original_producer_tool_session_id=64749, original_producer_tool_terminal_exit_code=None,
        original_producer_native_exit_code=0, original_producer_api_terminal_unknown_preserved=True,
        original_producer_whole_native_execution_closed=True, full_manifest_job_scope_checked=True,
        every_exported_atom_checked=True, source_cif_bytes_and_shards_rehashed_by_reader=True,
        raw_receipt=str(raw_path), raw_receipt_sha256=sha(raw_path),
        coordinate_tolerance=0.000501, occupancy_confidence_tolerance=0.005001,
        confidence_summary_tolerance=1e-10, tolerances_unchanged=True,
        shared_cif_lexical_parser=True, rejection_causes_independently_adjudicated=False,
        source_hashes=pins, independent_atom_readback_complete=True,
        scientific_eligibility=False, gpu=False, new_predictions=0, new_cost_usd=0,
        scope='All976357models/2454565intervals/977shards checked by unchanged separate atom reader. '
              'Complete original producer native custody gates this run while its missing API terminal '
              'remains null. Every exported atom/sequence/identity/coordinate/occupancy/confidence and '
              'full manifest/job/archive scope replayed; recorded rejection reasons retained. '
              'Reader original API/native/source closure remains required. No PAE calibration, '
              'biological boundary/homology acceptance, extraction retry or GPU prediction.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
