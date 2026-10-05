#!/usr/bin/env python3
"""Run the unchanged full SQL-union readback after the retained-output review closes."""
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
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text())
    pins = dict(plan['pins'])
    verify(pins)
    producer_path = Path(plan['producer_receipt'])
    producer = json.loads(producer_path.read_text())
    transport_path = Path(plan['producer_transport'])
    transport = json.loads(transport_path.read_text())
    assert producer['status'] == 'retained_full_domain_manifest_after_wrapper_reporting_failure_pending_independent_readback'
    assert producer['original_tool_terminal_exit_code'] == 1
    assert producer['original_native_terminal_exit_code'] is None
    assert producer['original_wrapper_failure_preserved'] and not producer['native_work_repeated']
    assert transport['status'] == 'verified_original_software_wait_and_whole_wrapper_payloads'
    assert transport['validation_sha256'] == sha(producer_path)
    assert transport['original_tool_terminal_exit_code'] == 0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    assert transport['manager_start_records'] == transport['manager_completion_records'] == 1
    for mapping in (producer['source_hashes'], transport['source_hashes']):
        for path, digest in mapping.items():
            bind(pins, path, digest)
    verify(pins)
    raw_path = Path(plan['raw_readback'])
    assert not raw_path.exists()
    subprocess.run([sys.executable, 'scripts/readback_domain_extraction_manifest_v2.py',
                    '--manifest', plan['manifest'], '--database', plan['database'],
                    '--output', str(raw_path)], check=True)
    raw = json.loads(raw_path.read_text())
    assert raw['status'] == 'passed_full_domain_interval_manifest_sql_union_readback'
    assert raw['source_receipt_sha256'] == producer['raw_receipt_sha256']
    assert raw['intervals'] == producer['unique_intervals']
    assert raw['models'] == producer['source_models']
    assert raw['residues'] == producer['unique_interval_residues']
    for path in (args.plan, producer_path, transport_path, raw_path, Path(__file__)):
        bind(pins, path)
    verify(pins)
    result = dict(raw, checked_utc=datetime.now(timezone.utc).isoformat(),
                  source_hashes=pins, scientific_eligibility=False, gpu=False,
                  new_predictions=0, original_producer_tool_session_id=producer['original_tool_session_id'],
                  original_producer_tool_terminal_exit_code=1,
                  original_producer_native_terminal_exit_code=None,
                  retained_output_review_tool_session_id=transport['original_tool_session_id'],
                  raw_readback_sha256=sha(raw_path),
                  scope='Unchanged independent full SQL union against every interval, source model, '
                        'sequence, path and inclusive bound. Original producer wrapper reporting '
                        'failure retained; actual original retained-output review closed first. '
                        'No original native exit assumed. No coordinate, confidence, PAE, physical boundary or '
                        'evolutionary-event qualification.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
