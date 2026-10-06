#!/usr/bin/env python3
"""Review a completed retained atom readback with explicit omitted-zero counts."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def normalize_counts(observed, expected):
    if not isinstance(observed, dict) or not isinstance(expected, dict):
        raise ValueError('Counts must be dictionaries')
    if set(observed) - set(expected):
        raise ValueError('Unknown readback counter')
    if set(expected) - set(observed) not in (set(), {'rejected'}):
        raise ValueError('Missing required readback counter')
    for key, value in expected.items():
        if type(value) is not int or value < 0:
            raise ValueError('Invalid producer counter')
        if key not in observed and value != 0:
            raise ValueError('Missing nonzero rejection counter')
        actual = observed.get(key, 0)
        if type(actual) is not int or actual < 0 or actual != value:
            raise ValueError('Readback counter differs from producer')
    return {key: observed.get(key, 0) for key in expected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('reader-plan', 'reader-execution', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    plan = json.loads(args.reader_plan.read_text())
    pins = dict(plan['pins'])
    verify(pins)
    execution = json.loads(args.reader_execution.read_text())
    assert execution['status'] == 'failed_stage_retained'
    assert execution['exit_code'] == 1 and execution['timed_out'] is False
    assert execution['receipt_sha256'] is None
    assert not Path(execution['receipt']).exists()
    assert execution['invocation_id'] == '07c83fa9e4b54a50800c1aebc9a6c2d8'
    assert execution['wrapper']['pid'] == 749172
    assert execution['wrapper']['created'] == 1791259825.33
    assert 'scripts/run_completed_domain_coordinate_atom_readback_v1.py' in execution['command']
    assert str(args.reader_plan) in execution['command']
    stderr = args.reader_execution.with_suffix('') / 'stderr.log'
    error = stderr.read_text()
    assert "assert raw['shards'] == 977 and raw['counts'] == closure['counts']" in error
    assert error.rstrip().endswith('AssertionError')
    raw_path = Path(plan['output']) / 'receipt.json'
    raw = json.loads(raw_path.read_text())
    assert raw['status'] == 'passed_all_exported_domain_atoms_and_full_disposition_scope'
    assert raw['plan_sha256'] == sha(args.reader_plan)
    assert raw['shards'] == plan['expected_shards'] == 977
    producer_plan = json.loads(Path(plan['producer_plan']).read_text())
    producer_path = Path(producer_plan['output']) / 'receipt.json'
    producer = json.loads(producer_path.read_text())
    assert raw['producer_receipt_sha256'] == sha(producer_path)
    assert producer['counts']['models'] == plan['expected_models'] == 976357
    assert producer['counts']['intervals'] == plan['expected_intervals'] == 2454565
    counts = normalize_counts(raw['counts'], producer['counts'])
    assert 'rejected' not in raw['counts'] and producer['counts']['rejected'] == 0
    for mapping in (execution['source_hashes'], execution['artifacts']):
        for path, digest in mapping.items():
            bind(pins, path, digest)
    for path in (args.reader_plan, args.reader_execution, raw_path, producer_path, Path(__file__)):
        bind(pins, path)
    verify(pins)
    result = dict(status='verified_retained_full_atom_readback_zero_count_reporting_review_pending_execution_closure',
        checked_utc=datetime.now(timezone.utc).isoformat(), counts=counts,
        original_raw_counts=raw['counts'], omitted_zero_counter='rejected',
        shards=raw['shards'], original_reader_wrapper_exit_code=1,
        original_reader_api_terminal_exit_code=None, original_reader_api_terminal_not_reconstructed=True,
        native_atom_readback_repeated=False, original_failed_adapter_unchanged=True,
        original_source_and_numeric_tolerances_unchanged=True,
        full_atom_readback_receipt_verified=True, source_hashes=pins,
        scientific_eligibility=False, gpu=False, new_predictions=0,
        scope='Retained complete all-atom reader receipt reviewed after the original adapter fails '
              'on its omitted zero rejection key. Only the explicit producer-zero rejected counter '
              'can be filled; other missing, changed or foreign counters fail. No atom check is '
              'repeated. Actual original API/native/journal closure remains a separate requirement.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
