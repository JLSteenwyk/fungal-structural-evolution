#!/usr/bin/env python3
"""Observe the original full domain atom reader without restarting or admitting it."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('launch', 'initial-tool', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    launch, payload = [json.loads(p.read_text()) for p in (args.launch, args.initial_tool)]
    assert payload['original_tool_session_id'] == payload['initial']['session_id']
    assert launch['invocation_id'] in payload['initial']['output']
    assert launch['unit'] in payload['initial']['output']
    plan_path = Path(launch['plan'])
    assert sha(plan_path) == launch['plan_sha256']
    plan = json.loads(plan_path.read_text())
    pins = dict(plan['pins'])
    verify(pins)
    handle = observe(str(args.launch), {'cpu.max': '800000 100000',
        'memory.max': str(64 * 2**30), 'memory.swap.max': '0'})
    state_path = Path(plan['output']) / 'state.json'
    state = dict(present=state_path.exists())
    if state['present']:
        blob = state_path.read_bytes()
        state['observed_bytes_sha256'] = hashlib.sha256(blob).hexdigest()
        try:
            state['value'] = json.loads(blob)
        except json.JSONDecodeError as error:
            state['partial_write_observation'] = str(error)
    for path in (args.launch, args.initial_tool, plan_path, Path(__file__)):
        bind(pins, path)
    result = dict(status='verified_original_full_domain_atom_readback_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), original_handle=handle,
        original_tool_session_id=payload['original_tool_session_id'],
        expected_models=plan['expected_models'], expected_intervals=plan['expected_intervals'],
        expected_shards=plan['expected_shards'], reader_state=state,
        reader_receipt_present=(Path(plan['output']) / 'receipt.json').exists(),
        source_hashes=pins, original_producer_api_terminal_unknown_preserved=True,
        independent_atom_readback_complete=False, scientific_eligibility=False,
        gpu=False, new_predictions=0,
        scope='Exact original wrapper/tool/plan identity and current resources or invocation-linked '
              'native terminal only. Mutable state is a partial observation, not terminal or full '
              'all-atom completion. Original producer API terminal remains unknown. No restart, '
              'prediction, boundary/homology acceptance or evolutionary event admission.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items()
                     if k not in ('source_hashes', 'original_handle')}, indent=2))


if __name__ == '__main__':
    main()
