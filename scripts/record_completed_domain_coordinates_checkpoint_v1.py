#!/usr/bin/env python3
"""Observe the exact original refreshed extraction without restarting or qualifying it."""
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
    parser.add_argument('--launch', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    launch = json.loads(args.launch.read_text())
    plan_path = Path(launch['plan'])
    assert sha(plan_path) == launch['plan_sha256']
    plan = json.loads(plan_path.read_text())
    pins = dict(plan['pins'])
    verify(pins)
    handle = observe(str(args.launch), {'cpu.max': '400000 100000',
                     'memory.max': str(32*2**30), 'memory.swap.max': '0'})
    root = Path(plan['output'])
    state_path = root/'state.json'
    state = dict(present=state_path.exists())
    if state['present']:
        blob = state_path.read_bytes()
        state['observed_bytes_sha256'] = hashlib.sha256(blob).hexdigest()
        try:
            state['value'] = json.loads(blob)
        except json.JSONDecodeError as error:
            state['partial_write_observation'] = str(error)
    for path in (args.launch, plan_path, Path(__file__)):
        bind(pins, path)
    result = dict(status='verified_original_full_refreshed_domain_extraction_runtime',
                  checked_utc=datetime.now(timezone.utc).isoformat(), original_handle=handle,
                  expected_models=plan['expected_models'], expected_intervals=plan['expected_intervals'],
                  output_root_present=root.exists(), producer_state=state,
                  immutable_shard_receipt_files=len(list((root/'shards').glob('*.receipt.json'))),
                  raw_producer_receipt_present=(root/'receipt.json').exists(),
                  source_hashes=pins, scientific_eligibility=False, gpu=False, new_predictions=0,
                  scope='Exact original identity/cgroup or invocation-linked terminal evidence; '
                        'frozen plan sources and current mutable state/receipt-file census only. '
                        'No actual original tool terminal, full archive/atom readback, confidence/PAE '
                        'or biological acceptance is inferred. No restart, new prediction or mutation.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_hashes', 'original_handle')}, indent=2))


if __name__ == '__main__':
    main()
