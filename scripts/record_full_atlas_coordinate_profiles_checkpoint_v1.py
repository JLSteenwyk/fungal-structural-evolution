#!/usr/bin/env python3
"""Observe original full-source profile audit identities, resources and progress only."""
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
        parser.add_argument('--'+name, type=Path, required=True)
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
                     'memory.max': str(64*2**30), 'memory.swap.max': '0'})
    root = Path(plan['output'])
    state_path = root/'state.json'
    state = dict(present=state_path.exists())
    if state['present']:
        blob = state_path.read_bytes()
        state.update(observed_bytes_sha256=hashlib.sha256(blob).hexdigest(), value=json.loads(blob))
    for path in (args.launch, args.initial_tool, plan_path, Path(__file__)):
        bind(pins, path)
    result = dict(status='verified_original_full_atlas_coordinate_profile_runtime',
                  checked_utc=datetime.now(timezone.utc).isoformat(), original_handle=handle,
                  original_tool_session_id=payload['original_tool_session_id'],
                  expected_models_by_source=plan['expected_models_by_source'],
                  expected_residues_by_source=plan['expected_residues_by_source'],
                  output_root_present=root.exists(), producer_state=state,
                  immutable_shard_receipt_files=len(list((root/'shards').glob('*.receipt.json'))),
                  raw_producer_receipt_present=(root/'receipt.json').exists(),
                  source_hashes=pins, scientific_eligibility=False, gpu=False, new_predictions=0,
                  scope='Exact original wrapper/tool identity and current cgroup resources, immutable '
                        'plan sources and partial producer state/file census. No original tool terminal, '
                        'independent full coordinate/array/archive/rejection readback, confidence/PAE '
                        'calibration, homology or biological inference is assumed. No restart or mutation.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_hashes', 'original_handle')}, indent=2))


if __name__ == '__main__':
    main()
