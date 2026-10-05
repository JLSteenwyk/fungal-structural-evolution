#!/usr/bin/env python3
"""Observe the exact original V10 grid without repeating inputs or restarting roles."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['launch', 'initial-tool', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    launch, payload = [json.loads(p.read_text()) for p in [args.launch, args.initial_tool]]
    assert launch['original_tool_session_id'] == payload['original_tool_session_id'] == payload['initial']['session_id']
    assert launch['unit'] in payload['initial']['output'] and launch['invocation_id'] in payload['initial']['output']
    verify(launch['source_hashes'])
    plan = json.loads(Path(launch['plan']).read_text())
    assert sha(launch['plan']) == launch['plan_sha256'] and plan['job_scope']['full_V10_matrix_roles'] == 24
    resources = plan['resources']
    expected = {'cpu.max': str(resources['cpus'] * 100000) + ' 100000',
                'memory.max': str(resources['memory_gib'] * 2**30), 'memory.swap.max': '0'}
    handle = observe(str(args.launch), expected)
    root = Path(plan['output'])
    chains = list((root / 'chains').glob('*.json'))
    outcomes = [json.loads(p.read_text()) for p in chains]
    assert len({r['chain_id'] for r in outcomes}) == len(outcomes) <= 24
    pins = dict(launch['source_hashes'])
    for path in [args.launch, args.initial_tool, Path(__file__), *chains]:
        bind(pins, path)
    result = dict(status='verified_original_V10_full_input_comparison_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), original_handle=handle,
        original_tool_session_id=launch['original_tool_session_id'], output_root_present=root.exists(),
        recorded_roles=len(outcomes), expected_roles=24,
        recorded_native_zero_roles=sum(r['exit_code'] == 0 for r in outcomes),
        recorded_native_failure_roles=sum(r['exit_code'] != 0 for r in outcomes),
        producer_receipt_present=(root / 'receipt.json').exists(), reader_receipt_present=(root / 'readback.json').exists(),
        source_hashes=pins, scientific_eligibility=False, posterior_qualified=False, gpu=False,
        scope='Exact original API identity, wrapper PID/create/command/invocation and current native/cgroup observations. '
              'Recorded chain census is partial producer state, not final paired-output or latent-value qualification. '
              'Does not repeat original source scan, restart attempts, use GPUs or admit posterior arrays.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'original_handle']}, indent=2))


if __name__ == '__main__':
    main()
