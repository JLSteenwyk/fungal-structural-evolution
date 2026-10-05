#!/usr/bin/env python3
"""Observe the exact original full genomic-CDS stage without repeating its source scan."""
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
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    launch, initial = [json.loads(p.read_text()) for p in (args.launch, args.initial_tool)]
    assert launch['original_tool_session_id'] == initial['original_tool_session_id'] == initial['initial']['session_id']
    assert launch['unit'] in initial['initial']['output'] and launch['invocation_id'] in initial['initial']['output']
    assert sha(launch['plan']) == launch['plan_sha256']
    verify(launch['source_hashes'])
    resources_path = Path(launch['cmdline'][launch['cmdline'].index('--resources') + 1])
    resources = json.loads(resources_path.read_text())
    expected = {'cpu.max': str(resources['cpus'] * 100000) + ' 100000',
                'memory.max': str(resources['memory_gib'] * 2**30), 'memory.swap.max': '0'}
    assert launch['actual_cgroup_limits'] == expected
    handle = observe(str(args.launch), expected)
    plan = json.loads(Path(launch['plan']).read_text())
    assert plan['expected_taxa'] == 526
    root = Path(plan['output'])
    state_path = root / 'state.json'
    state = dict(present=state_path.exists())
    if state['present']:
        blob = state_path.read_bytes()
        state.update(observed_bytes_sha256=hashlib.sha256(blob).hexdigest(), value=json.loads(blob))
    pins = dict(launch['source_hashes'])
    for path in [args.launch, args.initial_tool, resources_path, Path(__file__)]:
        bind(pins, path)
    result = dict(status='verified_original_full_genomic_cds_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), original_handle=handle,
        original_tool_session_id=launch['original_tool_session_id'], producer_state=state,
        output_root_present=root.exists(), immutable_taxon_receipts=len(list(root.glob('*/receipt.json'))),
        source_hashes=pins, scientific_eligibility=False, gpu=False, new_predictions=0,
        scope='Exact original launch/API identity and current native/cgroup state; immutable plan and '
              'execution-configuration bindings checked. Mutable counters/receipt census are partial '
              'observations, not full original API completion, independent genomic-CDS replay or biological '
              'acceptance. Original source hashes are checked by the producer; this observer does not '
              'repeat their full scan, restart jobs or alter original data.')
    with args.output.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes', 'original_handle')}, indent=2))


if __name__ == '__main__':
    main()
