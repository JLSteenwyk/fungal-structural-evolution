#!/usr/bin/env python3
"""Observe exact original assembly DNA acquisition without claiming full QC closure."""
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
    args = parser.parse_args(); assert not args.output.exists()
    launch, initial = [json.loads(p.read_text()) for p in (args.launch, args.initial_tool)]
    assert launch['original_tool_session_id'] == initial['original_tool_session_id'] == initial['initial']['session_id']
    assert launch['unit'] in initial['initial']['output'] and launch['invocation_id'] in initial['initial']['output']
    plan_path = Path(launch['plan']); assert sha(plan_path) == launch['plan_sha256']
    plan = json.loads(plan_path.read_text()); pins = dict(plan['pins']); verify(pins)
    resources_path = Path(launch['cmdline'][launch['cmdline'].index('--resources')+1])
    resources = json.loads(resources_path.read_text())
    limits = {'cpu.max': str(resources['cpus']*100000)+' 100000',
              'memory.max': str(resources['memory_gib']*2**30), 'memory.swap.max': '0'}
    assert launch['actual_cgroup_limits'] == limits
    handle = observe(str(args.launch), limits)
    root = Path(plan['output']); path = root/'state.json'; state = dict(present=path.exists())
    if state['present']:
        blob = path.read_bytes(); state.update(observed_bytes_sha256=hashlib.sha256(blob).hexdigest(), value=json.loads(blob))
    for path in (args.launch, args.initial_tool, plan_path, resources_path, Path(__file__)): bind(pins, path)
    result = dict(status='verified_original_full_assembly_dna_runtime', checked_utc=datetime.now(timezone.utc).isoformat(),
                  original_handle=handle, original_tool_session_id=launch['original_tool_session_id'],
                  producer_state=state, source_hashes=pins, scientific_eligibility=False, gpu=False, new_predictions=0,
                  scope='Exact original API initial identity, immutable plan/configuration sources and current '
                        'native process/cgroup; producer state is a partial observation. No original API terminal, '
                        'full genomic retrieval/statistic reconstruction, annotation or biology inferred. '
                        'This observation does not restart jobs or modify original publisher sources.')
    with args.output.open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes', 'original_handle')}, indent=2))


if __name__ == '__main__': main()
