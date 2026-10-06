#!/usr/bin/env python3
"""Observe original independent coding/model reader identity and partial progress."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('launch', 'initial-tool', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    args = p.parse_args()
    assert not args.output.exists()
    launch, payload = [json.loads(q.read_text()) for q in (args.launch, args.initial_tool)]
    assert launch['original_tool_session_id'] == payload['original_tool_session_id'] == payload['initial']['session_id']
    assert launch['invocation_id'] in payload['initial']['output']
    assert launch['unit'] in payload['initial']['output']
    verify(launch['source_hashes'])
    assert sha(launch['plan']) == launch['plan_sha256']
    plan = json.loads(Path(launch['plan']).read_text())
    handle = observe(str(args.launch), {'cpu.max': '400000 100000',
                     'memory.max': str(64 * 2**30), 'memory.swap.max': '0'})
    execution = Path('metadata/full_coding_structure_source_coupling_readback_execution_20261005_v1')
    log = execution / 'stdout.log'
    completed = []
    if log.exists():
        # Only whole immutable-at-observation lines; a partial trailing line is not progress.
        blob = log.read_bytes()
        lines = blob[:blob.rfind(b'\n') + 1].decode().splitlines()
        for line in lines:
            m = re.fullmatch(r'independent_coding_structure_source_replay (\d+) /526', line)
            if m:
                completed.append(int(m.group(1)))
        assert completed == list(range(1, len(completed) + 1))
    pins = dict(launch['source_hashes'])
    for path in (args.launch, args.initial_tool, Path(__file__)):
        bind(pins, path)
    result = dict(status='verified_original_full_coding_structure_readback_runtime',
                  checked_utc=datetime.now(timezone.utc).isoformat(), original_handle=handle,
                  original_tool_session_id=payload['original_tool_session_id'], expected_taxa=526,
                  reported_completed_taxa=len(completed),
                  output_root_present=Path(plan['reader_output']).exists(), source_hashes=pins,
                  scientific_eligibility=False, inherited_translation_independently_recomputed=False,
                  biological_codon_eligibility=False, gpu=False, new_predictions=0,
                  scope='Exact original API/wrapper/native identity and cgroup, launch source bindings '
                        'and observed complete progress lines. Log counts do not prove native completion, '
                        'original API zero or full independent source/output closure. No retranslation, '
                        'gene/copy/taxonomy/codon/model or evolutionary admission, restart or mutation.')
    with args.output.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_hashes', 'original_handle')}, indent=2))


if __name__ == '__main__':
    main()
