#!/usr/bin/env python3
"""Observe original full code-context producer or reader without restart."""
import argparse
from datetime import datetime, timezone
import hashlib
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
    assert launch['phase'] in ('producer', 'reader')
    assert launch['original_tool_session_id'] == payload['original_tool_session_id'] == payload['initial']['session_id']
    assert launch['invocation_id'] in payload['initial']['output'] and launch['unit'] in payload['initial']['output']
    verify(launch['source_hashes'])
    assert sha(launch['plan']) == launch['plan_sha256']
    plan = json.loads(Path(launch['plan']).read_text())
    handle = observe(str(args.launch), {'cpu.max': '400000 100000',
                     'memory.max': str(64 * 2**30), 'memory.swap.max': '0'})
    root = Path(plan['output' if launch['phase'] == 'producer' else 'reader_output'])
    state = dict(present=(root / 'state.json').exists())
    if state['present']:
        blob = (root / 'state.json').read_bytes()
        state.update(observed_bytes_sha256=hashlib.sha256(blob).hexdigest(), value=json.loads(blob))
    progress = []
    log = Path(launch['execution_root']) / 'stdout.log'
    if launch['phase'] == 'reader' and log.exists():
        blob = log.read_bytes()
        for line in blob[:blob.rfind(b'\n') + 1].decode().splitlines():
            m = re.fullmatch(r'independent_genetic_code_context_replay (\d+) /526', line)
            if m:
                progress.append(int(m.group(1)))
        assert progress == list(range(1, len(progress) + 1))
    pins = dict(launch['source_hashes'])
    for path in (args.launch, args.initial_tool, Path(__file__)):
        bind(pins, path)
    result = dict(status='verified_original_full_genetic_code_context_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), phase=launch['phase'], original_handle=handle,
        original_tool_session_id=payload['original_tool_session_id'], expected_taxa=526,
        expected_source_products=5927745, expected_target_records=5923039,
        output_root_present=root.exists(), producer_state=state, reported_reader_taxa=len(progress),
        final_receipt_present=Path(launch['receipt']).exists(), source_hashes=pins,
        scientific_eligibility=False, genetic_code_admission=False, biological_codon_eligibility=False,
        translations_recomputed=False, gpu=False, new_predictions=0,
        scope='Exact original API/wrapper/native identity and cgroup; immutable launch/configuration '
              'bindings and partial state/log observations. No full source replay by this observer, '
              'original API terminal, source/output closure, compartment inference, translation or '
              'genetic-code/selection/biological admission. No restart or mutation.')
    with args.output.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_hashes', 'original_handle')}, indent=2))


if __name__ == '__main__':
    main()
