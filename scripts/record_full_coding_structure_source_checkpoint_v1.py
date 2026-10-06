#!/usr/bin/env python3
"""Observe the original full coding/model source join without re-running it."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

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
    assert payload['original_tool_session_id'] == payload['initial']['session_id']
    assert launch['original_tool_session_id'] == payload['original_tool_session_id']
    assert launch['invocation_id'] in payload['initial']['output']
    assert launch['unit'] in payload['initial']['output']
    plan_path = Path(launch['plan'])
    assert sha(plan_path) == launch['plan_sha256']
    verify(launch['source_hashes'])
    plan = json.loads(plan_path.read_text())
    handle = observe(str(args.launch), {'cpu.max': '400000 100000',
                     'memory.max': str(64 * 2**30), 'memory.swap.max': '0'})
    root = Path(plan['output'])
    state_path = root / 'state.json'
    state = dict(present=state_path.exists())
    if state['present']:
        blob = state_path.read_bytes()
        state.update(observed_bytes_sha256=hashlib.sha256(blob).hexdigest(), value=json.loads(blob))
    pins = dict(launch['source_hashes'])
    for path in (args.launch, args.initial_tool, plan_path, Path(__file__)):
        bind(pins, path)
    result = dict(status='verified_original_full_coding_structure_source_coupling_runtime',
                  checked_utc=datetime.now(timezone.utc).isoformat(), original_handle=handle,
                  original_tool_session_id=payload['original_tool_session_id'],
                  expected_taxa=plan['expected_taxa'], source_products=plan['source_products'],
                  selected_representatives=plan['selected_representatives'],
                  target_records=plan['target_records'], producer_state=state,
                  output_root_present=root.exists(), source_hashes=pins,
                  scientific_eligibility=False, inherited_translation_independently_recomputed=False,
                  biological_codon_eligibility=False, gpu=False, new_predictions=0,
                  scope='Original API/wrapper/native identity and cgroup resources, immutable launch '
                        'and configuration bindings, partial source-join progress only. The full source '
                        'plan is bound here, not independently replayed by this observer. No original '
                        'API terminal, full source/output closure, independent every-record join replay, '
                        'translation, confidence, homology, taxonomy, copy or evolutionary admission. '
                        'No restart or mutation of original work.')
    with args.output.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_hashes', 'original_handle')}, indent=2))


if __name__ == '__main__':
    main()
