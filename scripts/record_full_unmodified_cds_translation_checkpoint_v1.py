#!/usr/bin/env python3
"""Observe the exact original full fixed-code translation job without restarting."""
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
    for name in ['launch', 'initial-tool', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    launch, payload = [json.loads(p.read_text()) for p in [args.launch, args.initial_tool]]
    assert launch['phase'] in ['producer', 'reader']
    assert payload['original_tool_session_id'] == payload['initial']['session_id'] == launch['original_tool_session_id']
    assert launch['invocation_id'] in payload['initial']['output'] and launch['unit'] in payload['initial']['output']
    verify(launch['source_hashes'])
    assert sha(launch['plan']) == launch['plan_sha256']
    plan = json.loads(Path(launch['plan']).read_text())
    assert plan['expected_taxa'] == 526 and plan['target_records'] == 5923039
    handle = observe(str(args.launch), {'cpu.max': '400000 100000',
                     'memory.max': str(64 * 2**30), 'memory.swap.max': '0'})
    root = Path(plan['output' if launch['phase'] == 'producer' else 'reader_output'])
    state = dict(present=(root / 'state.json').exists())
    if state['present']:
        blob = (root / 'state.json').read_bytes()
        state.update(observed_bytes_sha256=hashlib.sha256(blob).hexdigest(), value=json.loads(blob))
    pins = dict(launch['source_hashes'])
    for path in [args.launch, args.initial_tool, Path(__file__)]:
        bind(pins, path)
    result = dict(status='verified_original_full_unmodified_cds_translation_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), phase=launch['phase'], original_handle=handle,
        original_tool_session_id=launch['original_tool_session_id'], expected_taxa=526,
        expected_source_products=5927745, expected_target_records=5923039, output_root_present=root.exists(),
        observed_state=state, final_receipt_present=Path(launch['receipt']).exists(), source_hashes=pins,
        scientific_eligibility=False, genetic_code_admission=False, biological_codon_eligibility=False,
        gpu=False, new_predictions=0,
        scope='Exact original API/wrapper/native identity and cgroup; immutable launch/configuration '
              'bindings and partial state observations. No actual original API terminal, full independent '
              'source/output replay, genetic code/compartment adoption, codon alignment or biological result. '
              'Observation never restarts a run or changes its inputs.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'original_handle']}, indent=2))


if __name__ == '__main__':
    main()
