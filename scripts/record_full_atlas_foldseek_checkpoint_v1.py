#!/usr/bin/env python3
"""Observe the original complete-atlas native conversion and coordinate readback."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import psutil

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
        'memory.max': str(96 * 2**30), 'memory.swap.max': '0'})
    root = Path(plan['output'])
    state_path = root / 'state.json'
    state = dict(present=state_path.exists())
    if state['present']:
        blob = state_path.read_bytes()
        state.update(observed_bytes_sha256=hashlib.sha256(blob).hexdigest(), value=json.loads(blob))
    native_path = root / 'createdb-original-process.json'
    native = dict(original_identity_present=native_path.exists())
    if native['original_identity_present']:
        identity = json.loads(native_path.read_text())
        native['original_identity'] = identity
        try:
            process = psutil.Process(identity['pid'])
            live = process.create_time() == identity['created'] and process.status() != psutil.STATUS_ZOMBIE
            if live:
                assert process.cmdline() == identity['command']
                assert process.environ()['CUDA_VISIBLE_DEVICES'] == ''
                assert identity['command'][identity['command'].index('--gpu')+1] == '0'
                native.update(original_process_live=True, pid=process.pid, created=process.create_time(),
                              status=process.status(), cpu_seconds=sum(process.cpu_times()[:2]))
            else:
                native['original_process_live'] = False
        except psutil.NoSuchProcess:
            native['original_process_live'] = False
        execution_path = root / 'createdb-execution.json'
        if execution_path.exists():
            execution = json.loads(execution_path.read_text())
            assert all(execution[key] == identity[key] for key in ('pid', 'created', 'command'))
            native['saved_original_native_execution'] = execution
            bind(pins, execution_path)
        bind(pins, native_path)
    for path in (args.launch, args.initial_tool, plan_path, Path(__file__)):
        bind(pins, path)
    result = dict(status='verified_original_full_atlas_foldseek_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), original_handle=handle,
        original_tool_session_id=payload['original_tool_session_id'], native_createdb=native,
        expected_models_by_source=plan['expected_models_by_source'],
        expected_residues_by_source=plan['expected_residues_by_source'], producer_state=state,
        geometry_proof_files=len(list((root / 'geometry-proofs').glob('*.json'))),
        full_producer_receipt_present=Path('metadata/full_atlas_foldseek_database_20261006_v2.json').exists(),
        source_hashes=pins, independent_whole_database_reader_complete=False,
        scientific_eligibility=False, gpu=False, new_predictions=0,
        scope='Exact original wrapper/tool/plan/current resources and native createdb identity '
              'or saved native execution observed. Mutable stage/counts and proof-file census '
              'are partial evidence. No full source-coordinate or native/database closure, '
              '3Di accuracy, confidence/PAE calibration, clustering, homology or evolutionary '
              'admission inferred. No restart or mutation.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_hashes', 'original_handle')}, indent=2))


if __name__ == '__main__':
    main()
