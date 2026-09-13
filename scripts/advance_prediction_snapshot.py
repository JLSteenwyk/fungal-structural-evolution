#!/usr/bin/env python3
"""Audit and convert one full prediction queue after its pinned producer exits."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def identity(pid):
    try:
        raw = (Path('/proc') / str(pid) / 'stat').read_text()
    except FileNotFoundError:
        return None
    fields = raw[raw.rfind(')') + 2:].split()
    return None if fields[0] == 'Z' else fields[19]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    a = p.parse_args(); c = json.loads(a.config.read_text()); config_sha = sha(a.config)
    output = ROOT / c['control_output']; output.mkdir(parents=True, exist_ok=True)
    lock = (output / '.lock').open('w'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (output / 'launch.json').exists() or (output / 'receipt.json').exists():
        raise FileExistsError('Inspect previous execution before recovery')

    def verify():
        if sha(a.config) != config_sha:
            raise ValueError('Controller configuration changed')
        for name, expected in c['pinned_files'].items():
            if sha(ROOT / name) != expected:
                raise ValueError('Pinned dependency changed: ' + name)

    verify()
    print('Waiting for producer', c['producer_pid'], c['producer_start_ticks'], flush=True)
    while identity(c['producer_pid']) == c['producer_start_ticks']:
        time.sleep(20)
    verify()
    pred = ROOT / c['predictions']; chunk_path = pred / 'last_chunk.json'
    chunk = json.loads(chunk_path.read_text()); chunk_sha = sha(chunk_path)
    if (chunk.get('status') != 'production_chunk_finished'
            or chunk.get('config_sha256') != sha(pred / 'config.json')
            or chunk.get('interrupted') is not False or chunk.get('remaining_eligible') != 0
            or chunk.get('oom_deferred') != 0
            or chunk['cached_predictions'] + chunk['new_predictions'] != c['expected_predictions']):
        raise ValueError('Full prediction queue did not complete cleanly')
    if shutil.disk_usage(ROOT).free < c['minimum_free_disk_bytes']:
        raise ValueError('Insufficient output headroom')
    if any((ROOT / c[key]).exists() for key in ['audit_output', 'conversion_output']):
        raise FileExistsError('Use fresh immutable audit and conversion outputs')
    commands = [
        [sys.executable, 'scripts/audit_local_predictions.py', '--predictions', c['predictions'],
         '--inputs', c['inputs'], '--output', c['audit_output']],
        [sys.executable, 'scripts/convert_esmfold_snapshot.py', '--predictions', c['predictions'],
         '--audit', c['audit_output'], '--output', c['conversion_output']]]
    (output / 'launch.json').write_text(json.dumps({'commands': commands,
        'config_sha256': config_sha, 'prediction_chunk_sha256': chunk_sha,
        'resource_plan': c['resource_plan'], 'unix_time': time.time()}, indent=2) + '\n')
    env = os.environ.copy(); env['OPENBLAS_NUM_THREADS'] = '1'; env['OMP_NUM_THREADS'] = '1'
    for name, command in zip(['audit', 'conversion'], commands):
        verify()
        with (output / (name + '.log')).open('w') as log:
            subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    verify()
    if sha(chunk_path) != chunk_sha:
        raise ValueError('Prediction chunk changed during handoff')
    audit_path = ROOT / c['audit_output'] / 'receipt.json'
    conversion_path = ROOT / c['conversion_output'] / 'receipt.json'
    audit = json.loads(audit_path.read_text()); conversion = json.loads(conversion_path.read_text())
    if (audit['status'] != 'complete_artifact_readback' or audit['remaining_eligible'] != 0
            or audit['predictions'] != c['expected_predictions']
            or conversion['status'] != 'complete_audited_local_model_conversion'
            or conversion['models'] != c['expected_predictions']
            or conversion['source_audit_receipt_sha256'] != sha(audit_path)):
        raise ValueError('Audit/conversion universe mismatch')
    result = {'status': 'complete_prediction_snapshot_handoff', 'config_sha256': config_sha,
        'prediction_chunk_sha256': chunk_sha, 'audit_receipt_sha256': sha(audit_path),
        'conversion_receipt_sha256': sha(conversion_path), 'predictions': c['expected_predictions'],
        'interpretation': 'Full eligible queue audited and converted. Residue mapping, PAE binding, native features and evolutionary analyses remain separate.'}
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
