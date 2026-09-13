#!/usr/bin/env python3
"""Project and normalize paired-site ASA after the pinned full audit finishes."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from advance_prediction_snapshot import ROOT, identity, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    config_sha = sha(args.config)
    control = ROOT / config['control_output']
    control.mkdir(parents=True, exist_ok=True)
    lock = (control / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if any((control / name).exists() for name in ['launch.json', 'receipt.json']):
        raise FileExistsError('Inspect previous execution before recovery')

    def verify():
        if sha(args.config) != config_sha:
            raise ValueError('Changed configuration')
        for name, expected in config['pinned_files'].items():
            if sha(ROOT / name) != expected:
                raise ValueError('Changed dependency: ' + name)

    verify()
    print('Waiting for full accessibility audit', config['predecessor_pid'], flush=True)
    while identity(config['predecessor_pid']) == config['predecessor_start_ticks']:
        time.sleep(20)
    verify()
    audit_path = ROOT / config['audit'] / 'receipt.json'
    audit = json.loads(audit_path.read_text())
    if (audit['status'] != 'passed_full_accessibility_snapshot'
            or audit['models_audited'] != config['models']
            or audit['residues_audited'] != config['residues']
            or audit['assessment_receipt_sha256'] != sha(ROOT / config['accessibility'] / 'receipt.json')
            or audit['snapshot_receipt_sha256'] != sha(ROOT / config['snapshot'] / 'receipt.json')):
        raise ValueError('Full matching audit required')
    for key in ['projection_output', 'normalization_output']:
        if (ROOT / config[key]).exists():
            raise FileExistsError('Use fresh immutable outputs')
    available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                         if line.startswith('MemAvailable:'))) * 1024
    if (available < config['resource_plan']['minimum_available_memory_bytes']
            or shutil.disk_usage(ROOT).free < config['resource_plan']['output_allowance_bytes']):
        raise RuntimeError('Insufficient resource headroom')
    commands = [
        [sys.executable, 'scripts/link_paired_sites_accessibility.py', '--inputs', config['inputs'],
         '--snapshot', config['snapshot'], '--accessibility', config['accessibility'],
         '--audit', config['audit'], '--output', config['projection_output']],
        [sys.executable, 'scripts/normalize_paired_accessibility.py',
         '--projection', config['projection_output'], '--snapshot', config['snapshot'],
         '--config', config['normalization_config'], '--output', config['normalization_output']]]
    audit_sha = sha(audit_path)
    launch = {'commands': commands, 'config_sha256': config_sha,
              'audit_receipt_sha256': audit_sha, 'resource_plan': config['resource_plan']}
    (control / 'launch.json').write_text(json.dumps(launch, indent=2) + '\n')
    env = os.environ.copy()
    env['OPENBLAS_NUM_THREADS'] = env['OMP_NUM_THREADS'] = '1'
    for name, command in zip(['projection', 'normalization'], commands):
        verify()
        if sha(audit_path) != audit_sha:
            raise ValueError('Audit changed')
        with (control / (name + '.log')).open('w') as log:
            subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    verify()
    projection_path = ROOT / config['projection_output'] / 'receipt.json'
    normalization_path = ROOT / config['normalization_output'] / 'receipt.json'
    projection = json.loads(projection_path.read_text())
    normalization = json.loads(normalization_path.read_text())
    if (sha(audit_path) != audit_sha
            or projection['status'] != 'complete_audited_accessibility_projection'
            or projection['markers'] != config['ready_markers']
            or projection['source_receipts']['audit'] != audit_sha
            or normalization['status'] != 'complete_reference_normalization'
            or normalization['rows'] != projection['observed_sites_linked']
            or normalization['source_receipt_sha256'] != sha(projection_path)):
        raise ValueError('Output receipts disagree')
    receipt = {'status': 'complete_accessibility_handoff_pending_independent_readback',
               'config_sha256': config_sha, 'audit_receipt_sha256': audit_sha,
               'projection_receipt_sha256': sha(projection_path),
               'normalization_receipt_sha256': sha(normalization_path),
               'observed_sites': normalization['rows'],
               'interpretation': 'Projection and normalization commands completed; independent row readback and controlled evolutionary analyses remain separate.'}
    (control / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
