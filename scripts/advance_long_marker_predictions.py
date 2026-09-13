#!/usr/bin/env python3
"""Run longer-marker predictions after controls, with fresh reuse and resource checks."""
import argparse
import fcntl
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import ROOT, sha
from assess_pae_sensitivity import checked_receipt
from advance_experimental_control_predictions import process_identity, completed_chunk


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    args = p.parse_args(); config_hash = sha(args.config)
    c = json.loads(args.config.read_text())
    control = ROOT / c['control_output']; control.mkdir(parents=True, exist_ok=True)
    lock = (control / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (control / 'launch.json').exists() or (control / 'receipt.json').exists():
        raise FileExistsError('Inspect the previous execution before recovery')

    def verify():
        if sha(args.config) != config_hash:
            raise ValueError('Controller configuration changed')
        for name, digest in c['pinned_files'].items():
            if sha(ROOT / name) != digest:
                raise ValueError('Pinned dependency changed: ' + name)

    def wait_for_gpu():
        while True:
            apps = subprocess.run(['nvidia-smi', '--query-compute-apps=gpu_uuid,pid', '--format=csv,noheader'],
                                  check=True, text=True, capture_output=True).stdout
            occupied = {line.split(',')[0].strip() for line in apps.splitlines() if line.strip()}
            if c['gpu_uuid'] not in occupied:
                return apps
            time.sleep(20)

    verify()
    print('Waiting for pinned experimental-control controller', c['predecessor_pid'], flush=True)
    while process_identity(c['predecessor_pid']) == c['predecessor_start_ticks']:
        time.sleep(20)
    verify()
    previous_path = ROOT / c['predecessor_control'] / 'receipt.json'
    previous = json.loads(previous_path.read_text())
    if (previous['status'] != 'complete_experimental_control_prediction_chunk'
            or previous['controller_config_sha256'] != c['predecessor_config_sha256']):
        raise ValueError('Experimental controls did not finish successfully')
    previous_output = ROOT / c['predecessor_predictions']
    previous_chunk = previous_output / 'last_chunk.json'
    if (sha(previous_chunk) != previous['prediction_chunk_sha256']
            or not completed_chunk(json.loads(previous_chunk.read_text()), sha(previous_output / 'config.json'))):
        raise ValueError('Experimental-control prediction completion changed')
    wait_for_gpu()
    if shutil.disk_usage(ROOT).free < c['minimum_free_disk_bytes']:
        raise ValueError('Insufficient output headroom')
    inputs = ROOT / c['refreshed_inputs']
    refresh = [sys.executable, 'scripts/prepare_long_marker_predictions.py', '--inputs', *c['source_inputs'],
               '--predictions', *c['reuse_prediction_directories'], '--output', str(inputs)]
    with (control / 'input_refresh.log').open('w') as log:
        subprocess.run(refresh, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    receipt = checked_receipt(inputs)
    original = checked_receipt(ROOT / c['prepared_inputs'])
    if (receipt['status'] != 'complete_long_marker_input_preparation'
            or receipt['source_receipts'] != original['source_receipts']
            or receipt['selected_length_band_sequences'] != original['selected_length_band_sequences']
            or receipt['artifacts']['all_marker_links.tsv'] != original['artifacts']['all_marker_links.tsv']):
        raise ValueError('Refreshed input universe differs')
    records = list(SeqIO.parse(inputs / 'candidates.faa', 'fasta'))
    if len(records) != receipt['prediction_candidates'] or any(not 513 <= len(r.seq) <= 768 for r in records):
        raise ValueError('Refreshed prediction length/count mismatch')
    output = ROOT / c['prediction_output']
    if output.exists():
        raise FileExistsError('Use a fresh prediction output')
    if not records:
        result = {'status': 'complete_long_marker_queue_all_reusable',
                  'controller_config_sha256': config_hash, 'refreshed_input_receipt_sha256': sha(inputs / 'receipt.json')}
        (control / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
        return
    apps = wait_for_gpu(); verify()
    if shutil.disk_usage(ROOT).free < c['minimum_free_disk_bytes']:
        raise ValueError('Output headroom changed')
    command = [c['python'], str(ROOT / 'scripts/run_marker_predictions.py'), '--inputs', str(inputs),
               '--checkpoint', str(ROOT / c['checkpoint']), '--output', str(output),
               '--max-length', '768', '--limit', str(len(records))]
    launch = {'status': 'launching_long_marker_predictions', 'controller_config_sha256': config_hash,
              'predecessor_receipt_sha256': sha(previous_path), 'refreshed_input_receipt_sha256': sha(inputs / 'receipt.json'),
              'resource_plan_sha256': sha(inputs / 'resource_plan.json'), 'prediction_candidates': len(records),
              'command': command, 'gpu_uuid': c['gpu_uuid'], 'gpu_process_observation': apps, 'unix_time': time.time()}
    (control / 'launch.json').write_text(json.dumps(launch, indent=2) + '\n')
    print('Launching full longer-marker queue:', len(records), flush=True)
    env = os.environ.copy(); env['CUDA_VISIBLE_DEVICES'] = c['gpu_uuid']
    with (control / 'prediction.log').open('w') as log:
        run = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    chunk = output / 'last_chunk.json'
    success = run.returncode == 0 and chunk.exists() and completed_chunk(json.loads(chunk.read_text()), sha(output / 'config.json'))
    result = {'status': 'complete_long_marker_prediction_chunk' if success else 'long_marker_prediction_requires_review',
              'exit_code': run.returncode, 'controller_config_sha256': config_hash,
              'launch_sha256': sha(control / 'launch.json'), 'prediction_chunk_sha256': sha(chunk) if chunk.exists() else None,
              'interpretation': 'Clean execution accounting only; independent structure/PAE audit and biological analyses remain separate.'}
    (control / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)
    if not success:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
