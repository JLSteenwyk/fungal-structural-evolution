#!/usr/bin/env python3
"""Run a refreshed full longer-marker queue on an explicitly idle authorized GPU."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from Bio import SeqIO
from audit_busco_gene_copies import ROOT, sha
from assess_pae_sensitivity import checked_receipt
from advance_experimental_control_predictions import completed_chunk


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    a = p.parse_args(); c = json.loads(a.config.read_text()); config_sha = sha(a.config)
    control = ROOT / c['control_output']; control.mkdir(parents=True, exist_ok=True)
    lock = (control / '.lock').open('w'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (control / 'launch.json').exists() or (control / 'receipt.json').exists():
        raise FileExistsError('Review previous execution before recovery')
    for path, digest in c['pinned_files'].items():
        if sha(ROOT / path) != digest:
            raise ValueError('Changed launch dependency: ' + path)
    retired = json.loads((ROOT / c['superseded_controller_receipt']).read_text())
    if retired['status'] != 'superseded_before_prediction_launch':
        raise ValueError('Old queue must be retired before launch')
    inputs = ROOT / c['inputs']; original = ROOT / c['original_inputs']
    r, old = checked_receipt(inputs), checked_receipt(original)
    if r['source_receipts'] != old['source_receipts'] or r['selected_length_band_sequences'] != old['selected_length_band_sequences'] or r['artifacts']['all_marker_links.tsv'] != old['artifacts']['all_marker_links.tsv']:
        raise ValueError('Refreshed source universe differs')
    records = list(SeqIO.parse(inputs / 'candidates.faa', 'fasta'))
    if len(records) != r['prediction_candidates'] or not records or any(not 513 <= len(x.seq) <= 768 for x in records):
        raise ValueError('Invalid full longer-marker input grid')
    apps = subprocess.run(['nvidia-smi', '--query-compute-apps=gpu_uuid,pid', '--format=csv,noheader'], capture_output=True, text=True, check=True).stdout
    if c['gpu_uuid'] in {line.split(',')[0].strip() for line in apps.splitlines() if line.strip()}:
        raise RuntimeError('Selected GPU is occupied; no launch')
    if shutil.disk_usage(ROOT).free < c['minimum_free_disk_bytes']:
        raise RuntimeError('Insufficient output headroom')
    output = ROOT / c['prediction_output']
    if output.exists():
        raise FileExistsError('Use a fresh prediction output')
    command = [c['python'], str(ROOT / 'scripts/run_marker_predictions.py'), '--inputs', str(inputs),
               '--checkpoint', str(ROOT / c['checkpoint']), '--output', str(output),
               '--max-length', '768', '--limit', str(len(records))]
    launch = {'config_sha256': config_sha, 'command': command, 'prediction_candidates': len(records),
              'gpu_uuid': c['gpu_uuid'], 'gpu_process_observation': apps, 'unix_time': time.time(),
              'input_receipt_sha256': sha(inputs / 'receipt.json'), 'resource_plan_sha256': sha(inputs / 'resource_plan.json')}
    (control / 'launch.json').write_text(json.dumps(launch, indent=2) + '\n')
    print('Launching refreshed longer-marker queue', len(records), flush=True)
    env = os.environ.copy(); env['CUDA_VISIBLE_DEVICES'] = c['gpu_uuid']
    with (control / 'prediction.log').open('w') as log:
        run = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    chunk = output / 'last_chunk.json'
    success = run.returncode == 0 and chunk.exists() and completed_chunk(json.loads(chunk.read_text()), sha(output / 'config.json'))
    if sha(a.config) != config_sha or any(sha(ROOT / path) != digest for path, digest in c['pinned_files'].items()):
        raise ValueError('Pinned inputs changed during execution')
    result = {'status': 'complete_long_marker_prediction_chunk' if success else 'long_marker_prediction_requires_review',
              'exit_code': run.returncode, 'controller_config_sha256': config_sha,
              'launch_sha256': sha(control / 'launch.json'), 'prediction_chunk_sha256': sha(chunk) if chunk.exists() else None,
              'interpretation': 'Full queue execution accounting; model/PAE audits remain separate. OOM or interruption requires review.'}
    (control / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)
    if not success:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
