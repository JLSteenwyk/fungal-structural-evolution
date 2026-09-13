#!/usr/bin/env python3
"""Wait for the full resampling audit, then compute conditional joint tree-path summaries."""
import argparse
import fcntl
import json
import subprocess
import sys
import time
from pathlib import Path
from audit_busco_gene_copies import ROOT, sha


def process_identity(pid):
    try:
        value = Path(f'/proc/{pid}/stat').read_text()
    except FileNotFoundError:
        return None
    fields = value[value.rfind(')') + 2:].split()
    return {'start_ticks': fields[19], 'state': fields[0]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    a = p.parse_args(); c = json.loads(a.config.read_text()); config_hash = sha(a.config)
    out = ROOT / c['controller_output']; out.mkdir(parents=True, exist_ok=True)
    lock = (out / '.lock').open('w'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (out / 'receipt.json').exists():
        raise FileExistsError('Controller has already completed')
    def verify():
        if sha(a.config) != config_hash:
            raise ValueError('Controller configuration changed')
        for path, expected in c['pinned_files'].items():
            if sha(ROOT / path) != expected:
                raise ValueError('Pinned file changed: ' + path)
    verify()
    while True:
        current = process_identity(c['producer_pid'])
        if current is None or current['start_ticks'] != c['producer_start_ticks'] or current['state'] == 'Z':
            break
        print('waiting_for_producer', c['producer_pid'], flush=True)
        time.sleep(30)
    verify()
    receipt = ROOT / c['resampling_audit'] / 'receipt.json'
    if not receipt.exists():
        raise RuntimeError('Audit stopped without a complete receipt; path summaries not launched')
    producer = json.loads(receipt.read_text())
    if producer['status'] != 'complete_paired_resampling_audit':
        raise ValueError('Producer receipt does not establish completion')
    command = [sys.executable, str(ROOT / 'scripts/append_paired_path_uncertainty.py'),
               '--points', c['points'], '--inputs', c['inputs'], '--fits', c['fits'],
               '--resampling', c['resampling'], '--resampling-audit', c['resampling_audit'],
               '--output', c['path_output']]
    with (out / 'paths.log').open('w') as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    verify()
    audit = ROOT / c['path_output'] / 'receipt.json'
    result = json.loads(audit.read_text())
    if result['status'] != 'complete_joint_path_sampling_sensitivity':
        raise ValueError('Path summary did not establish completion')
    r = {'status': 'complete_joint_path_summary_handoff', 'config_sha256': config_hash,
         'script_sha256': sha(Path(__file__)), 'producer_receipt_sha256': sha(receipt),
         'path_receipt_sha256': sha(audit), 'command': command,
         'accepted_pairs': result['accepted_pairs'], 'excluded_pairs': result['excluded_pairs']}
    (out / 'receipt.json').write_text(json.dumps(r, indent=2) + '\n')
    print(json.dumps(r, indent=2), flush=True)


if __name__ == '__main__':
    main()
