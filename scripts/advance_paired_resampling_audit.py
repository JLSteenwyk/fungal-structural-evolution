#!/usr/bin/env python3
"""Wait for a specific producer to finish, then audit its complete paired branch resampling output."""
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
    receipt = ROOT / c['resampling'] / 'receipt.json'
    if not receipt.exists():
        raise RuntimeError('Producer stopped without a complete receipt; audit not launched')
    producer = json.loads(receipt.read_text())
    if producer['status'] != 'complete_paired_resampling_execution':
        raise ValueError('Producer receipt does not establish completion')
    command = [sys.executable, str(ROOT / 'scripts/summarize_paired_resampling.py'),
               '--inputs', c['inputs'], '--fits', c['fits'],
               '--resampling', c['resampling'], '--output', c['audit_output']]
    with (out / 'audit.log').open('w') as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    verify()
    audit = ROOT / c['audit_output'] / 'receipt.json'
    result = json.loads(audit.read_text())
    if result['status'] != 'complete_paired_resampling_audit':
        raise ValueError('Audit did not establish complete resampling coverage')
    r = {'status': 'complete_full_paired_resampling_audit', 'config_sha256': config_hash,
         'script_sha256': sha(Path(__file__)), 'producer_receipt_sha256': sha(receipt),
         'audit_receipt_sha256': sha(audit), 'command': command,
         'attempted_paired_draws': result['attempted_paired_draws'],
         'validated_fits': result['validated_fits'], 'unestimable_draws': result['unestimable_draws']}
    (out / 'receipt.json').write_text(json.dumps(r, indent=2) + '\n')
    print(json.dumps(r, indent=2), flush=True)


if __name__ == '__main__':
    main()
