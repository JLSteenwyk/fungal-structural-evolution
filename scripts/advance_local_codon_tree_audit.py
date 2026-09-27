#!/usr/bin/env python3
"""Wait for the identified local tree producer, then audit every planned fit."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)
    control = Path(plan['controller_output'])
    control.mkdir(parents=True, exist_ok=True)
    lock = (control / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (control / 'receipt.json').exists():
        raise FileExistsError('Audit controller already completed')

    def check_pins():
        if sha(args.plan) != plan_hash:
            raise ValueError('Controller plan changed')
        for path, expected in plan['pins'].items():
            if sha(path) != expected:
                raise ValueError('Pinned file changed: ' + path)

    check_pins()
    identity = plan['producer']
    while True:
        try:
            process = psutil.Process(identity['pid'])
            if process.create_time() != identity['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            if process.cmdline() != identity['cmdline']:
                raise ValueError('Producer command changed under the same identity')
        except psutil.NoSuchProcess:
            break
        print('waiting_for_identified_producer', identity['pid'], flush=True)
        time.sleep(30)
    check_pins()
    trees = Path(plan['trees'])
    producer_path = trees / 'receipt.json'
    result = json.loads(producer_path.read_text())
    if (result['status'] != 'complete_genus_nucleotide_tree_execution_pending_full_audit'
            or result['completed_cases'] != plan['expected_cases']
            or result['config_sha256'] != sha(trees / 'config.json')):
        raise ValueError('Producer lacks full expected completion')
    command = [sys.executable, 'scripts/audit_genus_codon_trees.py',
               '--trees', str(trees), '--inputs', plan['inputs'], '--output', plan['audit_output']]
    with (control / 'audit.log').open('x') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    check_pins()
    audit_root = Path(plan['audit_output'])
    audit = json.loads((audit_root / 'receipt.json').read_text())
    if (audit['status'] != 'passed_full_genus_tree_audit'
            or audit['audited_cases'] != plan['expected_cases']
            or audit['planned_cases'] != plan['expected_cases']
            or audit['pending_cases']
            or audit['bootstrap_trees_read'] != plan['expected_cases'] * 1000):
        raise ValueError('Incomplete independent tree audit')
    for path, expected in audit['artifacts'].items():
        if sha(audit_root / path) != expected:
            raise ValueError('Audit artifact changed: ' + path)
    receipt = dict(status='complete_full_local_codon_tree_audit_handoff',
                   plan_sha256=plan_hash, script_sha256=sha(__file__),
                   producer_receipt_sha256=sha(producer_path),
                   audit_receipt_sha256=sha(audit_root / 'receipt.json'),
                   audited_cases=audit['audited_cases'], bootstrap_trees_read=audit['bootstrap_trees_read'],
                   command=command,
                   scope='Full original audit applied to local alignments, without incomplete-case mode. SH-aLRT ranges checked, not recomputed. Topology sensitivity, model adequacy and selection eligibility remain separate.')
    (control / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
