#!/usr/bin/env python3
"""Run the full alignment/tree comparison only after the independent audit passes."""
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
    digest = sha(args.plan)
    control = Path(plan['controller_output'])
    control.mkdir(parents=True, exist_ok=True)
    lock = (control / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def verify():
        if sha(args.plan) != digest:
            raise ValueError('Plan changed')
        for path, expected in plan['pins'].items():
            if sha(path) != expected:
                raise ValueError('Pinned file changed: ' + path)

    verify()
    identity = plan['audit_controller']
    while True:
        try:
            process = psutil.Process(identity['pid'])
            if process.create_time() != identity['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            if process.cmdline() != identity['cmdline']:
                raise ValueError('Audit controller identity changed')
        except psutil.NoSuchProcess:
            break
        print('waiting_for_identified_audit_controller', identity['pid'], flush=True)
        time.sleep(30)
    verify()
    proof_path = Path(plan['audit_handoff_receipt'])
    proof = json.loads(proof_path.read_text())
    if (proof['status'] != 'complete_full_local_codon_tree_audit_handoff'
            or proof['audited_cases'] != 1632
            or proof['audit_receipt_sha256'] != sha(Path(plan['arguments']['local-audit']) / 'receipt.json')):
        raise ValueError('Full local audit has not passed')
    command = [sys.executable, 'scripts/compare_codon_alignment_trees.py']
    for key, value in plan['arguments'].items():
        command.extend(['--' + key, value])
    with (control / 'comparison.log').open('x') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    verify()
    output = Path(plan['arguments']['output'])
    result = json.loads((output / 'receipt.json').read_text())
    if result['cases'] != 1712 or result['dispositions'] != plan['expected_dispositions']:
        raise ValueError('Full comparison grid differs')
    receipt = dict(status='complete_codon_tree_comparison_execution_pending_readback',
                   plan_sha256=digest, audit_handoff_receipt_sha256=sha(proof_path),
                   comparison_receipt_sha256=sha(output / 'receipt.json'), command=command,
                   scope='Full comparison produced only after both source tree audits passed. Independent comparison-output readback remains required.')
    (control / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
