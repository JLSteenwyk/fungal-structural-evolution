#!/usr/bin/env python3
"""Wait for the verified full local-tree audit, then execute all planned MG94 fits."""
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
            raise ValueError('Controller plan changed')
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
                raise ValueError('Audit controller command changed')
        except psutil.NoSuchProcess:
            break
        print('waiting_for_full_tree_audit', identity['pid'], flush=True)
        time.sleep(30)
    verify()
    handoff_path = Path(plan['audit_handoff_receipt'])
    handoff = json.loads(handoff_path.read_text())
    audit_path = Path(plan['arguments']['tree-audit']) / 'receipt.json'
    if (handoff['status'] != 'complete_full_local_codon_tree_audit_handoff'
            or handoff['audited_cases'] != plan['cases']
            or handoff['audit_receipt_sha256'] != sha(audit_path)):
        raise ValueError('Complete matching tree audit required')
    command = [sys.executable, 'scripts/run_local_mg94_diagnostics.py']
    for key, value in plan['arguments'].items():
        command.extend(['--' + key, value])
    with (control / 'fits.log').open('x') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    verify()
    output = Path(plan['arguments']['output'])
    receipt = json.loads((output / 'receipt.json').read_text())
    if receipt['status'] != 'complete_full_genus_mg94_execution_pending_audit' or receipt['completed_cases'] != plan['cases']:
        raise ValueError('Incomplete MG94 refit grid')
    result = dict(status='complete_local_mg94_execution_handoff_pending_fit_audit',
                  plan_sha256=digest, tree_audit_handoff_sha256=sha(handoff_path),
                  fit_receipt_sha256=sha(output / 'receipt.json'), cases=plan['cases'], command=command,
                  scope='All eligible local codon diagnostic fits executed on audited fixed topologies. Saved-likelihood replay, parameter normalization, identifiability and selection eligibility remain downstream. No selection test.')
    (control / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
