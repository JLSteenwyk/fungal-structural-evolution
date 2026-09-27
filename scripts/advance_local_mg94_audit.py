#!/usr/bin/env python3
"""Audit every local MG94 fit after the identified full refit controller exits."""
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
    root = Path(plan['controller_output'])
    root.mkdir(parents=True, exist_ok=True)
    lock = (root / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def verify():
        if sha(args.plan) != digest:
            raise ValueError('Plan changed')
        for name, expected in plan['pins'].items():
            if sha(name) != expected:
                raise ValueError('Pinned file changed: ' + name)

    verify()
    identity = plan['producer']
    while True:
        try:
            process = psutil.Process(identity['pid'])
            if process.create_time() != identity['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            if process.cmdline() != identity['cmdline']:
                raise ValueError('Producer command changed')
        except psutil.NoSuchProcess:
            break
        print('waiting_for_full_local_mg94_execution', identity['pid'], flush=True)
        time.sleep(30)
    verify()
    handoff_path = Path(plan['producer_handoff'])
    handoff = json.loads(handoff_path.read_text())
    fits = Path(plan['fits'])
    batch = json.loads((fits / 'receipt.json').read_text())
    if (handoff['status'] != 'complete_local_mg94_execution_handoff_pending_fit_audit'
            or handoff['cases'] != plan['cases']
            or handoff['fit_receipt_sha256'] != sha(fits / 'receipt.json')
            or batch['status'] != 'complete_full_genus_mg94_execution_pending_audit'
            or batch['completed_cases'] != plan['cases']
            or batch['config_sha256'] != sha(fits / 'config.json')):
        raise ValueError('Full fit execution has not completed')
    expected = {r['case_id']: r['receipt_sha256'] for r in batch['case_receipts']}
    if len(expected) != len(batch['case_receipts']) or len(expected) != plan['cases']:
        raise ValueError('Invalid batch case grid')
    actual = {p.parent.name: sha(p) for p in fits.glob('*/receipt.json')}
    if expected != actual:
        raise ValueError('Full fit receipt grid changed')
    command = [sys.executable, 'scripts/audit_genus_mg94_fits.py', '--fits', str(fits),
               '--install', plan['install'], '--output', plan['audit_output']]
    with (root / 'audit.log').open('x') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    verify()
    audit_root = Path(plan['audit_output'])
    audit = json.loads((audit_root / 'receipt.json').read_text())
    if (audit['status'] not in ['passed_saved_fit_readback', 'complete_saved_fit_readback_with_review_flags']
            or audit['scope'] != 'all_executed_cases' or audit['audited_cases'] != plan['cases']
            or audit['planned_cases'] != plan['cases']
            or audit['source_config_sha256'] != sha(fits / 'config.json')
            or {r['case_id']: r['source_receipt_sha256'] for r in audit['source_cases']} != expected):
        raise ValueError('Full fit audit grid differs')
    for name, expected_hash in audit['artifacts'].items():
        if sha(audit_root / name) != expected_hash:
            raise ValueError('Changed audit artifact: ' + name)
    for row in audit['source_cases']:
        folder = audit_root / row['case_id']
        if (sha(folder / 'reevaluate.bf') != row['reevaluation_script_sha256']
                or sha(folder / 'reevaluate.log') != row['reevaluation_log_sha256']):
            raise ValueError('Changed likelihood replay artifact')
    result = dict(status='complete_full_local_mg94_audit_handoff',
                  fit_audit_status=audit['status'], plan_sha256=digest,
                  source_handoff_sha256=sha(handoff_path), audit_receipt_sha256=sha(audit_root / 'receipt.json'),
                  audited_cases=audit['audited_cases'], numerical_review_cases=audit['cases_with_numerical_review_flags'],
                  maximum_absolute_likelihood_difference=audit['maximum_absolute_likelihood_difference'], command=command,
                  scope='Full saved-likelihood replay and report consistency audit completed. Any numerical flags retained explicitly, never converted to a clean pass. Not evidence of optimality, model adequacy, normalized dS/dN or selection eligibility.')
    (root / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
