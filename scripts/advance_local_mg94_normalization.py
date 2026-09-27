#!/usr/bin/env python3
"""Bind and normalize every numerically verified local codon fit consistently."""
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
    root = Path(plan['controller_output'])
    root.mkdir(parents=True, exist_ok=True)
    lock = (root / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def verify():
        if sha(args.plan) != plan_hash:
            raise ValueError('Plan changed')
        for path, digest in plan['pins'].items():
            if sha(path) != digest:
                raise ValueError('Pinned input changed: ' + path)

    verify()
    identity = plan['producer']
    while True:
        try:
            process = psutil.Process(identity['pid'])
            if process.create_time() != identity['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            if process.cmdline() != identity['cmdline']:
                raise ValueError('Audit controller identity changed')
        except psutil.NoSuchProcess:
            break
        print('waiting_for_full_local_fit_audit', identity['pid'], flush=True)
        time.sleep(30)
    verify()
    handoff_path = Path(plan['audit_handoff'])
    handoff = json.loads(handoff_path.read_text())
    audit_path = Path(plan['audit']) / 'receipt.json'
    audit = json.loads(audit_path.read_text())
    if (handoff['status'] != 'complete_full_local_mg94_audit_handoff'
            or handoff['audit_receipt_sha256'] != sha(audit_path)
            or audit['scope'] != 'all_executed_cases' or audit['audited_cases'] != plan['cases']):
        raise ValueError('Complete matching fit audit required')
    if audit['status'] != 'passed_saved_fit_readback' or audit['cases_with_numerical_review_flags']:
        (root / 'numerical_review_required.json').write_text(json.dumps({
            'status': 'normalization_not_started_due_to_numerical_review_flags',
            'audit_receipt_sha256': sha(audit_path), 'review_cases': audit['cases_with_numerical_review_flags']}, indent=2) + '\n')
        raise ValueError('Numerical audit flags require review before normalization')
    if audit['branches'] != plan['branches']:
        raise ValueError('Audited branch grid differs')
    binding = root / 'tree_binding.json'
    command = [sys.executable, 'scripts/verify_genus_mg94_tree_bindings.py',
               '--fits', plan['fits'], '--trees', plan['trees'], '--audit', plan['tree_audit'],
               '--information', plan['information'], '--output', str(binding)]
    with (root / 'binding.log').open('x') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    binding_result = json.loads(binding.read_text())
    if binding_result['status'] != 'passed_full_mg94_to_audited_tree_binding' or binding_result['cases'] != plan['cases']:
        raise ValueError('Full tree binding failed')
    execution = dict(cases=plan['cases'], branches=plan['branches'], workers=2,
                     source_fit_receipt_sha256=sha(Path(plan['fits']) / 'receipt.json'),
                     source_audit_sha256=sha(audit_path), source_definition_sha256=sha(plan['source']),
                     binding_sha256=sha(binding), corrected_helper_sha256=sha(plan['helper']),
                     interpretation=plan['interpretation'])
    execution_path = root / 'normalization_plan.json'
    with execution_path.open('x') as handle:
        json.dump(execution, handle, indent=2); handle.write('\n')
    command = [sys.executable, 'scripts/normalize_genus_mg94_branches.py',
               '--fits', plan['fits'], '--audit', plan['audit'], '--binding', str(binding),
               '--source', plan['source'], '--plan', str(execution_path), '--install', plan['install'],
               '--opportunity-helper', plan['helper'], '--output', plan['output']]
    with (root / 'normalization.log').open('x') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    verify()
    output = Path(plan['output'])
    result = json.loads((output / 'receipt.json').read_text())
    if (result['status'] != 'complete_independently_checked_mg94_opportunity_normalization'
            or result['cases'] != plan['cases'] or result['branches'] != plan['branches']
            or result['per_fit_independent_checks'] != plan['cases']
            or result['sense_codon_opportunity_checks'] != 244):
        raise ValueError('Normalization incomplete')
    for path, digest in result['artifacts'].items():
        if sha(output / path) != digest:
            raise ValueError('Changed normalization artifact')
    for check in result['checks']:
        folder = output / check['case_id']
        if sha(folder / 'check.bf') != check['check_script_sha256'] or sha(folder / 'check.log') != check['check_log_sha256']:
            raise ValueError('Changed independent opportunity check')
    receipt = dict(status='complete_full_local_mg94_normalization_handoff', plan_sha256=plan_hash,
                   binding_receipt_sha256=sha(binding), normalization_receipt_sha256=sha(output / 'receipt.json'),
                   cases=plan['cases'], branches=plan['branches'], scope=plan['interpretation'])
    (root / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
