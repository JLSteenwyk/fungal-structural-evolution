#!/usr/bin/env python3
"""Wait for the identified divergence producer, then run the full independent audit."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import psutil


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    digest = sha(args.plan)
    root = Path(plan['controller_output'])
    root.mkdir(parents=True, exist_ok=False)
    lock = (root / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def verify():
        assert sha(args.plan) == digest
        for path, expected in plan['pins'].items():
            assert sha(path) == expected, path

    verify()
    identity = plan['producer']
    while True:
        try:
            process = psutil.Process(identity['pid'])
            if process.create_time() != identity['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            assert process.cmdline() == identity['cmdline'], 'Producer command changed'
        except psutil.NoSuchProcess:
            break
        print('waiting_for_identified_divergence_producer', identity['pid'], flush=True)
        time.sleep(30)
    verify()
    comparison_plan = Path(plan['comparison_plan'])
    source_root = Path(json.loads(comparison_plan.read_text())['output'])
    source = source_root / 'receipt.json'
    receipt = json.loads(source.read_text())
    assert receipt['status'] == 'complete_normalized_codon_alignment_divergence_comparison_pending_readback'
    assert receipt['plan_sha256'] == sha(comparison_plan)
    source_hash = sha(source)
    proof_path = root / 'readback.json'
    command = [sys.executable, 'scripts/readback_codon_alignment_divergence.py', '--plan', str(comparison_plan), '--proof', str(proof_path)]
    with (root / 'readback.log').open('x') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    verify()
    proof = json.loads(proof_path.read_text())
    assert proof['status'] == 'passed_full_normalized_codon_divergence_independent_readback'
    assert proof['source_receipt_sha256'] == source_hash == sha(source)
    assert proof['cases'] == plan['expected_cases'] and proof['dispositions'] == plan['expected_dispositions']
    result = dict(status='complete_full_codon_divergence_readback_handoff', plan_sha256=digest,
                  comparison_receipt_sha256=source_hash, readback_sha256=sha(proof_path), cases=proof['cases'],
                  dispositions=proof['dispositions'], split_rows=proof['split_rows'], pair_rows=proof['pair_rows'], command=command,
                  scope='Complete independent comparison arithmetic check; does not establish model adequacy, optimization uniqueness, or selection eligibility.')
    (root / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
