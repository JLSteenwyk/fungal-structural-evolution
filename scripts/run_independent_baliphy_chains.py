#!/usr/bin/env python3
"""Run the complete input/prior/seed grid under externally enforced resource caps."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import fcntl
import json
from pathlib import Path
import shutil
import subprocess

from ancestral_chain_attempt import run_attempt, sha, write_json
from readback_independent_baliphy_chain import readback


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    digest = sha(args.plan)
    def verify():
        assert sha(args.plan) == digest
        for path, expected in plan['pins'].items():
            assert sha(path) == expected, path
    verify()
    # Refuse an uncapped direct invocation; systemd cgroup limits cover descendants.
    state = dict(line.split('=', 1) for line in subprocess.check_output([
        'systemctl', '--user', 'show', plan['unit'], '-p', 'MemoryMax',
        '-p', 'MemorySwapMax', '-p', 'CPUQuotaPerSecUSec'], text=True).splitlines())
    assert int(state['MemoryMax']) == plan['memory_bytes'] and state['MemorySwapMax'] == '0'
    assert state['CPUQuotaPerSecUSec'] == str(plan['workers']) + 's', state
    group = next(line.split('::', 1)[1] for line in Path('/proc/self/cgroup').read_text().splitlines() if line.startswith('0::'))
    assert group.endswith('/' + plan['unit']), group
    for unit in plan['completed_units']:
        terminal = dict(line.split('=', 1) for line in subprocess.check_output([
            'systemctl', '--user', 'show', unit, '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        assert terminal == dict(ActiveState='inactive', Result='success', ExecMainStatus='0'), terminal
    jobs = json.loads(Path(plan['jobs']).read_text())
    assert len(jobs) == len({j['chain']['chain_id'] for j in jobs}) == 1620
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=True)
    lock = (out / 'run.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    binding = out / 'run_plan.json'
    if binding.exists():
        assert binding.read_bytes() == args.plan.read_bytes()
    else:
        binding.write_bytes(args.plan.read_bytes())
    def run(job):
        assert shutil.disk_usage(out).free >= plan['minimum_free_disk_bytes']
        chain, config = job['chain'], job['config']
        folder = out / chain['chain_id']
        receipt = run_attempt(folder, config)
        attempt = json.loads(receipt.read_text())
        if attempt['exit_code'] != 0 or attempt['status'] != 'exited_zero_pending_scientific_validation':
            return dict(chain_id=chain['chain_id'], status=attempt['status'], receipt=str(receipt), receipt_sha256=sha(receipt))
        audit_path = folder / (receipt.parent.name + '-sample-audit.json')
        if audit_path.exists():
            audit = json.loads(audit_path.read_text())
            assert audit['attempt_receipt_sha256'] == sha(receipt)
        else:
            audit = readback(chain, receipt, plan['iterations'], plan['mapping'])
            write_json(audit_path, audit)
        return dict(chain_id=chain['chain_id'], status=audit['status'],
                    receipt=str(receipt), receipt_sha256=sha(receipt),
                    sample_audit=str(audit_path), sample_audit_sha256=sha(audit_path))
    completed = []
    iterator = iter(jobs)
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        active = {pool.submit(run, next(iterator)) for _ in range(plan['workers'])}
        while active:
            done, active = wait(active, return_when=FIRST_COMPLETED)
            for future in done:
                result = future.result()
                completed.append(result)
                print(len(completed), '/1620', result['chain_id'], result['status'], flush=True)
                job = next(iterator, None)
                if job is not None:
                    active.add(pool.submit(run, job))
    verify()
    assert len(completed) == 1620
    write_json(out / 'receipt.json', dict(status='all_initial_horizon_dispositions_recorded',
        plan_sha256=digest, chains=sorted(completed, key=lambda r: r['chain_id']),
        counts=dict(Counter(r['status'] for r in completed)),
        scope='Full grid accounted for; failed/capped chains retained. Sampling horizon is not convergence.'))


if __name__ == '__main__':
    main()
