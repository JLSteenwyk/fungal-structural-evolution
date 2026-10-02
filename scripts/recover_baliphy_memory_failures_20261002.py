#!/usr/bin/env python3
"""Run separate same-seed attempts for all failed initial-horizon chains."""
import argparse
from collections import Counter
import copy
import fcntl
import json
from pathlib import Path
import shutil
import subprocess
from ancestral_chain_attempt import run_attempt, check_inputs, live_group
from readback_independent_baliphy_chain import readback
from record_project_runtime_checkpoint_v4 import journal_terminal
from run_ortholog_pair_guide_comparison import sha


def selected_jobs(plan):
    source_plan = json.loads(Path(plan['producer_plan']).read_text())
    source = Path(source_plan['output']) / 'receipt.json'; receipt = json.loads(source.read_text())
    assert receipt['status'] == 'all_initial_horizon_dispositions_recorded' and receipt['plan_sha256'] == sha(plan['producer_plan'])
    jobs = json.loads(Path(source_plan['jobs']).read_text()); index = {j['chain']['chain_id']: j for j in jobs}
    dispositions = {r['chain_id']: r for r in receipt['chains']}
    assert len(jobs) == len(index) == len(dispositions) == len(receipt['chains']) == 1620 and set(index) == set(dispositions)
    failed = {c for c, r in dispositions.items() if r['status'] != 'all_saved_alignments_and_candidate_nodes_checked'}
    assert failed == set(plan['failed_chain_ids']) and len(failed) == 3
    selected = []
    for cid in sorted(failed):
        old = index[cid]; record = dispositions[cid]; path = Path(record['receipt'])
        assert sha(path) == record['receipt_sha256']; attempt = json.loads(path.read_text())
        assert attempt['status'] == record['status'] == 'failed' and attempt['exit_code'] != 0
        for name, digest in attempt['artifacts'].items(): assert sha(path.parent / name) == digest
        assert 'std::bad_alloc' in (path.parent / 'stderr.log').read_text(errors='replace')
        process = json.loads((path.parent / 'process.json').read_text()); assert not live_group(process['pgid'])
        config = copy.deepcopy(old['config']); assert config['command'][0] == '/usr/bin/prlimit'
        flags = [i for i, x in enumerate(config['command']) if x.startswith('--as=')]; assert len(flags) == 1
        position = flags[0]; assert config['command'][position] == '--as=12884901888'
        config['command'][position] = '--as=' + str(plan['native_address_space_bytes'])
        # Only the computational limit changes. Priors, inputs, executable,
        # iteration horizon, seed, output flags and timeout retain original values.
        expected = copy.deepcopy(old['config']); expected['command'][position] = config['command'][position]
        assert config == expected and config['seed'] == old['chain']['seed']
        check_inputs(config)
        selected.append(dict(chain=old['chain'], config=config, original=record))
    return selected, source_plan


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); a = p.parse_args()
    plan = json.loads(a.plan.read_text()); digest = sha(a.plan)
    def verify():
        assert sha(a.plan) == digest
        for path, expected in plan['pins'].items(): assert sha(path) == expected, path
    verify()
    state = dict(line.split('=', 1) for line in subprocess.check_output(['systemctl', '--user', 'show', plan['unit'],
        '-p', 'MemoryMax', '-p', 'MemorySwapMax', '-p', 'CPUQuotaPerSecUSec'], text=True).splitlines())
    assert state == dict(MemoryMax=str(64 * 2**30), MemorySwapMax='0', CPUQuotaPerSecUSec='1s')
    cgroup = next(x.split('::', 1)[1] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
    assert cgroup.endswith('/' + plan['unit'])
    for launch_path in plan['original_launches']:
        r = json.loads(Path(launch_path).read_text()); r['launch'] = launch_path; journal_terminal(r)
    jobs, source = selected_jobs(plan); root = Path(plan['output']); root.mkdir(exist_ok=True)
    lock = (root / 'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (root / 'receipt.json').exists(), 'Completed recovery cannot be restarted'
    marker = root / 'stage_plan.json'; frozen = dict(plan_sha256=digest)
    if marker.exists(): assert json.loads(marker.read_text()) == frozen
    else: marker.write_text(json.dumps(frozen, indent=2) + '\n')
    results = []
    for job in jobs:
        verify(); assert shutil.disk_usage(root).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
        cid = job['chain']['chain_id']; rp = run_attempt(root / cid, job['config']); attempt = json.loads(rp.read_text())
        result = dict(chain_id=cid, original_failed_attempt=job['original']['receipt'],
            original_failed_attempt_sha256=job['original']['receipt_sha256'], new_attempt=str(rp), new_attempt_sha256=sha(rp),
            status=attempt['status'], native_bad_alloc='std::bad_alloc' in (rp.parent / 'stderr.log').read_text(errors='replace'))
        if attempt['exit_code'] == 0 and attempt['status'] == 'exited_zero_pending_scientific_validation':
            audit = readback(job['chain'], rp, source['iterations'], source['mapping'])
            ap = rp.parent.parent / (rp.parent.name + '-sample-audit.json')
            with ap.open('x') as handle: handle.write(json.dumps(audit, indent=2) + '\n')
            result.update(status=audit['status'], sample_audit=str(ap), sample_audit_sha256=sha(ap))
        results.append(result); print(json.dumps(result), flush=True)
    verify()
    receipt = dict(status='complete_baliphy_memory_recovery_pending_full_readback', plan_sha256=digest,
        full_original_grid_chains=1620, original_checked_chains=1617, recovery_chains=results,
        recovery_status_counts=dict(Counter(r['status'] for r in results)), original_failed_samples_concatenated=False,
        original_attempts_overwritten=False, scientific_eligibility=False, scope=plan['scope'])
    with (root / 'receipt.json').open('x') as handle: handle.write(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__': main()
