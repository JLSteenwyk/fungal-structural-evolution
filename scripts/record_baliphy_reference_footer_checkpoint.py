#!/usr/bin/env python3
"""Observe original follow-up handles and replay every currently saved startup."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_reference_startup_readback_v2 import inspect_v2
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from run_baliphy_reference_preflight import inspect, verify


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); path = Path('metadata/baliphy_reference_startup_footer_plan_20261003.json')
    plan = json.loads(path.read_text()); verify(plan); handles = []
    inventory = json.loads(Path('metadata/baliphy_reference_startup_footer_launches_20261003.json').read_text())
    assert inventory['source_plan_sha256'] == sha(path)
    for lp in inventory['launches']:
        r = json.loads(Path(lp).read_text()); r['launch'] = lp
        assert sha(r['plan']) == r['plan_sha256']; process = fingerprint(r)
        if process is None: observed = journal_terminal(r)
        else:
            observed = live_record(r, process)
            group = next(x[3:] for x in (Path('/proc') / str(r['pid']) / 'cgroup').read_text().splitlines() if x.startswith('0::'))
            root = Path('/sys/fs/cgroup') / group.lstrip('/')
            limits = {k: (root / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
            assert limits == r['actual_cgroup_limits'] == {
                'cpu.max': '200000 100000', 'memory.max': str(24 * 2**30), 'memory.swap.max': '0'}
            observed['current_cgroup_limits'] = limits
        handles.append(observed)
    native_path = Path(plan['native_plan']); native = json.loads(native_path.read_text()); verify(native)
    jobs = {j['chain']['chain_id']: j for j in json.loads(Path(native['jobs']).read_text())}
    paths = sorted((Path(native['output']) / 'chains').glob('*.json')); rows = []; bindings = {str(path): sha(path)}
    for cp in paths:
        old = json.loads(cp.read_text()); job = jobs[old['chain_id']]; receipt = Path(old['attempt_receipt'])
        assert old == inspect(job, receipt, sha(native_path))
        row = inspect_v2(job, receipt, sha(native_path), sha(path)); rows.append(row)
        bindings[str(cp)] = sha(cp); bindings[str(receipt)] = sha(receipt)
        for name, h in json.loads(receipt.read_text())['artifacts'].items(): bindings[str(receipt.parent / name)] = h
    result = dict(status='verified_original_footer_queue_and_all_available_startup_replay',
        checked_utc=datetime.now(timezone.utc).isoformat(), frozen_pins_checked=len(plan['pins']),
        observer_sha256=sha(__file__), original_handles=handles, available_original_startups=len(rows),
        revised_startup_status_counts=dict(Counter(r['status'] for r in rows)),
        reclassified_v1_dispositions=sum(r['status'] != r['original_startup_disposition'] for r in rows),
        checked_native_timing_footers=sum(r.get('stdout_timing_summary') is not None for r in rows),
        unresolved_available_startup_records=[r for r in rows if r['status'] != 'reference_startup_homology_density_and_representation_checked'],
        actual_checkpoint_and_native_artifact_hashes=bindings,
        full_footer_replay_completion_present=Path(plan['completion']).exists(),
        posterior_sampling_launched=False, scientific_eligibility=False,
        scope='Every original startup checkpoint available in this snapshot reread with both preserved and footer-aware decoders; actual input/runtime homology, representation, density and receipt/configuration/process/artifact links checked. Original follow-up PID/create/command and live caps observed. Full1620readback/source/two-journal closure is pending. No new native runs, allocation repair, converged posterior or project completion.')
    with a.output.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['original_handles', 'actual_checkpoint_and_native_artifact_hashes', 'unresolved_available_startup_records']}), flush=True)


if __name__ == '__main__': main()
