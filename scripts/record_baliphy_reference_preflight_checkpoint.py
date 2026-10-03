#!/usr/bin/env python3
"""Observe the original full startup preflight without restarting or sampling."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from run_baliphy_reference_preflight import verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, default=Path('metadata/baliphy_reference_preflight_plan_20261003.json'))
    p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    plan = json.loads(a.plan.read_text()); digest = sha(a.plan); verify(plan)
    inventory_path = Path(plan['launch_inventory']); inventory = json.loads(inventory_path.read_text())
    assert inventory['source_plan_sha256'] == digest
    handles = []
    for path in inventory['launches']:
        r = json.loads(Path(path).read_text()); r['launch'] = path
        assert sha(r['plan']) == r['plan_sha256']; proc = fingerprint(r)
        if proc is None: observed = journal_terminal(r)
        else:
            observed = live_record(r, proc)
            try:
                group = next(x[3:] for x in (Path('/proc') / str(r['pid']) / 'cgroup').read_text().splitlines()
                             if x.startswith('0::'))
                cg = Path('/sys/fs/cgroup') / group.lstrip('/')
                limits = {k: (cg / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
                assert limits == r['actual_cgroup_limits'] == {
                    'cpu.max': '200000 100000', 'memory.max': str(24 * 2**30), 'memory.swap.max': '0'}
                observed['current_cgroup_limits'] = limits
                observed['observed_memory_bytes'] = {k: (cg / k).read_text().strip()
                    for k in ['memory.current', 'memory.peak', 'memory.swap.current']}
                native_limits = []
                for child in observed['children']:
                    if not child['cmdline'] or not child['cmdline'][0].endswith('/bali-phy'): continue
                    try:
                        raw = (Path('/proc') / str(child['pid']) / 'limits').read_text()
                        actual = {}
                        for name in ['Max address space', 'Max cpu time', 'Max file size']:
                            line = next(x for x in raw.splitlines() if x.startswith(name))
                            actual[name] = line[len(name):].split()[:2]
                        assert actual == {'Max address space': [str(12 * 2**30)] * 2,
                                          'Max cpu time': ['600'] * 2,
                                          'Max file size': [str(256 * 2**20)] * 2}
                        native_limits.append(dict(pid=child['pid'], limits=actual))
                    except FileNotFoundError:
                        native_limits.append(dict(pid=child['pid'], status='Ephemeral native worker ended during observation'))
                observed['native_worker_limits'] = native_limits
            except FileNotFoundError:
                assert fingerprint(r) is None; observed = journal_terminal(r)
        handles.append(observed)
    root = Path(plan['output']); paths = sorted((root / 'chains').glob('*.json')); counts = Counter()
    for path in paths:
        r = json.loads(path.read_text()); assert r['chain_id'] == path.stem and r['plan_sha256'] == digest
        assert r['scientific_eligibility'] is False; counts[r['status']] += 1
    closure = None
    if Path(plan['completion']).exists():
        closure = json.loads(Path(plan['completion']).read_text())
        assert closure['status'] == 'complete_verified_full_reference_alignment_startup'
        assert sha(closure['full_hash_archive']) == closure['full_hash_archive_sha256']
    result = dict(status='verified_original_full_reference_startup_runtime_checkpoint',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan=str(a.plan), source_plan_sha256=digest,
        frozen_source_pins_checked=len(plan['pins']), observer_sha256=sha(__file__),
        launch_inventory_sha256=sha(inventory_path), original_handles=handles,
        unclosed_startup_checkpoints=len(paths), startup_status_counts=dict(counts), accounting_closure=closure,
        posterior_sampling_launched=False, gpu_inference_paused=True, all_eight_aims_incomplete=True,
        scope='All frozen startup pins, exact original handles or invocation-linked terminal journals, live cgroups and native worker limits observed. Checkpoints are unclosed progress. Compact completion archive hash checked when present; no full native-output replay by this observer. No prior inference restart or scientific/memory/mixing acceptance.')
    with a.output.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['original_handles', 'accounting_closure']}), flush=True)


if __name__ == '__main__': main()
