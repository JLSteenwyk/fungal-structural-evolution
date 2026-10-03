#!/usr/bin/env python3
"""Observe exact original proof controllers and verify completed proof bindings."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from reference_measurement_union_sources import verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    plan = json.loads(a.plan.read_text()); verify(plan['pins'])
    inventory = json.loads(Path(plan['launch_inventory']).read_text())
    assert inventory['source_plan_sha256'] == sha(a.plan)
    handles = []
    for path in inventory['launches']:
        row = json.loads(Path(path).read_text()); row['launch'] = path
        assert sha(row['plan']) == row['plan_sha256']
        process = fingerprint(row)
        if process is None:
            observed = journal_terminal(row)
        else:
            observed = live_record(row, process)
            group = next(s[3:] for s in Path(f'/proc/{process.pid}/cgroup').read_text().splitlines() if s.startswith('0::'))
            root = Path('/sys/fs/cgroup') / group.lstrip('/')
            limits = {k: (root / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
            assert limits == row['actual_cgroup_limits'] == {'cpu.max': '200000 100000', 'memory.max': str(16 * 2**30), 'memory.swap.max': '0'}
            observed['current_cgroup_limits'] = limits
        handles.append(observed)
    completion = None
    path = Path(plan['completion'])
    if path.exists():
        completion = json.loads(path.read_text())
        archive = Path(completion['full_hash_archive'])
        assert sha(archive) == completion['full_hash_archive_sha256']
        proof = json.loads(archive.read_text()); verify(proof['source_hashes'])
        assert len(proof['services']) == completion['exact_process_journals_checked'] == 2
        assert len(proof['source_hashes']) == completion['bound_source_hashes']
        assert proof['summary'] == {k: completion[k] for k in proof['summary']}
    result = dict(status='verified_original_full_exact_covariance_proof_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan_sha256=sha(a.plan),
        frozen_pins_checked=len(plan['pins']), original_handles=handles,
        completed_proof=completion, full_closure_bindings_rehashed=(completion['bound_source_hashes'] if completion else 0),
        raw_reml_basis_qualification_complete=False, production_fitting_launched=False,
        gpu=False, all_eight_aims_incomplete=True,
        scope='Original PID/create/command and invocation journals, live caps, immutable pins; every completed proof archive binding rehashed. No restart, source changes, fit, posterior or biological acceptance.')
    with a.output.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['original_handles', 'completed_proof']}))
    if completion: print('closed', completion['status'], completion['retained_basis_counts'])


if __name__ == '__main__': main()
