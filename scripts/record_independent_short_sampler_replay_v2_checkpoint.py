#!/usr/bin/env python3
"""Observe exact original replay and sampler handles without repeating inference."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from run_baliphy_reference_preflight import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); path = Path('metadata/independent_short_sampler_replay_v2_plan_20261003.json')
    plan = json.loads(path.read_text()); verify(plan)
    inventory = json.loads(Path(plan['launch_inventory']).read_text())
    assert inventory['source_plan_sha256'] == sha(path)
    native = json.loads(Path(plan['sampler_plan']).read_text())
    native_inventory = json.loads(Path(native['launch_inventory']).read_text())
    handles = []
    for i, launch in enumerate(inventory['launches'] + native_inventory['launches']):
        record = json.loads(Path(launch).read_text()); record['launch'] = launch
        assert sha(record['plan']) == record['plan_sha256']
        process = fingerprint(record)
        if process is None: observed = journal_terminal(record)
        else:
            observed = live_record(record, process)
            group = next(x[3:] for x in (Path('/proc')/str(record['pid'])/'cgroup').read_text().splitlines() if x.startswith('0::'))
            cg = Path('/sys/fs/cgroup')/group.lstrip('/')
            caps = {k:(cg/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
            cpu, memory = (2,16) if i < 3 else ((16,200) if i == 3 else (2,32))
            assert caps == record['actual_cgroup_limits'] == {'cpu.max':str(cpu*100000)+' 100000',
                'memory.max':str(memory*2**30),'memory.swap.max':'0'}
            observed['current_cgroup_limits'] = caps
        handles.append(observed)
    checkpoints = [json.loads(p.read_text()) for p in (Path(plan['output'])/'chains').glob('*.json')]
    assert all(r['scientific_eligibility'] is r['posterior_qualified'] is r['ancestral_categories_available'] is r['tip_logs_joint_trajectory_available'] is False for r in checkpoints)
    original = [json.loads(p.read_text()) for p in (Path(native['output'])/'chains').glob('*.json')]
    completed = None
    if Path(plan['completion']).exists():
        completed = json.loads(Path(plan['completion']).read_text())
        assert sha(completed['full_hash_archive']) == completed['full_hash_archive_sha256']
    result = dict(status='verified_original_full_independent_short_sampler_replay_v2_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan_sha256=sha(path), frozen_pins_checked=len(plan['pins']),
        original_handles=handles, unclosed_replay_checkpoints=len(checkpoints),
        replay_status_counts=dict(Counter(r['status'] for r in checkpoints)),
        original_sampler_checkpoints=len(original),original_sampler_status_counts=dict(Counter(r['status'] for r in original)),
        accounting_closure=completed,ancestral_categories_available=False,posterior_qualified=False,tip_logs_joint_trajectory_available=False,gpu=False,
        retained_unknown_draw_disagreements=sum(sum(f['tip_log_draw_disagreement_positions'] for f in r['frames']) for r in checkpoints),
        scope='Immutable pins, exact original six handles/journals and live caps checked. Checkpoint counts are unclosed progress; compact completion archive hash checked if present. No full artifact/native replay, sampler restart, inferred internal categories or scientific acceptance.')
    with args.output.open('x') as handle: handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['original_handles','accounting_closure']}))


if __name__ == '__main__': main()
