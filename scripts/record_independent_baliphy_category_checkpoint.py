#!/usr/bin/env python3
"""Observe original categorical handles, caps and unclosed checkpoint progress."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=Path('metadata/independent_baliphy_category_plan_20261002.json'))
    parser.add_argument('--output', type=Path, required=True); args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    pins = {str(Path(__file__)):sha(__file__), str(args.plan):sha(args.plan), **plan['pins']}
    assert all(sha(p) == d for p,d in pins.items())
    inventory_path = Path(plan['launch_inventory']); inventory = json.loads(inventory_path.read_text())
    assert inventory['source_plan_sha256'] == sha(args.plan)
    handles = []
    for path in inventory['launches']:
        record = json.loads(Path(path).read_text()); record['launch'] = path
        assert sha(record['plan']) == record['plan_sha256']
        proc = fingerprint(record)
        if proc is None:
            observation = journal_terminal(record)
        else:
            observation = live_record(record, proc)
            lines = (Path('/proc') / str(record['pid']) / 'cgroup').read_text().splitlines()
            unified = [line[3:] for line in lines if line.startswith('0::')]; assert len(unified) == 1
            root = Path('/sys/fs/cgroup') / unified[0].lstrip('/')
            limits = {k:(root / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
            assert limits == record['actual_cgroup_limits'] == {
                'cpu.max':'200000 100000', 'memory.max':str(32 * 2**30), 'memory.swap.max':'0'}
            observation['current_cgroup_limits'] = limits
            observation['observed_memory_bytes'] = {k:(root / k).read_text().strip()
                for k in ['memory.current', 'memory.peak', 'memory.swap.current']}
        handles.append(observation)
    root = Path(plan['output']); checkpoints = sorted((root / 'groups').glob('*.json'))
    patterns = indicators = complete = failed = 0
    for path in checkpoints:
        row = json.loads(path.read_text())
        assert row['group'] == path.stem and row['stage']['plan_sha256'] == sha(args.plan)
        result = row['result']; assert result['scientific_eligibility'] is False
        if result['status'] == 'unresolved_failed_native_chain_retained':
            assert result['cutoffs'] == {}; failed += 1
        else:
            assert result['status'] == 'complete_categorical_comparison_not_posterior_qualification'
            assert set(result['cutoffs']) == {'250', '500'}; complete += 1
            patterns += sum(r['patterns'] for r in result['cutoffs'].values())
            indicators += sum(r['declared_indicator_rows'] for r in result['cutoffs'].values())
    closure = None; cp = Path(plan['completion'])
    if cp.exists():
        closure = json.loads(cp.read_text())
        assert closure['status'] == 'complete_verified_full_independent_baliphy_categorical_comparison'
        assert sha(closure['full_hash_archive']) == closure['full_hash_archive_sha256']
    result = dict(status='verified_original_categorical_runtime_and_unclosed_checkpoint_observation',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan=str(args.plan), source_hashes=pins,
        original_handles=handles, expected=plan['expected'],
        unclosed_producer_checkpoints=dict(groups=len(checkpoints), complete_groups=complete,
            failed_groups=failed, pattern_cutoff_rows=patterns, declared_indicator_rows=indicators),
        accounting_closure=closure, gpu_inference_paused=True, all_eight_aims_incomplete=True,
        scope='Frozen small pins and exact original PID/create/command or actual completion journals observed; '
              'live cgroups and memory observed without restarting jobs. Checkpoint counts are producer '
              'progress, not full numerical/artifact/source/journal closure or scientific acceptance. '
              'A compact closure, when present, has its archive hash rechecked; this observer does not '
              'repeat full large-data hashes or native parsing.')
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes', 'original_handles', 'accounting_closure']}))


if __name__ == '__main__': main()
