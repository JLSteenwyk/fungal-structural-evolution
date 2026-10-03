#!/usr/bin/env python3
"""Observe the original alternate-reader/retained-basis handles without replay."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from reference_measurement_union_sources import verify


PLANS = {
    'parallel': 'metadata/parallel_covariance_readback_plan_20261003_v1.json',
    'retained': 'metadata/full_reduced_covariance_qualification_plan_20261003_v2.json',
}


def observe(path, first_expected=None):
    record = json.loads(Path(path).read_text())
    record['launch'] = path
    assert sha(record['plan']) == record['plan_sha256']
    process = fingerprint(record)
    if process is None:
        return journal_terminal(record)
    result = live_record(record, process)
    group = next(line[3:] for line in Path(f'/proc/{process.pid}/cgroup').read_text().splitlines()
                 if line.startswith('0::'))
    root = Path('/sys/fs/cgroup') / group.lstrip('/')
    limits = {k: (root / k).read_text().strip()
              for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert limits == record['actual_cgroup_limits']
    if first_expected is not None:
        assert limits == first_expected
    result['current_cgroup_limits'] = limits
    result['current_cgroup_memory_bytes'] = {
        k: int((root / k).read_text()) for k in ['memory.current', 'memory.peak', 'memory.swap.current']}
    result['current_cgroup_memory_events'] = {
        k: int(v) for k, v in (line.split() for line in (root / 'memory.events').read_text().splitlines())}
    # Sample descendants separately; controller memory is not a native peak.
    result['child_resource_observations'] = []
    for child in process.children(recursive=True):
        import psutil
        try:
            result['child_resource_observations'].append(dict(
                pid=child.pid, created=child.create_time(), cmdline=child.cmdline(),
                rss_bytes=child.memory_info().rss,
                process_limits=Path(f'/proc/{child.pid}/limits').read_text()))
        except (psutil.NoSuchProcess, psutil.ZombieProcess, FileNotFoundError):
            result['child_resource_observations'].append(dict(
                pid=child.pid, status='descendant_disappeared_during_observation'))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=PLANS, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    path = Path(PLANS[args.stage]); plan = json.loads(path.read_text())
    verify(plan['pins'])
    inventory = json.loads(Path(plan['launch_inventory']).read_text())
    assert inventory['source_plan_sha256'] == sha(path)
    own = inventory['launches']
    dependencies = ([plan['source_producer_launch']] if args.stage == 'parallel'
                    else plan['dependencies'])
    expected = {'cpu.max': str(plan['resources']['cpus'] * 100000) + ' 100000',
                'memory.max': str(plan['resources']['memory_gib'] * 2**30),
                'memory.swap.max': '0'}
    handles = [observe(fp, expected if i == 0 else None)
               for i, fp in enumerate(own + dependencies)]
    gates = {}
    for key in ['qualification_completion', 'exact_completion', 'completion']:
        if key not in plan:
            continue
        fp = Path(plan[key]); entry = dict(path=str(fp), present=fp.exists())
        if fp.exists():
            value = json.loads(fp.read_text()); archive = Path(value['full_hash_archive'])
            assert sha(archive) == value['full_hash_archive_sha256']
            entry.update(status=value['status'], sha256=sha(fp),
                         declared_source_bindings=value['bound_source_hashes'],
                         original_journals=value['exact_process_journals_checked'],
                         full_hash_archive_sha256=value['full_hash_archive_sha256'])
        gates[key] = entry
    root = Path(plan['output'])
    result = dict(status='verified_original_alternate_covariance_pipeline_runtime',
                  checked_utc=datetime.now(timezone.utc).isoformat(), stage=args.stage,
                  plan=str(path), plan_sha256=sha(path), frozen_pins_checked=len(plan['pins']),
                  original_handles=handles, completion_gates=gates,
                  output_root_present=root.exists(),
                  completed_cohort_files=len(list((root / 'cohorts').glob('*.json'))),
                  producer_receipt_present=(root / 'receipt.json').exists(),
                  independent_readback_present=(root / 'readback.json').exists(),
                  expected=inventory.get('expected', dict(cohorts=4340, audit_rows=1302000,
                                                        setting_audit_links=6220800)),
                  production_fitting_launched=False, gpu=False, all_eight_aims_incomplete=True,
                  scope='Original PID/create/command or invocation-linked terminal journals; live cgroup '
                  'limits/memory and immutable pins checked. Compact closure/archive hashes and available '
                  'file counts only: no full inherited source replay, independent artifact validation, '
                  'precise final native peaks, numerical qualification or biological acceptance. A live '
                  'source-checking child may precede creation of its fresh output root. No restart or mutation.')
    with args.output.open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'original_handles'}, indent=2))


if __name__ == '__main__':
    main()
