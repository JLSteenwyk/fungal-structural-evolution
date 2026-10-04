#!/usr/bin/env python3
"""Observe the exact protected debugger and inferior without changing either."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import psutil

from reference_measurement_union_sources import bind, verify


def exact(identity):
    process = psutil.Process(identity['pid'])
    assert process.create_time() == identity['created']
    assert process.cmdline() == identity['cmdline']
    assert process.status() != psutil.STATUS_ZOMBIE
    return process


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); assert not args.output.exists()
    plan_path = Path('metadata/retained_factor_protection_debugger_plan_20261004_v1.json')
    plan = json.loads(plan_path.read_text()); verify(plan['pins'])
    execution_root = Path('metadata/retained_factor_protection_debugger_execution_20261004_v1')
    config_path = execution_root/'configuration.json'; config = json.loads(config_path.read_text())
    wrapper = exact(config['wrapper'])
    debugger_root = Path(plan['debugger_output']); target_root = Path(plan['output'])
    debugger_path = debugger_root/'debugger_process.json'
    inferior_path = target_root/'protected_inferior_identity.json'
    debugger = exact(json.loads(debugger_path.read_text()))
    inferior_identity = json.loads(inferior_path.read_text()); inferior = exact(inferior_identity)
    assert inferior_identity['parent_pid'] == debugger.pid
    assert inferior_identity['parent_command'] == debugger.cmdline()
    assert all(value == '1' for value in inferior_identity['blas_environment'].values())
    group = next(row[3:] for row in Path(f'/proc/{inferior.pid}/cgroup').read_text().splitlines()
                 if row.startswith('0::'))
    cg = Path('/sys/fs/cgroup')/group.lstrip('/')
    caps = {name:(cg/name).read_text().strip()
            for name in ['cpu.max','memory.max','memory.swap.max']}
    assert caps == config['actual_cgroup_limits']
    assert caps == {'cpu.max':'200000 100000','memory.max':str(32*2**30),'memory.swap.max':'0'}
    ranges_path = target_root/'factor_protection_ranges.json'
    ranges = json.loads(ranges_path.read_text()); assert len(ranges) == 5
    maps = Path(f'/proc/{inferior.pid}/maps').read_text()
    mapped = []
    for line in maps.splitlines():
        interval, permission, *_ = line.split()
        start, end = [int(value,16) for value in interval.split('-')]
        mapped.append((start,end,permission))
    for interval in ranges.values():
        assert any(start <= interval['start'] < interval['end_exclusive'] <= end and permission == 'r--p'
                   for start,end,permission in mapped)
    pins = dict(plan['pins'])
    for path in [Path(__file__),plan_path,config_path,debugger_path,inferior_path,ranges_path]: bind(pins,path)
    verify(pins)
    result = dict(status='verified_live_original_protected_debugger_and_inferior',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_tool_session_id=40959,
        invocation_id=config['invocation_id'],wrapper_pid=wrapper.pid,debugger_pid=debugger.pid,
        inferior_pid=inferior.pid,inferior_created=inferior.create_time(),
        inferior_cpu_seconds=sum(inferior.cpu_times()[:2]),inferior_rss_bytes=inferior.memory_info().rss,
        actual_cgroup_limits=caps,protected_factor_ranges_verified=len(ranges),
        last_progress_lines=(debugger_root/'gdb_stdout.log').read_text().splitlines()[-5:],
        receipt_present=Path('metadata/retained_factor_protection_debugger_20261004_v1.json').exists(),
        source_hashes=pins,original_jobs_restarted=False,scientific_eligibility=False,gpu=False,
        scope='Read-only exact live wrapper/GDB/inferior identities and actual caps; all five factor '
              'interior ranges are mapped read-only in the actual inferior. No terminal proof, '
              'fault attribution, source repair, original restart or biological acceptance.')
    with args.output.open('x') as handle:
        json.dump(result,handle,indent=2); handle.write('\n')
    print(json.dumps({key:value for key,value in result.items() if key!='source_hashes'},indent=2))


if __name__ == '__main__':
    main()
