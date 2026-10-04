#!/usr/bin/env python3
"""Observe exact original hardware-watchpoint debugger and writable factor words."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import psutil
from reference_measurement_union_sources import bind, verify


def exact(identity):
    process = psutil.Process(identity['pid'])
    assert process.create_time() == identity['created'] and process.cmdline() == identity['cmdline']
    assert process.status() != psutil.STATUS_ZOMBIE
    return process


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    plan_path = Path('metadata/retained_factor_watchpoint_debugger_plan_20261004_v1.json')
    plan = json.loads(plan_path.read_text()); verify(plan['pins'])
    eroot = Path('metadata/retained_factor_watchpoint_debugger_execution_20261004_v1')
    config_path = eroot/'configuration.json'; config = json.loads(config_path.read_text())
    wrapper = exact(config['wrapper'])
    droot = Path(plan['debugger_output']); root = Path(plan['output'])
    dpath = droot/'debugger_process.json'; ipath = root/'watchpoint_inferior_identity.json'
    debugger = exact(json.loads(dpath.read_text())); inferior = exact(json.loads(ipath.read_text()))
    assert inferior.ppid() == debugger.pid
    assert set([debugger.pid, inferior.pid]) <= {q.pid for q in wrapper.children(recursive=True)}
    group = next(row[3:] for row in Path(f'/proc/{inferior.pid}/cgroup').read_text().splitlines() if row.startswith('0::'))
    cg = Path('/sys/fs/cgroup')/group.lstrip('/')
    caps = {key: (cg/key).read_text().strip() for key in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert caps == config['actual_cgroup_limits'] == {
        'cpu.max': '200000 100000', 'memory.max': str(32*2**30), 'memory.swap.max': '0'}
    environment = inferior.environ()
    assert all(environment[key] == '1' for key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'])
    metadata = Path(plan['watchpoint_metadata']); watch = json.loads(metadata.read_text())
    assert watch['inferior_pid'] == inferior.pid and watch['inferior_created'] == inferior.create_time()
    assert not watch['process_local_read_only_protection'] and not watch['artificial_control']
    mappings = []
    for line in Path(f'/proc/{inferior.pid}/maps').read_text().splitlines():
        interval, permissions, *_ = line.split(); start, end = [int(x, 16) for x in interval.split('-')]
        mappings.append((start, end, permissions))
    assert len(watch['targets']) == 2
    for name, word in watch['targets'].items():
        assert any(start <= word['address'] < word['address']+8 <= end and permissions == 'rw-p'
            for start, end, permissions in mappings)
    pins = dict(plan['pins'])
    for path in [Path(__file__), plan_path, config_path, dpath, ipath, metadata]: bind(pins, path)
    verify(pins)
    lines = (droot/'gdb_stdout.log').read_text().splitlines()
    assert sum(line.startswith('PROJECT_ARMED_WATCH=') for line in lines) == 2
    result = dict(status='verified_live_original_unprotected_factor_watchpoint_debugger',
        checked_utc=datetime.now(timezone.utc).isoformat(), original_tool_session_id=91413,
        invocation_id=config['invocation_id'], wrapper_pid=wrapper.pid, debugger_pid=debugger.pid,
        inferior_pid=inferior.pid, inferior_created=inferior.create_time(),
        inferior_cpu_seconds=sum(inferior.cpu_times()[:2]), inferior_rss_bytes=inferior.memory_info().rss,
        actual_cgroup_limits=caps, watched_words_mapped_writable=2,
        completed_group_files=len(list((root/'groups').glob('*.json'))),
        last_progress_lines=lines[-8:], source_hashes=pins,
        native_write_instruction_identified=False, original_corruption_resolved=False,
        original_jobs_restarted=False, scientific_eligibility=False, gpu=False,
        scope='Exact original live wrapper/GDB/inferior, caps and one BLAS thread observed. '
              'Both watched source words remain mapped writable; watchpoint installation recorded. '
              'No terminal proof, exact native write attribution, repair, retries or biological acceptance.')
    with a.output.open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__': main()
