#!/usr/bin/env python3
"""Observe the exact original process-local stack comparison and its native caps."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import psutil

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import verify


def main():
    pp = Path('metadata/scalar_native_stack_comparison_plan_20261004_v2.json')
    plan = json.loads(pp.read_text())
    verify(plan['pins'])
    unit = 'fungal-scalar-native-stack-comparison-20261004-v2.service'
    pid = int(subprocess.check_output(['systemctl', '--user', 'show', unit, '-p', 'MainPID', '--value'], text=True))
    inv = subprocess.check_output(['systemctl', '--user', 'show', unit, '-p', 'InvocationID', '--value'], text=True).strip()
    initial = json.loads(Path('metadata/scalar_native_stack_comparison_original_tool_initial_20261004_v2.json').read_text())
    assert initial['session_id'] == 99691 and inv in initial['output']
    controller = psutil.Process(pid)
    config = json.loads(Path('metadata/scalar_native_stack_comparison_execution_20261004_v2/configuration.json').read_text())
    assert config['wrapper']['pid'] == pid and config['wrapper']['created'] == controller.create_time()
    assert config['wrapper']['cmdline'] == controller.cmdline() and config['invocation_id'] == inv
    group = next(line[3:] for line in Path('/proc/'+str(pid)+'/cgroup').read_text().splitlines() if line.startswith('0::'))
    cg = Path('/sys/fs/cgroup')/group.lstrip('/')
    caps = {k: (cg/k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert caps == {'cpu.max': '200000 100000', 'memory.max': str(64*2**30), 'memory.swap.max': '0'}
    children = []
    native = []
    for process in controller.children(recursive=True):
        try:
            command = process.cmdline()
            record = dict(pid=process.pid, created=process.create_time(), command=command,
                status=process.status(), rss_bytes=process.memory_info().rss,
                address_space_limits=list(process.rlimit(psutil.RLIMIT_AS)),
                stack_limits=list(process.rlimit(psutil.RLIMIT_STACK)),
                cpu_limits=list(process.rlimit(psutil.RLIMIT_CPU)),
                file_limits=list(process.rlimit(psutil.RLIMIT_FSIZE)))
            if command == plan['native_command']:
                assert record['address_space_limits'] == [48*2**30, 48*2**30]
                assert record['stack_limits'] == [64*2**20, -1]
                assert record['cpu_limits'] == [5674, 5674] and record['file_limits'] == [2*2**30, 2*2**30]
                native.append(record)
            children.append(record)
        except psutil.NoSuchProcess:
            continue
    assert len(native) <= 1
    result = dict(status='verified_live_original_controlled_native_stack_comparison',
        checked_utc=datetime.now(timezone.utc).isoformat(), unit=unit,
        original_tool_session_id=99691, invocation_id=inv, wrapper=config['wrapper'],
        actual_cgroup_limits=caps, children=children, native_child_count=len(native),
        source_plan_sha256=sha(pp), production_stack_limits_changed=False,
        original_attempts_restarted=False, scientific_eligibility=False, gpu=False,
        scope='Exact live original wrapper and any present native target are observed; '
              'native64MiB soft/unlimited hard stack plus original AS/CPU/file limits '
              'are verified directly. Live status does not establish normal completion, '
              'root repair, output integrity or adequate posterior uncertainty.')
    with Path('metadata/scalar_native_stack_comparison_checkpoint_20261004_v2.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'children'}, indent=2))


if __name__ == '__main__':main()
