#!/usr/bin/env python3
"""Launch startup-only full grid after native/model and full-grid software gates."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import psutil

from launch_full_triad_sequence_geometry_followup import create
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from run_baliphy_joint_logger_preflight_v5 import SUMMARY_FIELDS, verify
from ancestral_chain_attempt import sha


def launch(name, command, dependencies, source):
    prefix = 'metadata/' + name.replace('-', '_') + '_v5'; unit = 'fungal-' + name + '-20261003-v5.service'
    wait_path = prefix + '_wait_plan_20261003.json'; launch_path = prefix + '_launch_20261003.json'
    assert not Path(launch_path).exists()
    create(wait_path, dict(dependencies=dependencies, command=command,
        pins={p: sha(p) for p in [source, *dependencies, 'scripts/run_after_verified_dependencies_v2.py']}))
    cmd = ['systemd-run', '--user', '--collect', '--unit=' + unit, '--working-directory=' + str(Path.cwd()),
        '--property=CPUQuota=200%', '--property=MemoryMax=32G', '--property=MemorySwapMax=0',
        '--setenv=PYTHONUNBUFFERED=1', '--setenv=OPENBLAS_NUM_THREADS=1', '--setenv=OMP_NUM_THREADS=1',
        '--setenv=MKL_NUM_THREADS=1', sys.executable, 'scripts/run_after_verified_dependencies_v2.py', '--plan', wait_path]
    subprocess.run(cmd, check=True); pid = 0
    for _ in range(20):
        pid = int(subprocess.check_output(['systemctl', '--user', 'show', unit, '-p', 'MainPID', '--value'], text=True))
        if pid: break
        time.sleep(.1)
    assert pid > 0; proc = psutil.Process(pid)
    expected = [sys.executable, 'scripts/run_after_verified_dependencies_v2.py', '--plan', wait_path]
    assert proc.cmdline() == expected
    cg = subprocess.check_output(['systemctl', '--user', 'show', unit, '-p', 'ControlGroup', '--value'], text=True).strip()
    limits = {k: (Path('/sys/fs/cgroup') / cg.lstrip('/') / k).read_text().strip()
              for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert limits == {'cpu.max': '200000 100000', 'memory.max': str(32 * 2**30), 'memory.swap.max': '0'}
    create(launch_path, dict(unit=unit, pid=pid, created=proc.create_time(), cmdline=expected,
        plan=wait_path, plan_sha256=sha(wait_path), checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_cgroup_limits=limits, systemd_launch_command=cmd))
    print(json.dumps(dict(launch=launch_path, pid=pid, unit=unit)), flush=True); return launch_path


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--grid-validation', type=Path, required=True); a = p.parse_args()
    plan = json.loads(a.plan.read_text()); verify(plan)
    grid = json.loads(a.grid_validation.read_text())
    assert grid['status'] == 'passed_full_joint_logger_preflight_software_contracts_v5'
    assert grid['full_roles'] == 1620 and grid['full_quartets'] == 405
    assert len(grid['altered_designs_rejected'])==11
    assert grid['native_startups']==3 and grid['full_producer_reader_serialization_checked']
    assert grid['artificial_failures_retained']==2
    assert grid['full_source_grid_checked']
    for path, h in grid['source_hashes'].items(): assert sha(path) == h, path
    assert a.grid_validation.as_posix() in plan['pins'], 'Grid software gate must be frozen in the execution plan'
    assert plan['resources']['posterior_sampling'] is False
    assert psutil.virtual_memory().available >= 32 * 2**30 and shutil.disk_usage('.').free >= 64 * 2**30
    assert not Path(plan['output']).exists()
    for path in plan['dependencies']:
        r = json.loads(Path(path).read_text()); r['launch'] = path
        assert fingerprint(r) is None; journal_terminal(r)
    execution_path=Path('metadata/baliphy_joint_logger_preflight_v5_software_execution_20261003.json')
    execution=json.loads(execution_path.read_text())
    assert execution['receipt_sha256']==sha(a.grid_validation)
    assert execution['status']=='exited_zero_with_receipt' and execution['exit_code']==0
    verify({'pins':execution['source_hashes']})
    root = Path(plan['output'])
    producer = launch('baliphy-joint-logger-preflight', [sys.executable, 'scripts/run_baliphy_joint_logger_preflight_v5.py',
        '--plan', str(a.plan)], plan['dependencies'], str(a.plan))
    reader = launch('baliphy-joint-logger-preflight-readback', [sys.executable, 'scripts/readback_baliphy_joint_logger_preflight_v5.py',
        '--plan', str(a.plan)], [producer], str(a.plan))
    completion = dict(source_plan=str(a.plan), producer_receipt=str(root / 'receipt.json'),
        independent_readback=str(root / 'readback.json'),
        producer_status='complete_full_joint_logger_startup_dispositions_pending_readback_v5',
        reader_status='passed_full_joint_logger_startup_serialized_readback_v5',
        completed_status='complete_verified_full_joint_logger_startup_v5', summary_fields=SUMMARY_FIELDS,
        launches=[producer, reader], pins={p: sha(p) for p in [str(a.plan), producer, reader,
        'scripts/close_full_triad_sequence_stage.py', 'scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'], scope=plan['scope'])
    create(plan['completion_plan'], completion)
    closer = launch('baliphy-joint-logger-preflight-closure', [sys.executable, 'scripts/close_full_triad_sequence_stage.py',
        '--plan', plan['completion_plan']], [producer, reader], plan['completion_plan'])
    create(plan['launch_inventory'], dict(status='launched_full_native_reference_startup_only',
        source_plan=str(a.plan), source_plan_sha256=sha(a.plan), launches=[producer, reader, closer],
        posterior_sampling_launched=False, gpu=False, new_cost_usd=0))


if __name__ == '__main__': main()
