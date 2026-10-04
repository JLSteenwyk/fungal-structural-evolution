#!/usr/bin/env python3
"""Launch the complete startup-only V6 grid and original-journal dependencies."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time

import psutil

from ancestral_chain_attempt import sha, write_json
from reference_measurement_union_sources import verify
from run_baliphy_scalar_v6_preflight_v1 import SUMMARY_FIELDS


def create(path, value):
    with Path(path).open('x') as f: json.dump(value, f, indent=2); f.write('\n')


def launch(label, command, dependencies, source):
    prefix = 'metadata/baliphy_scalar_v6_preflight_' + label + '_20261004_v1'
    unit = 'fungal-scalar-v6-preflight-' + label + '-20261004-v1.service'
    wait_path = prefix + '_wait_plan.json'; launch_path = prefix + '_launch.json'
    assert not Path(launch_path).exists()
    bindings = {path: sha(path) for path in [source, *dependencies, str(Path(__file__)),
        'scripts/run_after_verified_dependencies_v2.py']}
    create(wait_path, dict(dependencies=dependencies, command=command, pins=bindings))
    cmd = ['systemd-run', '--user', '--collect', '--unit=' + unit,
        '--working-directory=' + str(Path.cwd()), '-p', 'CPUQuota=200%',
        '-p', 'MemoryMax=32G', '-p', 'MemorySwapMax=0', '-p', 'TasksMax=128',
        '--setenv=PYTHONUNBUFFERED=1', '--setenv=OPENBLAS_NUM_THREADS=1',
        '--setenv=OMP_NUM_THREADS=1', '--setenv=MKL_NUM_THREADS=1',
        sys.executable, 'scripts/run_after_verified_dependencies_v2.py', '--plan', wait_path]
    subprocess.run(cmd, check=True)
    pid = 0
    for _ in range(20):
        pid = int(subprocess.check_output(['systemctl', '--user', 'show', unit,
            '-p', 'MainPID', '--value'], text=True))
        if pid: break
        time.sleep(.1)
    assert pid > 0
    proc = psutil.Process(pid)
    expected = [sys.executable, 'scripts/run_after_verified_dependencies_v2.py', '--plan', wait_path]
    assert proc.cmdline() == expected
    cg = subprocess.check_output(['systemctl', '--user', 'show', unit,
        '-p', 'ControlGroup', '--value'], text=True).strip()
    limits = {k: (Path('/sys/fs/cgroup') / cg.lstrip('/') / k).read_text().strip()
        for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert limits == {'cpu.max': '200000 100000', 'memory.max': str(32 * 2**30), 'memory.swap.max': '0'}
    create(launch_path, dict(unit=unit, pid=pid, created=proc.create_time(), cmdline=expected,
        plan=wait_path, plan_sha256=sha(wait_path), checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_cgroup_limits=limits, systemd_launch_command=cmd))
    print(json.dumps(dict(launch=launch_path, pid=pid, unit=unit)), flush=True)
    return launch_path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--validation', type=Path, required=True)
    p.add_argument('--transport', type=Path, required=True)
    p.add_argument('--validate-only', action='store_true')
    a = p.parse_args(); self_hash = sha(__file__)
    plan = json.loads(a.plan.read_text()); verify(plan['pins'])
    gate = json.loads(a.validation.read_text()); transport = json.loads(a.transport.read_text())
    assert gate['status'] == 'passed_full_scalar_v6_startup_software_contracts_v1'
    assert (gate['full_roles'], gate['full_quartets'], gate['native_startups']) == (1620, 405, 3)
    assert gate['full_source_grid_checked'] and gate['synthetic_full_producer_reader_serialization_checked']
    assert len(gate['altered_designs_rejected']) == 12
    assert len(gate['privately_rehashed_invalid_output_cases']) == 7
    assert len(gate['custody_cases_rejected']) == 4 and len(gate['serialization_rejections']) == 6
    assert gate['artificial_failures_retained'] == gate['artificial_unresolved_quartets_retained'] == 2
    assert transport['status'] == 'verified_original_scalar_v6_startup_software_wait_zero'
    assert transport['actual_tool_terminal_exit_code'] == 0 and transport['entire_terminal_payload_matched']
    assert transport['validation_sha256'] == sha(a.validation)
    assert transport['source_hashes'][str(a.validation)] == sha(a.validation)
    verify(transport['source_hashes'])
    for path in [a.validation, a.transport]: assert plan['pins'][str(path)] == sha(path)
    checked_jobs = [path for path in gate['source_hashes'] if path.endswith('/actual_full_grid_jobs.json')]
    assert len(checked_jobs) == 1 and gate['source_hashes'][checked_jobs[0]] == sha(plan['jobs'])
    resources = plan['resources']
    assert (resources['workers'], resources['cpus'], resources['memory_gib'], resources['swap_gib'],
            resources['blas_threads']) == (2, 2, 32, 0, 1)
    assert resources['native_startup_only'] and resources['posterior_sampling'] is False
    assert not plan['dependencies']
    assert psutil.virtual_memory().available >= resources['minimum_available_ram_gib'] * 2**30
    assert psutil.disk_usage('.').free >= resources['minimum_free_disk_gib'] * 2**30
    root = Path(plan['output']); assert not root.exists()
    if a.validate_only:
        print('V6 full startup source/job/software/resource gates passed; no launch.'); return
    producer = launch('producer', [sys.executable, 'scripts/run_baliphy_scalar_v6_preflight_v1.py',
        '--plan', str(a.plan)], [], str(a.plan))
    reader = launch('reader', [sys.executable, 'scripts/readback_baliphy_scalar_v6_preflight_v1.py',
        '--plan', str(a.plan)], [producer], str(a.plan))
    completion_path = 'metadata/baliphy_scalar_v6_preflight_completion_plan_20261004_v1.json'
    completion = dict(source_plan=str(a.plan), producer_receipt=str(root / 'receipt.json'),
        independent_readback=str(root / 'readback.json'),
        producer_status='complete_full_scalar_v6_startup_dispositions_pending_readback_v1',
        reader_status='passed_full_scalar_v6_startup_serialized_readback_v1',
        completed_status='complete_verified_full_scalar_v6_startup_v1', summary_fields=SUMMARY_FIELDS,
        launches=[producer, reader], pins={path: sha(path) for path in [str(a.plan), producer, reader,
            str(Path(__file__)), 'scripts/close_full_triad_sequence_stage.py',
            'scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'], scope=plan['scope'])
    create(completion_path, completion)
    closer = launch('closure', [sys.executable, 'scripts/close_full_triad_sequence_stage.py',
        '--plan', completion_path], [producer, reader], completion_path)
    create('metadata/baliphy_scalar_v6_preflight_launches_20261004_v1.json',
        dict(status='launched_full_scalar_v6_native_startup_only',
        source_plan=str(a.plan), source_plan_sha256=sha(a.plan), launches=[producer, reader, closer],
        launcher_sha256=self_hash, posterior_sampling_launched=False, gpu=False, new_cost_usd=0))
    assert sha(__file__) == self_hash


if __name__ == '__main__':
    main()
