#!/usr/bin/env python3
"""Queue independent missing-PAE readback after the exact downloader completes."""
import json
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
import time

import psutil

from run_after_verified_dependencies import live
from run_ortholog_pair_guide_comparison import sha


UNIT = 'fungal-full-atlas-missing-pae-readback-20261006-v1.service'
WAIT = Path('metadata/full_atlas_missing_pae_readback_wait_plan_20261006_v1.json')
LAUNCH = Path('metadata/full_atlas_missing_pae_readback_launch_20261006_v1.json')
PRODUCER_LAUNCH = Path('metadata/full_atlas_missing_pae_launch_20261005_v2.json')
PRODUCER_PLAN = Path('metadata/full_atlas_missing_pae_plan_20261005_v2.json')
PRODUCER_RECEIPT = Path('metadata/full_atlas_missing_pae_20261005_v2.json')
RESOURCES = Path('metadata/full_atlas_missing_pae_readback_resources_20261006_v1.json')
EXECUTION = Path('metadata/full_atlas_missing_pae_readback_execution_20261006_v1.json')
RECEIPT = Path('metadata/full_atlas_missing_pae_readback_20261006_v1.json')
OUTPUT = Path('results/structures/full-atlas-missing-pae-readback-20261006-v1')


def write_new(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def main():
    for path in (WAIT, LAUNCH, EXECUTION, RECEIPT):
        assert not path.exists(), f'Immutable launch artifact already exists: {path}'
    assert not OUTPUT.exists(), f'Fresh output required: {OUTPUT}'
    assert not PRODUCER_RECEIPT.exists(), 'Producer already completed; do not queue a delayed handoff'
    producer = json.loads(PRODUCER_LAUNCH.read_text())
    assert live(producer), 'Require the exact original downloader to be live when queueing the handoff'
    resource = json.loads(RESOURCES.read_text())
    assert (resource['cpus'], resource['memory_gib'], resource['swap_gib'], resource['gpu']) == (2, 16, 0, False)
    reader = [sys.executable, 'scripts/readback_full_atlas_missing_pae_v1.py',
              '--producer-receipt', str(PRODUCER_RECEIPT), '--producer-plan', str(PRODUCER_PLAN),
              '--output', str(OUTPUT), '--receipt', str(RECEIPT)]
    controller = [sys.executable, 'scripts/run_weighted_fit_controller_software_stage.py',
                  '--resources', str(RESOURCES), '--execution', str(EXECUTION), '--receipt', str(RECEIPT),
                  '--', *reader]
    pins = {str(path): sha(path) for path in (
        PRODUCER_LAUNCH, PRODUCER_PLAN, RESOURCES,
        Path('scripts/run_after_verified_dependencies_v2.py'),
        Path('scripts/run_weighted_fit_controller_software_stage.py'),
        Path('scripts/readback_full_atlas_missing_pae_v1.py'), Path(__file__))}
    write_new(WAIT, dict(dependencies=[str(PRODUCER_LAUNCH)], command=controller, pins=pins,
                         scope='Wait for the exact original complete missing-AFDB PAE downloader, then execute one bounded independent whole-queue readback. No download, cache mutation, prediction, structure search, or biological inference.'))
    command = ['systemd-run', '--user', '--collect', '--unit=' + UNIT,
               '--working-directory=' + str(Path.cwd()), '--property=CPUQuota=200%',
               '--property=MemoryMax=16G', '--property=MemorySwapMax=0',
               '--setenv=PYTHONUNBUFFERED=1', '--setenv=OPENBLAS_NUM_THREADS=1',
               '--setenv=OMP_NUM_THREADS=1', '--setenv=MKL_NUM_THREADS=1',
               sys.executable, 'scripts/run_after_verified_dependencies_v2.py', '--plan', str(WAIT)]
    subprocess.run(command, check=True)
    pid = 0
    for _ in range(20):
        pid = int(subprocess.check_output(['systemctl', '--user', 'show', UNIT, '-p', 'MainPID', '--value'], text=True))
        if pid:
            break
        time.sleep(.1)
    assert pid > 0
    process = psutil.Process(pid)
    expected = [sys.executable, 'scripts/run_after_verified_dependencies_v2.py', '--plan', str(WAIT)]
    assert process.cmdline() == expected
    cgroup = subprocess.check_output(['systemctl', '--user', 'show', UNIT, '-p', 'ControlGroup', '--value'], text=True).strip()
    limits = {key: (Path('/sys/fs/cgroup') / cgroup.lstrip('/') / key).read_text().strip()
              for key in ('cpu.max', 'memory.max', 'memory.swap.max')}
    assert limits == {'cpu.max': '200000 100000', 'memory.max': str(16 * 2**30), 'memory.swap.max': '0'}
    write_new(LAUNCH, dict(unit=UNIT, pid=pid, created=process.create_time(), cmdline=process.cmdline(),
                            plan=str(WAIT), plan_sha256=sha(WAIT), checked_utc=datetime.now(timezone.utc).isoformat(),
                            actual_cgroup_limits=limits, systemd_launch_command=command, scientific_eligibility=False,
                            scope='Bounded dependency waiter for independent full missing-PAE retrieval integrity readback; not a biological analysis.'))
    print(json.dumps(dict(unit=UNIT, pid=pid, launch=str(LAUNCH)), indent=2))


if __name__ == '__main__':
    main()
