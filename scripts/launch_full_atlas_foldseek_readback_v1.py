#!/usr/bin/env python3
"""Queue the full-atlas native-database readback after the exact producer exits.

The waiting unit is deliberately small: it uses the same declared ceiling as
the readback stage, consumes no GPU, and verifies the producer's original
systemd invocation before it starts the bounded reader.
"""
import json
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
import time

import psutil

from run_after_verified_dependencies import live
from run_ortholog_pair_guide_comparison import sha


ROOT = Path.cwd()
UNIT = 'fungal-full-atlas-foldseek-readback-20261006-v1.service'
WAIT_PLAN = Path('metadata/full_atlas_foldseek_readback_wait_plan_20261006_v1.json')
LAUNCH = Path('metadata/full_atlas_foldseek_readback_launch_20261006_v1.json')
PRODUCER_LAUNCH = Path('metadata/full_atlas_foldseek_launch_20261006_v2.json')
PRODUCER_PLAN = Path('metadata/full_atlas_foldseek_plan_20261006_v2.json')
PRODUCER_RECEIPT = Path('metadata/full_atlas_foldseek_database_20261006_v2.json')
RESOURCES = Path('metadata/full_atlas_foldseek_readback_resources_20261006_v1.json')
EXECUTION = Path('metadata/full_atlas_foldseek_readback_execution_20261006_v1.json')
RECEIPT = Path('metadata/full_atlas_foldseek_readback_20261006_v1.json')
OUTPUT = Path('results/structural_clusters/full-prediction-atlas-database-readback-20261006-v1')


def write_new(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def main():
    for path in (WAIT_PLAN, LAUNCH, EXECUTION, RECEIPT):
        assert not path.exists(), f'Immutable launch artifact already exists: {path}'
    assert not OUTPUT.exists(), f'Fresh output required: {OUTPUT}'
    assert not PRODUCER_RECEIPT.exists(), 'The producer already completed; do not create a delayed handoff'
    producer = json.loads(PRODUCER_LAUNCH.read_text())
    assert live(producer), 'Require the exact original producer to be live when queueing the handoff'
    resources = json.loads(RESOURCES.read_text())
    assert (resources['cpus'], resources['memory_gib'], resources['swap_gib'], resources['gpu']) == (2, 16, 0, False)
    reader = [sys.executable, 'scripts/readback_full_atlas_foldseek_database_v1.py',
              '--producer-receipt', str(PRODUCER_RECEIPT), '--producer-plan', str(PRODUCER_PLAN),
              '--output', str(OUTPUT), '--receipt', str(RECEIPT)]
    controller = [sys.executable, 'scripts/run_weighted_fit_controller_software_stage.py',
                  '--resources', str(RESOURCES), '--execution', str(EXECUTION), '--receipt', str(RECEIPT),
                  '--', *reader]
    pins = {str(path): sha(path) for path in (
        PRODUCER_LAUNCH, PRODUCER_PLAN, RESOURCES,
        Path('scripts/run_after_verified_dependencies_v2.py'),
        Path('scripts/run_weighted_fit_controller_software_stage.py'),
        Path('scripts/readback_full_atlas_foldseek_database_v1.py'),
        Path(__file__),
    )}
    write_new(WAIT_PLAN, dict(
        dependencies=[str(PRODUCER_LAUNCH)], command=controller, pins=pins,
        scope='Wait for the original live full mixed-source Foldseek database producer to complete successfully, then run one independent, bounded, whole-database native-record readback. No conversion, prediction, search, clustering, or biological inference.'
    ))
    command = ['systemd-run', '--user', '--collect', '--unit=' + UNIT,
               '--working-directory=' + str(ROOT), '--property=CPUQuota=200%',
               '--property=MemoryMax=16G', '--property=MemorySwapMax=0',
               '--setenv=PYTHONUNBUFFERED=1', '--setenv=OPENBLAS_NUM_THREADS=1',
               '--setenv=OMP_NUM_THREADS=1', '--setenv=MKL_NUM_THREADS=1',
               sys.executable, 'scripts/run_after_verified_dependencies_v2.py', '--plan', str(WAIT_PLAN)]
    subprocess.run(command, check=True)
    pid = 0
    for _ in range(20):
        pid = int(subprocess.check_output(['systemctl', '--user', 'show', UNIT, '-p', 'MainPID', '--value'], text=True))
        if pid:
            break
        time.sleep(.1)
    assert pid > 0, 'Waiting unit did not expose a main process'
    process = psutil.Process(pid)
    expected = [sys.executable, 'scripts/run_after_verified_dependencies_v2.py', '--plan', str(WAIT_PLAN)]
    assert process.cmdline() == expected
    cgroup = subprocess.check_output(['systemctl', '--user', 'show', UNIT, '-p', 'ControlGroup', '--value'], text=True).strip()
    limits = {key: (Path('/sys/fs/cgroup') / cgroup.lstrip('/') / key).read_text().strip()
              for key in ('cpu.max', 'memory.max', 'memory.swap.max')}
    assert limits == {'cpu.max': '200000 100000', 'memory.max': str(16 * 2**30), 'memory.swap.max': '0'}
    write_new(LAUNCH, dict(unit=UNIT, pid=pid, created=process.create_time(), cmdline=process.cmdline(),
                            plan=str(WAIT_PLAN), plan_sha256=sha(WAIT_PLAN),
                            checked_utc=datetime.now(timezone.utc).isoformat(), actual_cgroup_limits=limits,
                            systemd_launch_command=command, scientific_eligibility=False,
                            scope='Bounded dependency waiter for an independent database integrity readback; not a biological analysis.'))
    print(json.dumps({'unit': UNIT, 'pid': pid, 'launch': str(LAUNCH)}, indent=2))


if __name__ == '__main__':
    main()
