#!/usr/bin/env python3
"""Launch isolated memory recovery with one CPU and a verified 64 GiB cap."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import psutil
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); a = p.parse_args()
    plan = json.loads(a.plan.read_text()); lp = Path('metadata/baliphy_memory_recovery_launch_20261002.json')
    assert not lp.exists(), 'Existing original recovery launch must remain immutable'
    for path, digest in plan['pins'].items(): assert sha(path) == digest, path
    assert psutil.virtual_memory().available >= 64 * 2**30
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    command = [sys.executable, 'scripts/recover_baliphy_memory_failures_20261002.py', '--plan', str(a.plan)]
    launch = ['systemd-run', '--user', '--collect', '--unit=' + plan['unit'], '--working-directory=' + str(Path.cwd()),
        '--property=CPUQuota=100%', '--property=MemoryMax=64G', '--property=MemorySwapMax=0',
        '--setenv=PYTHONUNBUFFERED=1', '--setenv=OPENBLAS_NUM_THREADS=1', '--setenv=OMP_NUM_THREADS=1', '--setenv=MKL_NUM_THREADS=1', *command]
    subprocess.run(launch, check=True)
    pid = 0
    for _ in range(20):
        pid = int(subprocess.check_output(['systemctl', '--user', 'show', plan['unit'], '-p', 'MainPID', '--value'], text=True))
        if pid: break
        time.sleep(.1)
    assert pid > 0; process = psutil.Process(pid); assert process.cmdline() == command
    cgroup = subprocess.check_output(['systemctl', '--user', 'show', plan['unit'], '-p', 'ControlGroup', '--value'], text=True).strip()
    limits = {k: (Path('/sys/fs/cgroup') / cgroup.lstrip('/') / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert limits == {'cpu.max': '100000 100000', 'memory.max': str(64 * 2**30), 'memory.swap.max': '0'}
    record = dict(unit=plan['unit'], pid=pid, created=process.create_time(), cmdline=command,
        plan=str(a.plan), plan_sha256=sha(a.plan), checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_cgroup_limits=limits, systemd_launch_command=launch)
    with lp.open('x') as handle: handle.write(json.dumps(record, indent=2) + '\n')
    print(json.dumps(dict(launch=str(lp), pid=pid, limits=limits)), flush=True)


if __name__ == '__main__': main()
