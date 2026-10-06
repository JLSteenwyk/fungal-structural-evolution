#!/usr/bin/env python3
"""Queue the independent full coordinate-profile readback after the atlas build."""
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
UNIT = "fungal-full-atlas-coordinate-profiles-readback-20261006-v1.service"
RESOURCES = Path("metadata/full_atlas_coordinate_profiles_readback_resources_20261006_v1.json")
PLAN = Path("metadata/full_atlas_coordinate_profiles_readback_plan_20261006_v1.json")
WAIT = Path("metadata/full_atlas_coordinate_profiles_readback_wait_plan_20261006_v1.json")
LAUNCH = Path("metadata/full_atlas_coordinate_profiles_readback_launch_20261006_v1.json")
EXECUTION = Path("metadata/full_atlas_coordinate_profiles_readback_execution_20261006_v1.json")
RECEIPT = Path("metadata/full_atlas_coordinate_profiles_readback_20261006_v1.json")
OUTPUT = Path("results/structures/full-atlas-coordinate-profiles-readback-20261006-v1")
PRODUCER_PLAN = Path("metadata/full_atlas_coordinate_profiles_plan_20261005_v1.json")
PRODUCER_RECEIPT = Path("metadata/full_atlas_coordinate_profiles_20261005_v1.json")
PRODUCER_TRANSPORT = Path("metadata/full_atlas_coordinate_profiles_transport_20261006_v1.json")
ATLAS_LAUNCH = Path("metadata/full_atlas_foldseek_launch_20261006_v2.json")


def write_new(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")


def main():
    for path in (PLAN, WAIT, LAUNCH, EXECUTION, RECEIPT):
        if path.exists():
            raise FileExistsError(path)
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    atlas = json.loads(ATLAS_LAUNCH.read_text())
    if not live(atlas):
        raise RuntimeError("Queue only while the exact Foldseek producer is live")
    resources = json.loads(RESOURCES.read_text())
    if (resources["cpus"], resources["memory_gib"], resources["swap_gib"], resources["gpu"]) != (8, 64, 0, False):
        raise ValueError("Unexpected independent-reader resource envelope")
    sources = (PRODUCER_PLAN, PRODUCER_RECEIPT, PRODUCER_TRANSPORT, RESOURCES,
               Path("scripts/readback_full_atlas_coordinate_profiles_v1.py"), Path(__file__))
    pins = {str(path): sha(path) for path in sources}
    write_new(PLAN, {
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "output": str(OUTPUT),
        "producer_plan": str(PRODUCER_PLAN),
        "producer_receipt": str(PRODUCER_RECEIPT),
        "producer_transport": str(PRODUCER_TRANSPORT),
        "resources": {"cpu": resources["cpus"]},
        "pins": pins,
        "scope": "Independent complete coordinate-profile readback. The job starts only after the active full Foldseek producer exits, avoiding two simultaneous complete scans of the source CIF atlas.",
    })
    reader = [sys.executable, "scripts/readback_full_atlas_coordinate_profiles_v1.py", "--plan", str(PLAN), "--receipt", str(RECEIPT)]
    controller = [sys.executable, "scripts/run_weighted_fit_controller_software_stage.py", "--resources", str(RESOURCES), "--execution", str(EXECUTION), "--receipt", str(RECEIPT), "--", *reader]
    wait_pins = {str(path): sha(path) for path in (ATLAS_LAUNCH, PLAN, RESOURCES,
                 Path("scripts/run_after_verified_dependencies_v2.py"),
                 Path("scripts/run_weighted_fit_controller_software_stage.py"),
                 Path("scripts/readback_full_atlas_coordinate_profiles_v1.py"), Path(__file__))}
    write_new(WAIT, {"dependencies": [str(ATLAS_LAUNCH)], "command": controller, "pins": wait_pins,
                     "scope": "Wait for the exact live full-atlas Foldseek producer, then execute one bounded, independent complete coordinate-profile readback. No prediction, PAE retrieval, source/archive mutation, clustering or biological inference."})
    command = ["systemd-run", "--user", "--collect", "--unit=" + UNIT, "--working-directory=" + str(ROOT),
               "--property=CPUQuota=800%", "--property=MemoryMax=64G", "--property=MemorySwapMax=0",
               "--setenv=PYTHONUNBUFFERED=1", "--setenv=OPENBLAS_NUM_THREADS=1", "--setenv=OMP_NUM_THREADS=1",
               "--setenv=MKL_NUM_THREADS=1", sys.executable, "scripts/run_after_verified_dependencies_v2.py", "--plan", str(WAIT)]
    subprocess.run(command, check=True)
    pid = 0
    for _ in range(20):
        pid = int(subprocess.check_output(["systemctl", "--user", "show", UNIT, "-p", "MainPID", "--value"], text=True))
        if pid:
            break
        time.sleep(0.1)
    if not pid:
        raise RuntimeError("Waiting unit did not expose a main process")
    process = psutil.Process(pid)
    expected = [sys.executable, "scripts/run_after_verified_dependencies_v2.py", "--plan", str(WAIT)]
    if process.cmdline() != expected:
        raise RuntimeError("Waiting process differs")
    cgroup = subprocess.check_output(["systemctl", "--user", "show", UNIT, "-p", "ControlGroup", "--value"], text=True).strip()
    limits = {key: (Path("/sys/fs/cgroup") / cgroup.lstrip("/") / key).read_text().strip()
              for key in ("cpu.max", "memory.max", "memory.swap.max")}
    expected_limits = {"cpu.max": "800000 100000", "memory.max": str(64 * 2**30), "memory.swap.max": "0"}
    if limits != expected_limits:
        raise RuntimeError("Waiting cgroup limits differ")
    write_new(LAUNCH, {"unit": UNIT, "pid": pid, "created": process.create_time(), "cmdline": process.cmdline(),
                       "plan": str(WAIT), "plan_sha256": sha(WAIT), "checked_utc": datetime.now(timezone.utc).isoformat(),
                       "actual_cgroup_limits": limits, "systemd_launch_command": command, "scientific_eligibility": False,
                       "scope": "Bounded dependency waiter for independent full coordinate-profile archive and source-atom readback; not a biological analysis."})
    print(json.dumps({"unit": UNIT, "pid": pid, "launch": str(LAUNCH)}, indent=2))


if __name__ == "__main__":
    main()
