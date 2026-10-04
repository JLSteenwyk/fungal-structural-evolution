#!/usr/bin/env python3
"""Prepare one read-only independent reader gated on the original comparison producer."""
import json
from pathlib import Path
import sys

import psutil

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import verify


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n')


def main():
    prefix = 'metadata/baliphy_joint_fasta_v7_failure_grid_'
    pp = Path(prefix + 'plan_20261004_v1.json')
    config_path = Path(prefix + 'execution_20261004_v1/configuration.json')
    initial_path = Path(prefix + 'original_tool_initial_20261004_v1.json')
    plan, config, initial = [json.loads(p.read_text()) for p in [pp, config_path, initial_path]]
    verify(plan['pins'])
    assert initial['session_id'] == 72605 and config['invocation_id'] in initial['output']
    wrapper = config['wrapper']; process = psutil.Process(wrapper['pid'])
    assert process.create_time() == wrapper['created'] and process.cmdline() == wrapper['cmdline']
    root = Path(plan['output'])
    assert not (root / 'readback.json').exists()
    launch_path = Path(prefix + 'producer_launch_20261004_v1.json')
    write(launch_path, dict(unit='fungal-joint-fasta-v7-failure-grid-20261004-v1.service',
        **wrapper, plan=str(pp), plan_sha256=sha(pp), invocation_id=config['invocation_id'],
        original_tool_session_id=72605, actual_cgroup_limits=config['actual_cgroup_limits'],
        source_hashes={str(p): sha(p) for p in [pp, config_path, initial_path]}, scientific_eligibility=False))
    resources_path = Path(prefix + 'reader_resources_20261004_v1.json')
    resources = dict(cpus=2, memory_gib=32, swap_gib=0, blas_threads=1,
        address_space_gib=24, cpu_seconds_per_stage=7200, wall_seconds_per_stage=10800,
        per_file_limit_mib=2048, minimum_available_ram_gib=32, minimum_free_disk_gib=256,
        expected_roles=24, native_prediction_workers=0, gpu=False, new_cost_usd=0,
        output_allowance_gib=1, prospective_reader_wall_hours=[0.1, 3],
        runtime_estimate_uncalibrated=True,
        scope='Independent read-only reconstruction of every completed comparison outcome; '
              '2CPU/32GiB/noSwap/24GiB AS, unchanged old readers and no native MCMC or export creation. '
              'Gated on original producer invocation success; failures and review outcomes remain explicit.')
    write(resources_path, resources)
    execution_path = Path(prefix + 'reader_execution_20261004_v1.json')
    wait_path = Path(prefix + 'reader_wait_plan_20261004_v1.json')
    command = [sys.executable, 'scripts/run_weighted_fit_controller_software_stage.py',
               '--resources', str(resources_path), '--execution', str(execution_path),
               '--receipt', str(root / 'readback.json'), '--', sys.executable,
               'scripts/run_baliphy_joint_fasta_v7_failure_grid_v1.py', '--plan', str(pp), '--reader']
    bound = [pp, launch_path, resources_path, Path(plan['jobs']), Path(__file__),
             Path('scripts/run_after_verified_dependencies_v2.py'),
             Path('scripts/run_weighted_fit_controller_software_stage.py'),
             Path('scripts/run_baliphy_joint_fasta_v7_failure_grid_v1.py')]
    pins = {str(p): sha(p) for p in bound}; verify(pins)
    write(wait_path, dict(dependencies=[str(launch_path)], command=command, pins=pins,
        scientific_eligibility=False, no_native_attempts_restarted=True))
    print(json.dumps(dict(wait_plan=str(wait_path), resources=str(resources_path),
                         expected_roles=24, reader_launched=False, scientific_eligibility=False), indent=2))


if __name__ == '__main__':
    main()
