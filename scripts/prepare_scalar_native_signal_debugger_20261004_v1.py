#!/usr/bin/env python3
"""Freeze a separate OG0000972 native-stop diagnostic after actual controls."""
from datetime import datetime, timezone
import json
from pathlib import Path
import statistics

import psutil

from ancestral_chain_attempt import sha
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify


def write(path, value):
    with Path(path).open('x') as handle:json.dump(value, handle, indent=2);handle.write('\n')


def main():
    tp = Path('metadata/scalar_native_signal_debugger_software_transport_20261004_v1.json')
    transport = json.loads(tp.read_text())
    assert transport['original_tool_session_id'] == 77247 and transport['original_tool_terminal_exit_code'] == 0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    verify(transport['source_hashes'])
    failures_path = Path('metadata/baliphy_scalar_v6_failed_attempts_recheck_20261004_v1.json')
    failures = json.loads(failures_path.read_text())
    verify(failures['source_hashes'])
    assert failures['failed_attempts'] == 12 and failures['exit_code_counts'] == {'-11': 12}
    row = next(r for r in failures['attempts'] if r['prior'] == 'broad' and r['role'] == 1)
    command = row['original_identity']['command']
    assert command[:2] == ['/usr/bin/prlimit', '--as='+str(48*2**30)]
    native = command[command.index('--')+1:]
    sp = Path('metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json')
    sampler = json.loads(sp.read_text())
    pins = dict(sampler['pins'])
    verify(pins)
    for path in [tp, failures_path, sp, Path(native[0]), Path('/usr/bin/gdb'), Path('/usr/bin/prlimit')]:
        bind(pins, path)
    project_sources(pins, [Path(__file__), Path('scripts/run_scalar_native_signal_debugger_v1.py')])
    caps = dict(address_space_bytes=48*2**30, cpu_seconds=5674, per_file_bytes=2*2**30)
    assert '--cpu=5674' in command and '--fsize='+str(2*2**30) in command
    resources = dict(cpus=2, memory_gib=64, swap_gib=0, blas_threads=1,
        address_space_gib=64, cpu_seconds_per_stage=7200, wall_seconds_per_stage=7500,
        per_file_limit_mib=2048, minimum_available_ram_gib=64, minimum_free_disk_gib=256,
        observed_original_role_wall_seconds=row['elapsed_worker_seconds'],
        median_original_failed_role_wall_seconds=statistics.median(r['elapsed_worker_seconds'] for r in failures['attempts']),
        prospective_wall_minutes=[2,15], timing_uncalibrated_under_debugger=True,
        native_caps=caps, output_allowance_gib=4, new_cost_usd=0, gpu=False,
        scope='One separate debugger inferior on unchanged failed input/prior/seed, with original '
              '48GiB AS/5674CPU-second/2GiB-file native limits. Debugger/wrapper share64GiB '
              'cgroup and2CPU/no swap. Original failures took about two minutes; debugger '
              'execution time is uncalibrated. No original retry, production change, paid service or posterior acceptance.')
    assert psutil.virtual_memory().available >= 64*2**30 and psutil.disk_usage('.').free >= 256*2**30
    rp = Path('metadata/scalar_native_signal_debugger_resources_20261004_v1.json')
    write(rp, resources)
    bind(pins, rp)
    plan = dict(status='separate_native_signal_diagnostic_prepared_after_actual_controls',
        prepared_utc=datetime.now(timezone.utc).isoformat(),
        output='results/ancestral/scalar-native-signal-debugger-20261004-v1',
        native_command=native, native_caps=caps, debugger_wall_seconds=7200,
        original_failed_receipt=row['native_receipt'], original_chain_id=row['chain_id'],
        resources=resources, artificial_control=False, pins=pins,
        scope=resources['scope']+' Signal stop captures PC, siginfo, backtrace, registers, '
            'instructions, shared libraries, process limits/status/maps. Debugger changes '
            'execution context; failure or nonrecurrence is not itself a root-cause repair. '
            'No core dump is requested. All remaining full-grid roles continue independently.')
    verify(pins)
    pp = Path('metadata/scalar_native_signal_debugger_plan_20261004_v1.json')
    write(pp, plan)
    print(json.dumps(dict(plan=str(pp), plan_sha256=sha(pp), original_chain=row['chain_id'], resources=resources), indent=2))


if __name__ == '__main__':main()
