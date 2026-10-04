#!/usr/bin/env python3
"""Prepare one observer-instrumented full-input diagnostic after paired native controls."""
from datetime import datetime, timezone
import json
from pathlib import Path

import psutil

from ancestral_chain_attempt import sha
from baliphy_logger_phase_trace_v1 import reverse, transform
from reference_measurement_union_sources import bind, verify


def save(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def main():
    plan_path = Path('metadata/baliphy_logger_phase_diagnostic_plan_20261004_v1.json')
    resources_path = Path('metadata/baliphy_logger_phase_diagnostic_resources_20261004_v1.json')
    assert not plan_path.exists() and not resources_path.exists()
    gate_path = Path('metadata/baliphy_logger_phase_trace_software_validation_20261004_v1.json')
    transport_path = Path('metadata/baliphy_logger_phase_trace_software_transport_20261004_v1.json')
    gate, transport = [json.loads(p.read_text()) for p in [gate_path, transport_path]]
    assert gate['status'] == 'passed_reversible_logger_phase_trace_and_paired_native_controls_v1'
    assert gate['full_programs_reversibly_checked'] == 405 and gate['native_prior_controls'] == 3
    assert transport['validation_sha256'] == sha(gate_path)
    assert transport['original_tool_terminal_exit_code'] == 0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    pins = dict(transport['source_hashes'])
    verify(pins)
    driver_gate_path = Path('metadata/scalar_native_signal_debugger_software_transport_20261004_v2.json')
    driver_gate = json.loads(driver_gate_path.read_text())
    assert driver_gate['original_tool_terminal_exit_code'] == 0
    verify(driver_gate['source_hashes'])
    pins.update(driver_gate['source_hashes'])
    baseline_path = Path('metadata/scalar_native_signal_debugger_plan_20261004_v1.json')
    baseline = json.loads(baseline_path.read_text())
    verify(baseline['pins'])
    pins.update(baseline['pins'])
    original = Path(baseline['native_command'][baseline['native_command'].index('run') + 1])
    prepared = Path('results/ancestral/prepared-logger-phase-diagnostic-20261004-v1')
    prepared.mkdir(exist_ok=False)
    program = (prepared / original.name).resolve()
    source = original.read_text()
    program.write_text(transform(source))
    assert reverse(program.read_text()) == source
    command = list(baseline['native_command'])
    command[command.index('run') + 1] = str(program)
    assert [(i, a, b) for i, (a, b) in enumerate(zip(command, baseline['native_command'])) if a != b] == [
        (command.index('run') + 1, str(program), str(original))]
    assert command[command.index('--seed') + 1] == '1322110685'
    assert command[command.index('--iterations') + 1] == '20'
    original_attempt = Path(baseline['original_failed_receipt'])
    assert json.loads(original_attempt.read_text())['exit_code'] == -11
    native_caps = dict(baseline['native_caps'], stack_soft_bytes=8 * 2**20)
    assert native_caps['address_space_bytes'] == 48 * 2**30
    assert native_caps['cpu_seconds'] == 5674 and native_caps['per_file_bytes'] == 2 * 2**30
    scope = ('One separate full-input failed-family diagnostic, with reversible IO phase markers only. '
             'Original probability expressions, alignment, tree, prior, seed1322110685,20iterations, '
             '48GiB AS/5674CPU-seconds/2GiB file/8MiB soft-unlimited hard stack remain. '
             'No original attempt is restarted and no full-sampling scope is reduced. '
             'Marker noninterference is qualified on three existing tiny controls, not assumed on full input. '
             'Debugger exit0 and null callback status do not prove a completed native sampler; '
             'original signal text, output horizon, and independent readback must be reviewed. '
             'No posterior or biological acceptance, GPU, paid resource, or installed-library change.')
    resources = dict(prepared_utc=datetime.now(timezone.utc).isoformat(),cpus=2,memory_gib=64,
                     swap_gib=0,blas_threads=1,address_space_gib=64,cpu_seconds_per_stage=7200,
                     wall_seconds_per_stage=7500,per_file_limit_mib=2048,
                     minimum_available_ram_gib=64,minimum_free_disk_gib=256,
                     available_ram_gib=psutil.virtual_memory().available/2**30,
                     free_disk_gib=psutil.disk_usage('.').free/2**30,
                     output_allowance_gib=4,prospective_wall_minutes=[2,15],
                     prospective_debugger_runtime_uncalibrated=True,
                     prior_original_failed_role_wall_seconds=baseline['resources']['observed_original_role_wall_seconds'],
                     prior_64_mib_stack_diagnostic_minutes=94.7,
                     maximum_native_cpu_seconds=5674,new_cost_usd=0,gpu=False,scope=scope)
    assert resources['available_ram_gib'] >= 64 and resources['free_disk_gib'] >= 256
    save(resources_path, resources)
    for path in [Path(__file__), Path('scripts/baliphy_logger_phase_trace_v1.py'),
                 Path('scripts/run_scalar_native_signal_debugger_v2.py'), program, original,
                 gate_path, transport_path, driver_gate_path, baseline_path,
                 original_attempt, resources_path, Path('/usr/bin/gdb'), Path('/usr/bin/prlimit')]:
        bind(pins, path)
    verify(pins)
    plan = dict(status='prepared_separate_qualified_logger_phase_diagnostic',
                prepared_utc=datetime.now(timezone.utc).isoformat(),
                output='results/ancestral/baliphy-logger-phase-diagnostic-20261004-v1',
                native_command=command,native_caps=native_caps,debugger_wall_seconds=7200,
                original_model=str(original),original_model_sha256=sha(original),
                diagnostic_model=str(program),diagnostic_model_sha256=sha(program),
                exact_model_reversibility_checked=True,original_failed_receipt=str(original_attempt),
                artificial_control=False,resources=resources,pins=pins,scope=scope)
    assert not Path(plan['output']).exists()
    save(plan_path, plan)
    print(json.dumps({k:v for k,v in plan.items() if k not in ['pins','resources']}, indent=2))


if __name__ == '__main__':
    main()
