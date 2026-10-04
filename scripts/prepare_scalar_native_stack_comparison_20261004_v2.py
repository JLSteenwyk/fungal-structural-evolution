#!/usr/bin/env python3
"""Change only the separate native diagnostic's process-local stack soft limit."""
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify


def write(path, value):
    with Path(path).open('x') as handle:json.dump(value, handle, indent=2);handle.write('\n')


def main():
    original = Path('metadata/scalar_native_signal_debugger_plan_20261004_v1.json')
    plan = json.loads(original.read_text())
    pins = dict(plan['pins'])
    verify(pins)
    software_path = Path('metadata/scalar_native_signal_debugger_software_transport_20261004_v2.json')
    software = json.loads(software_path.read_text())
    assert software['original_tool_session_id'] == 92249 and software['original_tool_terminal_exit_code'] == 0
    assert software['whole_wrapper_initial_and_terminal_payloads_matched']
    verify(software['source_hashes'])
    fault_path = Path('metadata/scalar_native_stack_fault_20261004_v1.json')
    fault = json.loads(fault_path.read_text())
    verify(fault['source_hashes'])
    assert fault['stack_exhaustion_supported_in_this_diagnostic'] and not fault['validated_fix']
    assert fault['stack_soft_limit_bytes'] == 8*2**20
    for path in [original, software_path, fault_path]:bind(pins, path)
    project_sources(pins, [Path(__file__), Path('scripts/run_scalar_native_signal_debugger_v2.py')])
    caps = dict(plan['native_caps'], stack_soft_bytes=64*2**20)
    resources = dict(plan['resources'], native_caps=caps,
        scope='One separate controlled stack comparison: unchanged original native command, '
              'input/prior/seed/20iterations and48GiB AS/5674CPU-second/2GiB-file limits; '
              'only native stack soft limit increases8to64MiB, hard limit remains unlimited. '
              'TwoCPU/64GiB/no swap/BLAS1 debugger wrapper. No global limits, sources, '
              'priors, production attempts or failed outputs changed. Runtime uncalibrated; '
              'a normal native exit still requires full independent output validation.')
    rp = Path('metadata/scalar_native_stack_comparison_resources_20261004_v2.json')
    write(rp, resources)
    bind(pins, rp)
    plan.update(status='controlled_process_local_stack_comparison_prepared',
        prepared_utc=datetime.now(timezone.utc).isoformat(),
        output='results/ancestral/scalar-native-stack-comparison-20261004-v2',
        native_caps=caps, resources=resources, pins=pins, scope=resources['scope'],
        preceding_captured_fault=str(fault_path), production_stack_limits_changed=False)
    verify(pins)
    pp = Path('metadata/scalar_native_stack_comparison_plan_20261004_v2.json')
    write(pp, plan)
    print(json.dumps(dict(plan=str(pp), plan_sha256=sha(pp), native_caps=caps,
        original_native_command_unchanged=True, production_stack_limits_changed=False), indent=2))


if __name__ == '__main__':main()
