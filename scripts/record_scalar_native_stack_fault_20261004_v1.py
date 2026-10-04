#!/usr/bin/env python3
"""Verify the captured instruction, stack mapping, fault address and original caps."""
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    tp = Path('metadata/scalar_native_signal_debugger_transport_20261004_v1.json')
    transport = json.loads(tp.read_text())
    assert transport['original_tool_session_id'] == 19850 and transport['original_tool_terminal_exit_code'] == 0
    verify(transport['source_hashes'])
    rp = Path('metadata/scalar_native_signal_debugger_20261004_v1.json')
    receipt = json.loads(rp.read_text())
    verify(receipt['source_hashes'])
    verify(receipt['artifacts'])
    assert receipt['native_stop']['signal'] == 'SIGSEGV' and receipt['native_stop']['live_stop']
    root = Path('results/ancestral/scalar-native-signal-debugger-20261004-v1')
    log = (root/'gdb_stdout.log').read_text()
    fault = int(re.search(r'_sigfault = \{si_addr = (0x[0-9a-f]+)', log).group(1), 16)
    rsp = int(re.search(r'^rsp\s+(0x[0-9a-f]+)', log, re.M).group(1), 16)
    current = next(line for line in log.splitlines() if line.startswith('=> '))
    assert '\tcall ' in current and 'si_code = 1' in log
    mapping = next(line for line in (root/'native-stop-maps.txt').read_text().splitlines() if '[stack]' in line)
    lo, hi = [int(v, 16) for v in mapping.split()[0].split('-')]
    limits = (root/'native-stop-limits.txt').read_text()
    stack = next(line for line in limits.splitlines() if line.startswith('Max stack size'))
    assert stack.split()[3:5] == [str(8*2**20), 'unlimited']
    assert hi-lo == 8*2**20 and fault == rsp-8 and fault < lo
    status = (root/'native-stop-status.txt').read_text()
    values = {key: int(re.search(r'^'+key+r':\s+(\d+) kB', status, re.M).group(1))*1024
              for key in ['VmPeak', 'VmSize', 'VmHWM', 'VmRSS', 'VmStk']}
    assert values['VmStk'] == 8*2**20 and values['VmPeak'] < receipt['native_caps']['address_space_bytes']
    pins = {}
    for path in [Path(__file__), tp, rp, *map(Path, receipt['artifacts'])]:bind(pins, path)
    verify(pins)
    result = dict(status='verified_captured_native_stack_boundary_fault',
        checked_utc=datetime.now(timezone.utc).isoformat(), original_tool_session_id=19850,
        native_signal='SIGSEGV', program_counter=receipt['native_stop']['pc'],
        current_instruction=current, fault_address=hex(fault), stack_pointer=hex(rsp),
        fault_equals_call_return_address_write=True, stack_mapping=[hex(lo), hex(hi)],
        fault_bytes_below_stack_mapping=lo-fault, stack_soft_limit_bytes=8*2**20,
        stack_mapping_bytes=hi-lo, captured_memory_bytes=values,
        original_address_space_cap_bytes=receipt['native_caps']['address_space_bytes'],
        stack_exhaustion_supported_in_this_diagnostic=True,
        all_original_twelve_causes_proven=False, validated_fix=False,
        source_hashes=pins, original_attempts_restarted=False, posterior_qualified=False,
        scientific_eligibility=False, gpu=False,
        scope='Captured x86-64 call writes at rsp-8 below the fully grown8MiB stack '
              'mapping with8MiB soft limit and SIGSEGV MAPERR. Repeated evaluator '
              'frames support stack exhaustion in this separate diagnostic. This '
              'does not prove all original twelve crashes share that cause or '
              'that a larger process-local stack resolves full sampling. No global '
              'limit/source/prior change or posterior acceptance.')
    with Path('metadata/scalar_native_stack_fault_20261004_v1.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':main()
