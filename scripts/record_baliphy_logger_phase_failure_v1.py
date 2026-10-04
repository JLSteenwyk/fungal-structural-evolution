#!/usr/bin/env python3
"""Attribute the captured diagnostic stop to its last completed logger phase only."""
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from ancestral_chain_attempt import sha
from baliphy_logger_phase_trace_v1 import events
from reference_measurement_union_sources import bind, verify


def main():
    output=Path('metadata/baliphy_logger_phase_failure_review_20261004_v1.json')
    assert not output.exists()
    transport_path=Path('metadata/baliphy_logger_phase_diagnostic_transport_20261004_v1.json')
    transport=json.loads(transport_path.read_text())
    pins=dict(transport['source_hashes'])
    verify(pins)
    assert transport['original_tool_terminal_exit_code']==0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    receipt_path=Path('metadata/baliphy_logger_phase_diagnostic_20261004_v1.json')
    receipt=json.loads(receipt_path.read_text())
    assert receipt['native_stop']['live_stop'] and receipt['native_stop']['signal']=='SIGSEGV'
    root=Path('results/ancestral/baliphy-logger-phase-diagnostic-20261004-v1')
    trace=events((root/'gdb_stderr.log').read_text())
    assert trace==[
        dict(logger='joint-states',phase='context',boundary='begin',iteration=0),
        dict(logger='joint-states',phase='context',boundary='end',iteration=0),
        dict(logger='joint-states',phase='write',boundary='begin',iteration=0)]
    text=(root/'gdb_stdout.log').read_text()
    fault=int(re.search(r'_sigfault = \{si_addr = (0x[0-9a-f]+)',text).group(1),16)
    rsp=int(re.search(r'^rsp\s+(0x[0-9a-f]+)',text,re.M).group(1),16)
    mapping=next(line for line in (root/'native-stop-maps.txt').read_text().splitlines() if line.endswith('[stack]'))
    low,high=[int(value,16) for value in mapping.split()[0].split('-')]
    assert fault==rsp-8 and fault<low and high-low==8*2**20
    limit=next(line for line in (root/'native-stop-limits.txt').read_text().splitlines() if line.startswith('Max stack size'))
    assert limit.split()[3:5]==[str(8*2**20),'unlimited']
    scalar=root/'independent-chain-1/C1.log'
    joint=root/'independent-chain-1/C1.P1.site-property-samples.jsonl'
    assert scalar.stat().st_size==joint.stat().st_size==0
    for path in [Path(__file__),transport_path,receipt_path]:bind(pins,path)
    verify(pins)
    result=dict(status='verified_instrumented_full_input_stack_fault_inside_joint_write',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_tool_session_id=64683,
        debugger_wrapper_exit_code=0,native_signal='SIGSEGV',native_exit_zero_proven=False,
        last_phase_events=trace,fault_address_hex=hex(fault),stack_pointer_hex=hex(rsp),
        fault_bytes_below_stack_mapping=low-fault,stack_mapping_bytes=high-low,
        actual_stack_soft_bytes=8*2**20,actual_stack_hard_unlimited=True,
        scalar_rows_saved=0,joint_frames_saved=0,planned_iterations=20,
        planned_horizon_completed=False,original_uninstrumented_failure_iteration_equivalent=False,
        exact_original_model_reversibility_verified=True,
        root_cause_of_all_original_failures_proven=False,validated_repair=False,
        source_hashes=pins,original_jobs_restarted=False,scientific_eligibility=False,
        posterior_qualified=False,gpu=False,
        scope='The separate instrumented full-input run exhausts its original8MiB stack '
              'during iteration0 joint logger write, before any joint frame or scalar row. '
              'Original uninstrumented diagnostics reached later iterations; IO instrumentation '
              'and debugger context can change evaluation order. This localizes this diagnostic '
              'only, without attributing a particular joint field or all24original failures.')
    with output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
