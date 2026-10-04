#!/usr/bin/env python3
"""Separate debugger success from a killed, incomplete native sampler."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from read_baliphy_scalar_json_v6b import compare_tsv
from record_user_error_recheck_20261004_1122 import close_execution
from reference_measurement_union_sources import bind, verify


def main():
    output = Path('metadata/scalar_native_stack_comparison_terminal_review_20261004_v3.json')
    assert not output.exists()
    ep = Path('metadata/scalar_native_stack_comparison_execution_20261004_v2.json')
    pp = Path('metadata/scalar_native_stack_comparison_original_tool_payloads_20261004_v2.json')
    payload = json.loads(pp.read_text())
    execution, transport = close_execution(ep,
        'fungal-scalar-native-stack-comparison-20261004-v2.service',payload)
    rp = Path(execution['receipt']); receipt = json.loads(rp.read_text())
    verify(receipt['source_hashes']); verify(receipt['artifacts'])
    assert receipt['debugger_exit_code']==0
    assert receipt['native_stop']==dict(live_stop=False,native_exit_code=None,pid=0,signal=None)
    root = Path('results/ancestral/scalar-native-stack-comparison-20261004-v2')
    text = (root/'gdb_stdout.log').read_text()
    assert text.count('Program terminated with signal SIGKILL, Killed.')==1
    assert 'The program no longer exists.' in text
    scalar = compare_tsv(root/'independent-chain-1')
    assert scalar['rows']==10 and scalar['nonfinite_reviews']==[]
    lines = (root/'independent-chain-1'/'C1.log').read_text().splitlines()
    iterations = [int(line.split('\t',1)[0]) for line in lines[1:]]
    assert iterations==list(range(10))
    live_path = Path('metadata/scalar_native_stack_comparison_checkpoint_20261004_goal_1300.json')
    live = json.loads(live_path.read_text()); native = live['native'][0]
    assert native['pid']==3215026 and native['native_limits']['cpu']==[5674,5674]
    assert native['native_limits']['stack']==[64*2**20,-1]
    pins = dict(receipt['source_hashes'])
    for p,digest in receipt['artifacts'].items():bind(pins,p,digest)
    for p in [ep,rp,pp,live_path,Path(__file__)]:bind(pins,p)
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl','--user','-u','fungal-scalar-native-stack-comparison-20261004-v2.service',
         '-o','json','--no-pager'],text=True).splitlines()]
    rows = [r for r in rows if execution['invocation_id'] in
            [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    journal = ep.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as handle:
        for row in rows:handle.write(json.dumps(row,sort_keys=True)+'\n')
    bind(pins,journal); verify(pins)
    result = dict(status='verified_terminal_stack_comparison_native_sigkill_incomplete',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_tool_session_id=99691,
        original_debugger_wrapper_wait_exit_code=0,debugger_exit_code=0,
        original_wrapper_transport=transport,native_signal_from_original_debugger_text='SIGKILL',
        native_exit_code_recorded_by_original_callback=None,callback_did_not_capture_kill_signal=True,
        last_saved_iteration=9,planned_iterations=20,planned_horizon_completed=False,
        partial_scalar_readback=scalar,native_cpu_hard_cap_seconds=5674,
        termination_consistent_with_native_cpu_hard_cap=True,termination_cause_independently_proven=False,
        stack_soft_bytes=64*2**20,stack_hard_unlimited=True,
        standalone_stack_change_validated_as_fix=False,source_hashes=pins,
        original_jobs_restarted=False,production_or_global_limits_changed=False,
        scientific_eligibility=False,posterior_qualified=False,gpu=False,
        scope='Actual original wait and whole wrapper/manager payloads are verified; native '
        'SIGKILL is explicit in the unchanged original GDB text. Partial scalar rows0..9 '
        'pass independent paired readback, but iteration20 never completes. The native '
        'CPU cap and prior live caps are retained; CPU-cap attribution is an inference, '
        'not a separately captured kernel cause. GDB exit0 and callback-null signal/code '
        'must not be treated as native success. No posterior or repaired production run.')
    with output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','original_wrapper_transport']},indent=2))


if __name__=='__main__':main()
