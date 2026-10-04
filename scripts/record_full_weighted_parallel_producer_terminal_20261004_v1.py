#!/usr/bin/env python3
"""Verify the original complete numerical producer without claiming readback."""
import ast
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from full_weighted_covariance_qualification import PRODUCER
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from reference_measurement_union_sources import bind, verify


def main():
    output = Path('metadata/full_weighted_parallel_producer_terminal_verified_20261004_v1.json')
    assert not output.exists()
    lp = Path('metadata/full_weighted_parallel_producer_launch_20261004_v1.json')
    launch = json.loads(lp.read_text()); launch['launch'] = str(lp)
    assert fingerprint(launch) is None
    original = json.loads(Path('metadata/full_weighted_parallel_producer_original_tool_payloads_20261004_v1.json').read_text())
    assert original['original_tool_session_id']==original['initial']['session_id']==45225
    assert original['terminal']['exit_code']==0
    observed = journal_terminal(launch)
    plan_path = Path('metadata/full_weighted_parallel_qualification_plan_20261004_v1.json')
    plan = json.loads(plan_path.read_text()); root = Path(plan['output'])
    receipt_path = root/'receipt.json'; receipt = json.loads(receipt_path.read_text())
    assert receipt['status']==PRODUCER and receipt['plan_sha256']==sha(plan_path)
    assert receipt['cohorts']==4340 and receipt['numerical_audit_rows']==5208000
    assert receipt['setting_audit_links']==24883200 and receipt['working_model_fits_computed']==0
    assert receipt['cached_numeric_inputs_preserved_for_all_cohorts'] is True
    assert receipt['scientific_eligibility'] is receipt['nonuniform_weighting_accepted'] is False
    pins = dict(receipt['source_hashes'])
    for name,digest in receipt['artifacts'].items():bind(pins,root/name,digest)
    checkpoints = sorted((root/'checkpoints').glob('*.producer.json'))
    assert len(checkpoints)==4340
    for number,p in enumerate(checkpoints):
        cp = json.loads(p.read_text())
        assert cp['cohort_index']==number and cp['cached_numeric_inputs_preserved'] is True
        assert cp['cached_input_array_bindings_checked']>0
        assert cp['worker_address_space_limit_bytes']==12*2**30
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl','--user','-u',launch['unit'],'-o','json','--no-pager'],text=True).splitlines()]
    inv = launch['invocation_id']; assert inv in original['initial']['output']
    rows = [r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    wrapper = [r for r in rows if r.get('_PID')==str(launch['pid'])
               and r.get('_CMDLINE')==' '.join(launch['cmdline'])]
    headers = [r for r in wrapper if r.get('MESSAGE','').startswith('running_original_journal_verified_command ')]
    assert len(headers)==1
    wait = json.loads(Path(launch['plan']).read_text())
    assert ast.literal_eval(headers[0]['MESSAGE'].split('running_original_journal_verified_command ',1)[1])==wait['command']
    summary = {k:v for k,v in receipt.items() if k not in ['source_hashes','artifacts','scope']}
    terminals = []
    for row in wrapper:
        try:message = json.loads(row['MESSAGE'])
        except (ValueError,TypeError):continue
        if isinstance(message,dict) and message.get('status')==PRODUCER:terminals.append(message)
    assert terminals==[summary]
    starts = [r for r in rows if r.get('USER_INVOCATION_ID')==inv and 'Started ' in r.get('MESSAGE','')]
    ends = [r for r in rows if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts)==len(ends)==1
    journal = Path('metadata/full_weighted_parallel_producer_original_journal_20261004_v1.jsonl')
    with journal.open('x') as handle:
        for row in rows:handle.write(json.dumps(row,sort_keys=True)+'\n')
    for p in [lp,Path(launch['plan']),plan_path,receipt_path,journal,Path(__file__),
              Path('metadata/full_weighted_parallel_producer_original_tool_payloads_20261004_v1.json')]:bind(pins,p)
    verify(pins)
    result = dict(status='verified_original_full_parallel_numerical_producer_only',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_tool_session_id=45225,
        original_tool_terminal_exit_code=0,invocation_id=inv,wrapper=launch,
        whole_original_command_header_and_native_terminal_summary_matched=True,
        original_manager_start_records=1,original_manager_completion_records=1,
        producer_receipt=str(receipt_path),producer_receipt_sha256=sha(receipt_path),
        cohorts=4340,numerical_audit_rows=5208000,setting_audit_links=24883200,
        original_terminal_handle=observed,audit_status_counts=receipt['audit_status_counts'],
        full_declared_bindings_freshly_verified=len(pins),source_hashes=pins,
        independent_readback_completed=False,full_numerical_closure_verified=False,
        scientific_eligibility=False,biological_fits=0,gpu=False,
        scope='All declared original producer source/output files freshly rehashed; all4340 '
        'cached-array checkpoints checked; original wait/header/native-summary/manager journals '
        'matched. This verifies producer termination and arithmetic accounting only. Full '
        'independent numerical readback/closure, actual timing and biological inference remain required.')
    with output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper','original_terminal_handle']},indent=2))


if __name__=='__main__':main()
