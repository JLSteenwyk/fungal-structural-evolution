#!/usr/bin/env python3
"""Preserve original software waits and compare whole wrapper journal payloads."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execution',type=Path,required=True)
    p.add_argument('--unit',required=True)
    p.add_argument('--tool-payloads',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--expected-exit',type=int,required=True)
    a=p.parse_args();assert not a.output.exists()
    e=json.loads(a.execution.read_text());payload=json.loads(a.tool_payloads.read_text())
    assert e['exit_code']==a.expected_exit and not e['timed_out']
    assert payload['terminal']['exit_code']==a.expected_exit
    assert payload['initial']['session_id']>0 and e['invocation_id'] in payload['initial']['output']
    assert ('code=exited/status='+str(a.expected_exit)) in payload['terminal']['output']
    receipt=Path(e['receipt']);assert sha(receipt)==e['receipt_sha256']
    pins={}
    for mapping in [e['source_hashes'],e['artifacts']]:
        for path,digest in mapping.items():bind(pins,path,digest)
    if a.expected_exit==0:
        r=json.loads(receipt.read_text())
        for path,digest in r['source_hashes'].items():bind(pins,path,digest)
    rows=[json.loads(line) for line in subprocess.check_output(['journalctl','--user','-u',a.unit,'-o','json','--no-pager'],text=True).splitlines()]
    inv=e['invocation_id'];rows=[r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    exact=[r for r in rows if r.get('_PID')==str(e['wrapper']['pid']) and r.get('_CMDLINE')==' '.join(e['wrapper']['cmdline'])]
    assert len(exact)==2 and all(r['_SYSTEMD_INVOCATION_ID']==inv for r in exact)
    assert json.loads(exact[0]['MESSAGE'])==dict(original_wrapper=e['wrapper'],invocation_id=inv)
    terminal={k:v for k,v in e.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}
    assert json.loads(exact[1]['MESSAGE'])==terminal
    starts=[r for r in rows if r.get('USER_INVOCATION_ID')==inv and 'Started ' in r.get('MESSAGE','')]
    ends=[r for r in rows if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts)==len(ends)==1
    journal=a.execution.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as f:
        for row in rows:f.write(json.dumps(row,sort_keys=True)+'\n')
    for path in [Path(__file__),a.execution,a.tool_payloads,receipt,journal]:bind(pins,path)
    verify(pins)
    result=dict(status='verified_original_weighted_parallel_software_wait_zero' if a.expected_exit==0 else 'verified_original_weighted_parallel_fixture_report_failure_retained',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_tool_session_id=payload['initial']['session_id'],
        original_tool_initial_payload=payload['initial'],original_tool_terminal_payload=payload['terminal'],
        original_tool_terminal_exit_code=a.expected_exit,invocation_id=inv,wrapper=e['wrapper'],
        whole_wrapper_initial_and_terminal_payloads_matched=True,manager_start_records=1,
        manager_completion_records=1,source_hashes=pins,scientific_eligibility=False,
        scope='Original software tool wait, wrapper identity, whole wrapper initial/terminal payloads, '
              'source hashes and original manager journal are verified. Report failures remain '
              'retained. Generic wrapper scope is superseded by its exact command; no production '
              'numerical closure, native repair or biological inference.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper','original_tool_initial_payload','original_tool_terminal_payload']},indent=2))


if __name__=='__main__':main()
