#!/usr/bin/env python3
"""Retain exact original software wait, complete wrapper payload and native artifacts."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execution',type=Path,required=True);p.add_argument('--validation',type=Path,required=True)
    p.add_argument('--native-root',type=Path,required=True);p.add_argument('--unit',required=True)
    p.add_argument('--actual-tool-session',type=int,required=True);p.add_argument('--actual-tool-exit',type=int,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
    e=json.loads(a.execution.read_text());assert e['exit_code']==a.actual_tool_exit and not e['timed_out']
    pins={}
    for d in [e['source_hashes'],e['artifacts']]:
        for n,h in d.items():bind(pins,n,h)
    if a.actual_tool_exit==0:
        assert e['status']=='exited_zero_with_receipt' and e['receipt_sha256']==sha(a.validation)
        v=json.loads(a.validation.read_text());verify(v['source_hashes'])
        assert v['status']=='passed_scalar_json_v6_full_source_and_paired_native_qualification'
        assert (v['full_programs_checked'],v['full_roles_checked'],v['native_runs'],v['scalar_rows_checked'])==(405,1620,7,63)
        assert (v['native_finite_value_roundtrips'],v['native_nonfinite_tags_checked'],v['native_literal_null_tags_checked'])==(18,6,2)
        assert len(v['malformed_record_cases_rejected_by_both_readers'])==22
        assert v['scientific_eligibility'] is v['posterior_qualified'] is False
        for n,h in v['source_hashes'].items():bind(pins,n,h)
        bind(pins,a.validation)
        status='verified_original_scalar_v6_software_wait_zero'
    else:
        assert a.actual_tool_exit==1 and e['status']=='failed_stage_retained'
        assert e['receipt_sha256'] is None and not a.validation.exists()
        status='retained_original_failed_scalar_v6_software_attempt'
    # Snapshot every original nested native receipt/configuration/output;
    # native receipt contents establish the files existed before child exit.
    receipts=list(a.native_root.rglob('receipt.json'));assert len(receipts)==(7 if a.actual_tool_exit==0 else 1)
    for path in receipts:
        r=json.loads(path.read_text());bind(pins,path)
        config=path.parent.parent/'configuration.json';bind(pins,config)
        for n,h in json.loads(config.read_text())['pins'].items():bind(pins,n,h)
        for n,h in r['artifacts'].items():bind(pins,path.parent/n,h)
    rows=[json.loads(l) for l in subprocess.check_output(['journalctl','--user','-u',a.unit,'-o','json','--no-pager'],text=True).splitlines()]
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
        for r in rows:f.write(json.dumps(r,sort_keys=True)+'\n')
    for path in [a.execution,journal,Path(__file__)]:bind(pins,path)
    verify(pins)
    result=dict(status=status,checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_tool_session_id=a.actual_tool_session,actual_tool_terminal_exit_code=a.actual_tool_exit,
        invocation_id=inv,unit=a.unit,wrapper=e['wrapper'],exact_wrapper_pid_journal_entries=2,
        original_start_records=1,original_completion_records=1,entire_terminal_payload_matched=True,
        source_hashes=pins,native_receipts_retained=len(receipts),gpu=False,scientific_eligibility=False,
        scope='Original actual tool wait and exact manager/wrapper PID/create/command/invocation custody. '
              'The reused generic wrapper scope refers to an older kernel checker; its actual command '
              'and this receipt identify the scalar V6 software check. Original failed attempts and '
              'their byte-identical programs remain retained, never restarted or overwritten.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper']},indent=2))


if __name__=='__main__':main()
