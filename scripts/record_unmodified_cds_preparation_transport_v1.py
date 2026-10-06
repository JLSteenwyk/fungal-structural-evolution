#!/usr/bin/env python3
"""Close the actual full original-CDS translation source preparation and immutable plan."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from full_weighted_fit_exports import atomic
from reference_measurement_union_sources import bind, verify


def run(execution_path, validation_path, payload_path, unit, output):
    ep, vp, pp = map(Path, [execution_path, validation_path, payload_path])
    e, v, payload = [json.loads(p.read_text()) for p in [ep, vp, pp]]
    assert e['status']=='exited_zero_with_receipt' and e['exit_code']==0 and not e['timed_out']
    assert e['receipt_sha256']==sha(vp) and v['scientific_eligibility'] is False
    assert payload['initial']['session_id']==payload['original_tool_session_id']
    assert payload['terminal']['exit_code']==0
    assert v['status'] == 'prepared_full526_unmodified_original_cds_fixed_code_translation'
    inv = e['invocation_id']
    assert inv in payload['initial']['output'] and unit in payload['initial']['output']
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl','--user','-u',unit,'--all','-o','json','--no-pager'],text=True).splitlines()]
    rows = [r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    exact = [r for r in rows if r.get('_PID')==str(e['wrapper']['pid'])
             and r.get('_CMDLINE')==' '.join(e['wrapper']['cmdline'])]
    assert len(exact)==2
    assert json.loads(exact[0]['MESSAGE'])==dict(original_wrapper=e['wrapper'],invocation_id=inv)
    terminal = {k:v for k,v in e.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}
    assert json.loads(exact[1]['MESSAGE'])==terminal
    starts = [r for r in rows if r.get('USER_INVOCATION_ID')==inv and 'Started ' in r.get('MESSAGE','')]
    ends = [r for r in rows if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts)==len(ends)==1
    journal = ep.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as handle:
        for r in rows:handle.write(json.dumps(r,sort_keys=True)+'\n')
    pins = {}
    for mapping in [e['source_hashes'],e['artifacts'],v['pins']]:
        for p,h in mapping.items():bind(pins,p,h)
    for p in [ep,vp,pp,journal,Path(__file__)]:bind(pins,p)
    verify(pins)
    result = dict(status='verified_original_software_wait_and_whole_wrapper_payloads',
        checked_utc=datetime.now(timezone.utc).isoformat(),unit=unit,invocation_id=inv,wrapper=e['wrapper'],
        original_tool_session_id=payload['original_tool_session_id'],original_tool_terminal_exit_code=0,
        actual_tool_session_id=payload['original_tool_session_id'],actual_tool_terminal_exit_code=0,
        validation_sha256=sha(vp),whole_wrapper_initial_and_terminal_payloads_matched=True,
        manager_start_records=1,manager_completion_records=1,exact_wrapper_pid_journal_entries=2,
        original_start_records=1,original_completion_records=1,source_hashes=pins,
        scientific_eligibility=False,scope='Exact original software wait and complete wrapper payloads, '
        'original invocation manager start/end and all declared software/output source bindings. '
        'This closes complete original-CDS translation source preparation, not translation, genetic-code adoption or biological analysis.')
    atomic(Path(output),result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper']},indent=2))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['execution','validation','payloads','unit','output']:p.add_argument('--'+name,required=True)
    a=p.parse_args();run(a.execution,a.validation,a.payloads,a.unit,a.output)
