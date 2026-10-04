#!/usr/bin/env python3
"""Match a bounded software receipt to its original terminal journal payload."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from full_weighted_fit_exports import atomic
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execution',type=Path,required=True);p.add_argument('--transport',type=Path,required=True)
    p.add_argument('--journal',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    e=json.loads(a.execution.read_text());t=json.loads(a.transport.read_text())
    verify(t['source_hashes']);assert t['source_hashes'][str(a.execution)]==sha(a.execution)
    assert t['source_hashes'][str(a.journal)]==sha(a.journal)
    assert t['wrapper']==e['wrapper'] and t['invocation_id']==e['invocation_id']
    assert t['actual_tool_terminal_exit_code']==e['exit_code']==0 and e['status']=='exited_zero_with_receipt' and not e['timed_out']
    assert sha(e['receipt'])==e['receipt_sha256']==t['validation_sha256']
    rows=[json.loads(l) for l in a.journal.read_text().splitlines()]
    exact=[r for r in rows if r.get('_PID')==str(e['wrapper']['pid']) and r.get('_CMDLINE')==' '.join(e['wrapper']['cmdline'])]
    assert len(exact)==2 and all(r['_SYSTEMD_INVOCATION_ID']==e['invocation_id'] for r in exact)
    initial=json.loads(exact[0]['MESSAGE']);assert initial==dict(original_wrapper=e['wrapper'],invocation_id=e['invocation_id'])
    expected={k:v for k,v in e.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}
    terminal=json.loads(exact[1]['MESSAGE']);assert terminal==expected
    bindings=dict(t['source_hashes'])
    for path in [a.transport,a.execution,a.journal,Path(__file__)]:bind(bindings,path)
    verify(bindings)
    result=dict(status='verified_exact_original_software_terminal_journal_payload',checked_utc=datetime.now(timezone.utc).isoformat(),
        execution=str(a.execution),transport=str(a.transport),actual_tool_session_id=t['actual_tool_session_id'],
        actual_tool_terminal_exit_code=0,wrapper=e['wrapper'],invocation_id=e['invocation_id'],
        original_terminal_receipt_sha256=terminal['receipt_sha256'],entire_terminal_payload_matched=True,
        source_hashes=bindings,scientific_eligibility=False,
        scope='Original wrapper initial identity and entire original terminal stdout payload match saved execution, including receipt hash/status/exit/caps and timing. Complements original manager start/completion and actual tool-wait proof; no additional scientific execution or acceptance.')
    atomic(a.output,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper']},indent=2))


if __name__=='__main__':main()
