#!/usr/bin/env python3
"""Bind an observed original tool wait to exact software invocation journals."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from full_weighted_fit_exports import atomic
from reference_measurement_union_sources import verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--version',required=True)
    p.add_argument('--actual-tool-session',required=True,type=int);p.add_argument('--actual-tool-exit',required=True,type=int)
    a=p.parse_args();prefix='metadata/full_weighted_timing_software_';suffix='20261004_'+a.version+'.json'
    ep=Path(prefix+'execution_'+suffix);e=json.loads(ep.read_text());assert e['exit_code']==a.actual_tool_exit
    unit='fungal-full-weighted-timing-software-20261004-'+a.version+'.service'
    rows=[json.loads(line) for line in subprocess.check_output(['journalctl','--user','-u',unit,'-o','json','--no-pager'],text=True).splitlines()]
    invocation=e['invocation_id'];rows=[r for r in rows if invocation in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    exact=[r for r in rows if r.get('_PID')==str(e['wrapper']['pid']) and r.get('_CMDLINE')==' '.join(e['wrapper']['cmdline'])]
    assert len(exact)==2 and all(r['_SYSTEMD_INVOCATION_ID']==invocation for r in exact)
    first=json.loads(exact[0]['MESSAGE']);assert first['original_wrapper']==e['wrapper'] and first['invocation_id']==invocation
    starts=[r for r in rows if r.get('USER_INVOCATION_ID')==invocation and 'Started ' in r.get('MESSAGE','')]
    ends=[r for r in rows if r.get('USER_INVOCATION_ID')==invocation and r.get('CPU_USAGE_NSEC')]
    assert len(starts)==len(ends)==1
    journal=ep.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as f:
        for row in rows:f.write(json.dumps(row,sort_keys=True)+'\n')
    bindings={str(ep):sha(ep),str(journal):sha(journal),str(Path(__file__)):sha(__file__)}
    bindings.update(e['source_hashes']);bindings.update(e['artifacts'])
    success=a.actual_tool_exit==0
    if success:
        assert e['status']=='exited_zero_with_receipt' and not e['timed_out']
        gp=Path(e['receipt']);assert sha(gp)==e['receipt_sha256'];gate=json.loads(gp.read_text())
        assert gate['scientific_eligibility'] is False and gate['fits_computed']==0
        verify(gate['source_hashes']);bindings.update(gate['source_hashes']);bindings[str(gp)]=sha(gp)
    else:assert e['status']=='failed_stage_retained' and e['receipt_sha256'] is None
    verify(bindings)
    proof=dict(status='verified_original_complete_weighted_timing_software_actual_wait_zero' if success else
        'preserved_original_complete_weighted_timing_software_failure',checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_tool_session_id=a.actual_tool_session,actual_tool_terminal_exit_code=a.actual_tool_exit,
        unit=unit,invocation_id=invocation,wrapper=e['wrapper'],exact_wrapper_pid_journal_entries=len(exact),
        original_start_records=len(starts),original_completion_records=len(ends),source_hashes=bindings,
        scope='Actual session/exit supplied from the original tool wait, exact wrapper PID/create/command and original invocation start/end journals. Failure and partial bytes preserved; successful gate bindings freshly verified. No real-data timing, fit, scientific acceptance, restart or resource change.')
    out=Path(prefix+('transport_' if success else 'failure_')+suffix);atomic(out,proof)
    print(json.dumps({k:v for k,v in proof.items() if k not in ['source_hashes','wrapper']},indent=2))


if __name__=='__main__':main()
