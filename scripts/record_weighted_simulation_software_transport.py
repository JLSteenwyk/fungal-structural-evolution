#!/usr/bin/env python3
"""Bind the original simulation software tool wait to invocation-specific journals."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from full_weighted_fit_exports import atomic
from reference_measurement_union_sources import bind, verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--actual-tool-session',type=int,required=True);p.add_argument('--actual-tool-exit',type=int,required=True)
    a=p.parse_args();assert a.actual_tool_exit==0
    prefix='metadata/weighted_shared_entity_simulation_software_';suffix='20261004_v1.json'
    ep=Path(prefix+'execution_'+suffix);gp=Path(prefix+'validation_'+suffix)
    e=json.loads(ep.read_text());v=json.loads(gp.read_text())
    assert e['status']=='exited_zero_with_receipt' and e['exit_code']==0 and not e['timed_out']
    assert e['receipt_sha256']==sha(gp) and v['fits_computed']==0 and v['scientific_eligibility'] is False
    unit='fungal-weighted-shared-entity-simulation-software-20261004-v1.service'
    rows=[json.loads(l) for l in subprocess.check_output(['journalctl','--user','-u',unit,'-o','json','--no-pager'],text=True).splitlines()]
    inv=e['invocation_id'];rows=[r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    exact=[r for r in rows if r.get('_PID')==str(e['wrapper']['pid']) and r.get('_CMDLINE')==' '.join(e['wrapper']['cmdline'])]
    assert len(exact)==2 and all(r['_SYSTEMD_INVOCATION_ID']==inv for r in exact)
    first=json.loads(exact[0]['MESSAGE']);assert first['original_wrapper']==e['wrapper'] and first['invocation_id']==inv
    starts=[r for r in rows if r.get('USER_INVOCATION_ID')==inv and 'Started ' in r.get('MESSAGE','')]
    ends=[r for r in rows if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts)==len(ends)==1
    journal=ep.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as f:
        for r in rows:f.write(json.dumps(r,sort_keys=True)+'\n')
    bindings={}
    for mapping in [e['source_hashes'],e['artifacts'],v['source_hashes']]:
        for path,d in mapping.items():bind(bindings,path,d)
    for path in [ep,gp,journal,Path(__file__),'metadata/weighted_shared_entity_simulation_environment_20261004_v1.json']:bind(bindings,path)
    verify(bindings)
    proof=dict(status='verified_original_weighted_shared_entity_simulation_software_actual_wait_zero',
        checked_utc=datetime.now(timezone.utc).isoformat(),actual_tool_session_id=a.actual_tool_session,
        actual_tool_terminal_exit_code=a.actual_tool_exit,validation_sha256=sha(gp),unit=unit,invocation_id=inv,
        wrapper=e['wrapper'],exact_wrapper_pid_journal_entries=len(exact),original_start_records=len(starts),
        original_completion_records=len(ends),source_hashes=bindings,
        scope='Actual original software wait, exact wrapper PID/create/command and invocation-specific start/end journals. All declared generator/parent/fixture sources hashed, plus a separate post-execution installed-environment observation. Gaussian software responses are neither native refits nor a global joint-null/biological calibration. No real fit or full-data simulation launched.')
    atomic(Path(prefix+'transport_'+suffix),proof)
    print(json.dumps({k:v for k,v in proof.items() if k not in ['source_hashes','wrapper']},indent=2))


if __name__=='__main__':main()
