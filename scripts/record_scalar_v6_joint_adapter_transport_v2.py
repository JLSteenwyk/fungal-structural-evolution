#!/usr/bin/env python3
"""Verify exact original wait and whole payload for V6 scalar admission software."""
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
    p.add_argument('--actual-tool-session',type=int,required=True);p.add_argument('--actual-tool-exit',type=int,required=True)
    p.add_argument('--unit',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists();e=json.loads(a.execution.read_text())
    assert e['exit_code']==a.actual_tool_exit==0 and not e['timed_out'] and e['status']=='exited_zero_with_receipt'
    assert sha(a.validation)==e['receipt_sha256'];v=json.loads(a.validation.read_text())
    assert v['status']=='passed_v6_scalar_joint_admission_adapter_software_v2'
    assert v['inherited_v3_plan_provider_pins_merged'] is True and v['transitive_project_source_modules']>0
    assert (v['full_source_roles'],v['full_source_models'],v['retained_actual_native_roles'],
            v['retained_actual_scalar_rows'],v['retained_actual_mapped_values'],v['retained_actual_joint_frames'])==(1620,405,3,63,2709,9)
    assert len(v['design_cases_rejected'])==12 and len(v['rehashed_private_output_cases'])==11
    assert v['readback_missing_export_cases_rejected']==3 and v['new_native_runs']==0
    assert all(v[k] is False for k in ['full_grid_native_execution_launched','full_workflow_serialization_qualified',
                                     'full_grid_v6_startup_qualified','scientific_eligibility','posterior_qualified','gpu'])
    pins={}
    for mapping in [e['source_hashes'],e['artifacts'],v['source_hashes']]:
        for n,h in mapping.items():bind(pins,n,h)
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
    for path in [Path(__file__),a.execution,a.validation,journal]:bind(pins,path)
    verify(pins)
    result=dict(status='verified_original_v6_scalar_joint_adapter_v2_software_wait_zero',checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_tool_session_id=a.actual_tool_session,actual_tool_terminal_exit_code=0,
        invocation_id=inv,unit=a.unit,wrapper=e['wrapper'],validation_sha256=sha(a.validation),
        entire_terminal_payload_matched=True,original_start_records=1,original_completion_records=1,
        source_hashes=pins,new_native_runs=0,scientific_eligibility=False,gpu=False,
        scope='Original bounded software wait, exact wrapper PID/create/command/invocation, whole '
              'initial/terminal payloads and manager start/end verified. All source/output bindings '
              'rehash. The generic wrapper scope mentions an older kernel check; its actual command '
              'identifies the V6 scalar admission adapter. Three native fixture attempts are inherited '
              'and unchanged; rehashed negative copies and full-grid summary are explicitly synthetic. '
              'No production startup/controller/closure, native failure repair or posterior acceptance.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper']},indent=2))


if __name__=='__main__':main()
