#!/usr/bin/env python3
"""Verify original V6 sampler-controller software wait and complete terminal payload."""
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
    assert v['status']=='passed_full_scalar_v6_sampler_controller_software_v1'
    assert (v['full_roles'],v['full_quartets'],v['retained_actual_native_roles'],
            v['retained_actual_scalar_rows'],v['retained_actual_mapped_values'],v['retained_actual_joint_frames'])==(1620,405,3,63,2709,9)
    assert v['full_job_file_matches_qualified_admission_adapter'] and v['inherited_design_rejections']==12
    assert v['synthetic_full_controller_serialization_checked'] and v['synthetic_checked_roles']==1617
    assert v['synthetic_unresolved_roles']==3 and v['synthetic_unresolved_quartets']==3
    assert v['synthetic_scalar_review_roles']==1 and v['synthetic_joint_frames']==4851
    assert len(v['serialization_cases_rejected'])==len(v['synthetic_startup_metadata_cases_rejected'])==11
    assert len(v['abort_admission_cases_checked'])==2 and v['missing_actual_export_readbacks_rejected']==3
    assert v['actual_prerequisite_gate'] in ['actual_unclosed_v6_startup_refused_before_native_admission',
                                           'all_original_prerequisite_closures_verified']
    assert v['new_native_runs']==0
    assert all(v[k] is False for k in ['full_grid_sampler_launched','native_resource_observer_qualified',
                                     'posterior_qualified','scientific_eligibility','gpu'])
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
    result=dict(status='verified_original_scalar_v6_sampler_controller_software_wait_zero',checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_tool_session_id=a.actual_tool_session,actual_tool_terminal_exit_code=0,
        invocation_id=inv,unit=a.unit,wrapper=e['wrapper'],validation_sha256=sha(a.validation),
        entire_terminal_payload_matched=True,original_start_records=1,original_completion_records=1,
        source_hashes=pins,new_native_runs=0,new_mcmc_runs=0,scientific_eligibility=False,gpu=False,
        scope='Original capped V6 sampler-controller software wait; exact wrapper PID/create/command/'
              'invocation, whole initial/terminal payloads and manager start/end verified. All '
              'source/artifact bindings freshly rehashed. Actual native fixtures are inherited; '
              'complete1620serialization explicitly mocks native/admission/gates/caps and retains '
              'invalid/review/failure outcomes. Generic wrapper kernel scope is superseded by its '
              'actual controller-checker command. No future native sampler/resource observation, '
              'posterior qualification, original restart or biological acceptance.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper']},indent=2))


if __name__=='__main__':main()
