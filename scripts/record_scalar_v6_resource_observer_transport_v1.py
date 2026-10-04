#!/usr/bin/env python3
"""Verify original V6 resource-observer software wait and complete terminal payload."""
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
    assert v['status']=='passed_full_scalar_v6_sampler_resource_observer_software_v1'
    assert v['full_roles_checked']==v['full_v6_configurations_checked']==1620
    assert v['historical_qualified_native_probe_readings']>0 and v['actual_live_native_observations']==0
    assert v['full_producer_and_reader_serialization_checked']
    assert v['artificial_failed_roles_retained']==v['artificial_missing_live_observations_retained']==2
    assert len(v['serialization_rejection_cases'])==8 and len(v['design_rejection_cases'])==6
    assert len(v['malformed_or_altered_observations_rejected'])==9
    assert v['new_native_inference_runs']==0 and v['future_native_observer_launched'] is False
    assert v['scientific_eligibility'] is False
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
    result=dict(status='verified_original_scalar_v6_resource_observer_software_wait_zero',checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_tool_session_id=a.actual_tool_session,actual_tool_terminal_exit_code=0,
        invocation_id=inv,unit=a.unit,wrapper=e['wrapper'],validation_sha256=sha(a.validation),
        entire_terminal_payload_matched=True,original_start_records=1,original_completion_records=1,
        source_hashes=pins,new_native_runs=0,new_mcmc_runs=0,scientific_eligibility=False,gpu=False,
        scope='Original capped V6 resource-observer software wait; exact wrapper PID/create/command/'
              'invocation, whole initial/terminal payloads and manager start/end verified. All '
              'source/artifact bindings freshly rehashed. Actual shared procfs decoder fixtures are inherited; '
              'complete1620telemetry serialization explicitly mocks native/capture/journals/caps and retains '
              'missing-observation/failure outcomes. Generic wrapper kernel scope is superseded by its '
              'actual controller-checker command. No future native sampler/resource observation, '
              'posterior qualification, original restart or biological acceptance.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper']},indent=2))


if __name__=='__main__':main()
