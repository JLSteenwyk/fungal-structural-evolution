#!/usr/bin/env python3
"""Close every original failed-role comparison, retaining all special-value reviews."""
import argparse
import ast
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists();pins={};journals=[]
    root=Path('results/ancestral/baliphy-joint-fasta-v7-all-failed-roles-20261004-v1')
    producer=json.loads((root/'receipt.json').read_text());reader=json.loads((root/'readback.json').read_text())
    assert producer['status']=='complete_all24_v7_failure_comparison_dispositions_pending_readback'
    assert reader['status']=='passed_independent_all24_v7_failure_comparison_readback'
    assert reader['producer_receipt_sha256']==sha(root/'receipt.json')
    assert producer['summary']=={k:reader[k] for k in producer['summary']}
    assert producer['roles']==reader['roles']==24 and producer['native_zero_exit_roles']==reader['native_zero_exit_roles']==24
    assert producer['integrity_checked_roles']==reader['integrity_checked_roles']==22
    assert not producer['full_original_failure_grid_output_integrity_checked'] and not reader['posterior_qualified']
    for prefix,session,unit,rp in [
        ('baliphy_joint_fasta_v7_failure_grid',72605,'fungal-joint-fasta-v7-failure-grid-20261004-v1.service',root/'receipt.json'),
        ('baliphy_joint_fasta_v7_failure_grid_reader',84274,'fungal-joint-fasta-v7-failure-grid-reader-20261004-v1.service',root/'readback.json')]:
        ep=Path('metadata',prefix+'_execution_20261004_v1.json')
        tp=Path('metadata',prefix+'_transport_20261005_v1.json')
        pp=Path('metadata',prefix+'_original_tool_payloads_20261005_v1.json')
        e,t,tool=[json.loads(q.read_text()) for q in [ep,tp,pp]]
        assert e['exit_code']==0 and not e['timed_out'] and e['receipt_sha256']==sha(rp)==t['validation_sha256']
        assert t['original_tool_session_id']==tool['original_tool_session_id']==session
        assert tool['initial']['session_id']==session and tool['terminal']['exit_code']==0
        assert t['whole_wrapper_initial_and_terminal_payloads_matched'] and t['manager_start_records']==t['manager_completion_records']==1
        for mapping in [e['source_hashes'],e['artifacts'],t['source_hashes'],json.loads(rp.read_text())['source_hashes']]:
            for q,h in mapping.items():bind(pins,q,h)
        for q in [ep,tp,pp,rp]:bind(pins,q)
        text=subprocess.check_output(['journalctl','--user','-u',unit,'--all','-o','json','--no-pager'],text=True)
        rows=[json.loads(line) for line in text.splitlines()]
        inv=e['invocation_id'];rows=[r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
        journal=Path('metadata',prefix+'_original_whole_journal_20261005_v1.jsonl')
        with journal.open('x') as f:
            for r in rows:f.write(json.dumps(r,sort_keys=True)+'\n')
        bind(pins,journal)
        exact=[r for r in rows if r.get('_PID')==str(e['wrapper']['pid']) and r.get('_CMDLINE')==' '.join(e['wrapper']['cmdline'])]
        assert len(exact)==2
        assert json.loads(exact[0]['MESSAGE'])==dict(original_wrapper=e['wrapper'],invocation_id=inv)
        assert json.loads(exact[1]['MESSAGE'])=={k:v for k,v in e.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}
        native_log=ep.with_suffix('')/'stdout.log';last=[line for line in native_log.read_text().splitlines() if line.strip()][-1]
        assert json.loads(last)==producer['summary']
        queue_verified=False
        if session==84274:
            lp=Path('metadata/baliphy_joint_fasta_v7_failure_grid_reader_launch_20261004_v1.json');launch=json.loads(lp.read_text())
            wait=Path(launch['plan']);plan=json.loads(wait.read_text());assert sha(wait)==launch['plan_sha256']
            assert launch['invocation_id']==inv and launch['original_tool_session_id']==session
            outer=[r for r in rows if r.get('_PID')==str(launch['pid']) and r.get('_CMDLINE')==' '.join(launch['cmdline'])]
            commands=[r['MESSAGE'] for r in outer if r.get('MESSAGE','').startswith('running_original_journal_verified_command ')]
            assert len(commands)==1 and ast.literal_eval(commands[0].split(' ',1)[1])==plan['command']
            assert any(r.get('MESSAGE','').startswith('waiting_for_original_dependency_completion ') for r in outer)
            for q,h in plan['pins'].items():bind(pins,q,h)
            for q in [lp,wait]:bind(pins,q)
            queue_verified=True
        assert len([r for r in rows if r.get('USER_INVOCATION_ID')==inv and 'Started ' in r.get('MESSAGE','')])==1
        assert len([r for r in rows if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')])==1
        assert not any('Failed with result' in r.get('MESSAGE','') or 'Main process exited' in r.get('MESSAGE','') for r in rows)
        journals.append(dict(unit=unit,invocation_id=inv,original_tool_session_id=session,original_exit_code=0,
            native_exit_code=0,whole_wrapper_and_native_summary_verified=True,queued_outer_controller_verified=queue_verified,
            journal=str(journal),journal_sha256=sha(journal)))
    dispositions=json.loads((root/'dispositions.json').read_text());assert len(dispositions)==24
    reviews=[r for r in dispositions if r['status']=='explicit_nonfinite_or_literal_null_scalar_output_retained_for_review']
    assert len(reviews)==2 and all(r['scalar_integrity_accepted'] is False and not r['joint_frames'] for r in reviews)
    special=[v for r in reviews for v in r['scalar_v6_audit']['nonfinite_reviews']]
    assert len(special)==8 and all(v['kind']=='positive_infinity' and v['path']==['S1/','ASRV.Gamma:alpha'] for v in special)
    assert all(not r['scalar_v6_audit']['literal_null_iterations'] for r in reviews)
    bind(pins,Path(__file__));verify(pins)
    result=dict(status='complete_verified_all24_v7_failure_comparison_dispositions',checked_utc=datetime.now(timezone.utc).isoformat(),
        original_roles=24,native_zero_exits=24,fully_finite_integrity_roles=22,explicit_special_value_review_roles=2,
        complete_integrity_quartets=4,prior_quartets=6,scalar_mapped_values=19866,accepted_joint_frames=66,
        accepted_ancestral_residue_category_pairs=21017330,positive_infinity_observations=8,
        special_parameter=['S1/','ASRV.Gamma:alpha'],original_journal_closures=journals,
        source_hashes=pins,all_bound_sources_and_outputs=len(pins),full_all24_output_integrity_accepted=False,
        installed_formatter_repaired=False,scientific_eligibility=False,posterior_qualified=False,new_mcmc_runs=0,
        gpu=False,all_eight_aims_incomplete=True,
        scope='All24 original failed-role comparisons and every disposition independently replayed and closed against '
              'actual original waits, complete wrapper/native summaries, queued outer controller and manager start/end. '
              'All native exits zero;22fully finite output checks and two explicit special-value reviews retained. '
              'Eight positive-infinity ASRV alpha observations remain excluded from accepted arrays. Original24SIGSEGV '
              'attempts and all seeds/priors/input/caps unchanged; this closes comparison custody, not a global crash '
              'or formatter repair, sufficient posterior, chosen prior, model calibration or completed biological aim.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes']},indent=2))


if __name__=='__main__':main()
