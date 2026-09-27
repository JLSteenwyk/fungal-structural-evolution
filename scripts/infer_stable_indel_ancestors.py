#!/usr/bin/env python3
"""Propagate selected audited gap models to every fitted-tree node."""
import hashlib
import json
import subprocess
import time
from pathlib import Path
import numpy as np
import psutil
from Bio import Phylo,SeqIO
from scipy.special import gammainc
from scipy.stats import gamma
from threadpoolctl import threadpool_limits
from stable_indel_posteriors import infer
from map_ancestral_fitted_nodes import signatures
from prepare_case_ancestral_neighborhoods import sha


def produce(job,audit_path,out):
    audit=json.loads(audit_path.read_text())
    source=Path(job['output'])/'receipt.json'
    assert sha(source)==audit['source_receipt_sha256']
    fit=json.loads(source.read_text());assert fit['job']==job
    for p,h in job['pins'].items():
        assert sha(p)==h,p
    folder=out/job['job_id'];folder.mkdir()
    result=dict(job_id=job['job_id'],source_receipt_sha256=sha(source),audit_sha256=sha(audit_path))
    if not job['character_count']:
        assert audit['status']=='empty_input_verified'
        result['status']='no_coded_characters_no_posterior'
    else:
        assert audit['status']=='all_nine_candidate_likelihoods_independently_verified',job['job_id']
        tree=Phylo.read(job['tree'],'newick');node_rows=[]
        sigs=signatures(tree)
        for i,node in enumerate(tree.find_clades()):
            original=node.name
            if not node.is_terminal():
                node.name=f'indel_internal_{i}'
            row=dict(name=node.name,original_name=original,is_tip=node.is_terminal())
            if node in sigs:
                row['incident_tip_partition_sha256']=hashlib.sha256(json.dumps(sigs[node],separators=(',',':')).encode()).hexdigest()
                row['incident_edge_count']=len(sigs[node])
            node_rows.append(row)
        seq={r.id:str(r.seq) for r in SeqIO.parse(job['characters'],'fasta')}
        g,l,a=fit['best']['parameters']
        rates=4*np.diff(gammainc(a+1,a*gamma.ppf([0,.25,.5,.75,1],a=a,scale=1/a)))
        probs,logl=infer(tree,seq,g,l,rates)
        names=[r['name'] for r in node_rows]
        values=np.array([probs[n] for n in names])
        assert values.shape==(len(names),job['character_count'])
        assert np.isfinite(values).all() and values.min()>=0 and values.max()<=1+1e-12
        max_tip_error=0.
        for n,s in seq.items():
            for col,state in enumerate(s):
                if state!='?':
                    max_tip_error=max(max_tip_error,abs(probs[n][col]-int(state)))
        assert max_tip_error<1e-12
        np.savez_compressed(folder/'posterior.npz',probability_gap=values,node_names=np.array(names),
            character=np.arange(1,job['character_count']+1),unconditioned_site_log_likelihood=logl)
        (folder/'nodes.json').write_text(json.dumps(node_rows,indent=2)+'\n')
        result.update(status='conditional_gap_posteriors_produced_pending_independent_audit',
            nodes=len(names),characters=job['character_count'],probabilities=int(values.size),
            maximum_observed_tip_error=max_tip_error,selected_fit=fit['best'],
            selected_source_start=fit['selected_source_start'],
            optimization_diagnostics={k:audit[k] for k in ['best_optimizer_success','best_near_bound','start_likelihood_spread']},
            correction=job['correction'],terminal_policy=job['terminal_policy'])
    result['artifacts']={p.name:sha(p) for p in folder.iterdir() if p.is_file()}
    (folder/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    pp=Path('metadata/stable_indel_ancestor_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():
        assert sha(p)==h,p
    source_plan=json.loads(Path(plan['fit_plan']).read_text())
    identity=json.loads(Path(plan['auditor_launch']).read_text())
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    done={};terminal_seen=False
    with threadpool_limits(limits=1):
        while True:
            for job in source_plan['jobs']:
                if job['job_id'] in done:
                    continue
                ap=Path(plan['audit_output'])/(job['job_id']+'.json')
                if not ap.exists():
                    continue
                try:
                    json.loads(ap.read_text())
                except json.JSONDecodeError:
                    continue
                result=produce(job,ap,out);done[job['job_id']]=result
                print(job['job_id'],result['status'],flush=True)
            state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',identity['unit'],'--property=ActiveState,MainPID,Result,ExecMainStatus'],text=True).splitlines())
            if state['ActiveState']=='inactive':
                assert state['Result']=='success' and state['ExecMainStatus']=='0',state
                if len(done)<len(source_plan['jobs']) and not terminal_seen:
                    terminal_seen=True;continue
                assert len(done)==len(source_plan['jobs']);break
            assert state['ActiveState']=='active',state
            p=psutil.Process(identity['pid'])
            assert p.create_time()==identity['created'] and p.cmdline()==identity['cmdline'] and int(state['MainPID'])==p.pid
            print('waiting_for_verified_auditor',len(done),flush=True);time.sleep(30)
    audit_receipt=Path(plan['audit_output'])/'receipt.json'
    ar=json.loads(audit_receipt.read_text())
    assert ar['starts_outside_tolerance']==0 and ar['starts_checked']==2754
    for r in done.values():
        assert ar['artifacts'][r['job_id']+'.json']==r['audit_sha256']
    result=dict(status='all_selected_fit_gap_posteriors_produced_pending_independent_validation',
        input_dispositions=len(done),nonempty_models=sum('probabilities' in r for r in done.values()),
        probabilities=sum(r.get('probabilities',0) for r in done.values()),
        plan_sha256=sha(pp),audit_receipt_sha256=sha(audit_receipt),auditor_terminal_state=state,
        job_receipts={str(p.relative_to(out)):sha(p) for p in out.glob('*/receipt.json')},scope=plan['scope'])
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
