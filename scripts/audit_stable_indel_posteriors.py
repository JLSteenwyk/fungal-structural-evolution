#!/usr/bin/env python3
"""Check every produced gap marginal with independent log-space messages."""
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
from logspace_indel_posteriors import infer
from prepare_case_ancestral_neighborhoods import sha


def audit(job,path):
    r=json.loads(path.read_text())
    for p,h in job['pins'].items():
        assert sha(p)==h
    assert r['source_receipt_sha256']==sha(Path(job['output'])/'receipt.json')
    for name,digest in r['artifacts'].items():
        assert sha(path.parent/name)==digest
    result=dict(job_id=job['job_id'],source_receipt_sha256=sha(path))
    if not job['character_count']:
        assert r['status']=='no_coded_characters_no_posterior'
        result['status']='empty_input_verified'
        return result
    fitted=json.loads((Path(job['output'])/'receipt.json').read_text())
    assert r['selected_fit']==fitted['best']
    nodes=json.loads((path.parent/'nodes.json').read_text())
    tree=Phylo.read(job['tree'],'newick')
    assert len(list(tree.find_clades()))==len(nodes)
    for node,row in zip(tree.find_clades(),nodes):
        assert node.name==row['original_name'] and node.is_terminal()==row['is_tip']
        node.name=row['name']
    seq={s.id:str(s.seq) for s in SeqIO.parse(job['characters'],'fasta')}
    g,l,a=r['selected_fit']['parameters']
    rates=4*np.diff(gammainc(a+1,a*gamma.ppf([0,.25,.5,.75,1],a=a,scale=1/a)))
    probs,site=infer(tree,seq,g,l,rates)
    with np.load(path.parent/'posterior.npz',allow_pickle=False) as output:
        assert list(output['node_names'])==[x['name'] for x in nodes]
        assert np.array_equal(output['character'],np.arange(1,job['character_count']+1))
        observed=output['probability_gap']
        expected=np.array([probs[x['name']] for x in nodes])
        assert observed.shape==expected.shape
        error=float(np.max(abs(observed-expected)))
        ll_error=float(np.max(abs(site-output['unconditioned_site_log_likelihood'])))
        assert error<1e-9 and ll_error<1e-8,(job['job_id'],error,ll_error)
        result.update(status='all_node_gap_probabilities_independently_verified',
            probabilities=int(observed.size),maximum_probability_error=error,
            maximum_site_log_likelihood_error=ll_error)
    return result


def main():
    pp=Path('metadata/stable_indel_posterior_audit_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():
        assert sha(p)==h
    jobs=json.loads(Path(plan['fit_plan']).read_text())['jobs']
    identity=json.loads(Path(plan['producer_launch']).read_text())
    root=Path(plan['posterior_output']);out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    done={};terminal_seen=False
    with threadpool_limits(limits=1):
        while True:
            for job in jobs:
                if job['job_id'] in done:
                    continue
                p=root/job['job_id']/'receipt.json'
                if not p.exists():
                    continue
                try:
                    json.loads(p.read_text())
                except json.JSONDecodeError:
                    continue
                r=audit(job,p);done[job['job_id']]=r
                (out/(job['job_id']+'.json')).write_text(json.dumps(r,indent=2)+'\n')
                print(job['job_id'],r['status'],flush=True)
            state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',identity['unit'],'--property=ActiveState,MainPID,Result,ExecMainStatus'],text=True).splitlines())
            if state['ActiveState']=='inactive':
                assert state['Result']=='success' and state['ExecMainStatus']=='0',state
                if len(done)<len(jobs) and not terminal_seen:
                    terminal_seen=True;continue
                assert len(done)==len(jobs);break
            assert state['ActiveState']=='active',state
            p=psutil.Process(identity['pid'])
            assert p.create_time()==identity['created'] and p.cmdline()==identity['cmdline'] and int(state['MainPID'])==p.pid
            print('waiting_for_verified_posterior_producer',len(done),flush=True);time.sleep(30)
    rp=root/'receipt.json';source=json.loads(rp.read_text())
    for jid,r in done.items():
        assert source['job_receipts'][jid+'/receipt.json']==r['source_receipt_sha256']
    nonempty=[r for r in done.values() if 'probabilities' in r]
    receipt=dict(status='all_selected_fit_gap_probabilities_independently_verified',
        input_dispositions=len(done),nonempty_models=len(nonempty),
        probabilities=sum(r['probabilities'] for r in nonempty),
        maximum_probability_error=max(r['maximum_probability_error'] for r in nonempty),
        maximum_site_log_likelihood_error=max(r['maximum_site_log_likelihood_error'] for r in nonempty),
        producer_terminal_state=state,source_receipt_sha256=sha(rp),plan_sha256=sha(pp),
        artifacts={p.name:sha(p) for p in out.glob('*.json')},scope=plan['scope'])
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
