#!/usr/bin/env python3
"""Independently replay every completed indel fit while retaining diagnostics."""
import json
import subprocess
import time
from pathlib import Path
import numpy as np
import psutil
from Bio import Phylo, SeqIO
from scipy.special import gammainc
from scipy.stats import gamma
from threadpoolctl import threadpool_limits
from prepare_case_ancestral_neighborhoods import sha
from replay_fastml_indel_probabilities import infer, analytic_check
from diagnose_indel_extreme_rates import high_precision


def audit(job, source):
    receipt = json.loads(source.read_text())
    assert receipt['job'] == job
    for p,h in job['pins'].items():
        assert sha(p) == h,p
    if not job['character_count']:
        assert receipt['status'] == 'no_coded_characters_no_inference'
        return dict(job_id=job['job_id'],status='empty_input_verified',source_receipt_sha256=sha(source))
    tree = Phylo.read(job['tree'],'newick')
    for i,node in enumerate(tree.get_nonterminals()):
        node.name = f'audit_internal_{i}'
    seq = {r.id:str(r.seq) for r in SeqIO.parse(job['characters'],'fasta')}
    width = job['character_count']
    assert len(seq) == job['proteins'] and {len(s) for s in seq.values()} == {width}
    if job['correction'] == 'all_taxa':
        augmented = {n:s+'0' for n,s in seq.items()}
    else:
        assert job['correction'] == 'observed_mask'
        augmented = {n:s+''.join('?' if x=='?' else '0' for x in s) for n,s in seq.items()}
    starts = []
    for index,record in enumerate(receipt['starts']):
        g,l,a = record['parameters']
        assert np.max(abs(np.log(record['parameters'])-record['log_parameters'])) < 1e-12
        rates = 4*np.diff(gammainc(a+1,a*gamma.ppf([0,.25,.5,.75,1],a=a,scale=1/a)))
        _,logp = infer(tree,augmented,g,l,rates)
        excluded = logp[-1] if job['correction']=='all_taxa' else logp[width:]
        corrected = float((logp[:width]-np.log(-np.expm1(excluded))).sum())
        error = abs(corrected-record['log_likelihood'])
        double_error = error
        precision = 'double_expm'
        if error >= 1e-6:
            corrected = float(high_precision(tree,seq,record['parameters'],job['correction']))
            error = abs(corrected-record['log_likelihood'])
            precision = '70_digit_recursive_pruning'
        # Retain genuine discrepancies and finish the entire audit; never label them verified.
        starts.append(dict(start=index,independent_log_likelihood=corrected,absolute_error=error,
            precision=precision,double_expm_error=double_error,within_tolerance=bool(error<1e-6),
            optimizer_success=record['optimizer_success'],near_bound=record['near_bound'],
            improvement_over_start=corrected-record['initial_log_likelihood'],
            reported_gradient=record['gradient']))
    assert len(starts)==3
    best=max(receipt['starts'],key=lambda r:r['log_likelihood'])
    assert receipt['best']==best
    spread=max(r['log_likelihood'] for r in receipt['starts'])-min(r['log_likelihood'] for r in receipt['starts'])
    assert abs(spread-receipt['fitted_log_likelihood_spread'])<1e-12
    return dict(job_id=job['job_id'],status=('all_three_fitted_likelihoods_independently_verified' if all(s['within_tolerance'] for s in starts) else 'likelihood_discrepancies_retained'),
        starts=starts,best_optimizer_success=best['optimizer_success'],
        best_near_bound=best['near_bound'],start_likelihood_spread=spread,
        source_receipt_sha256=sha(source))


def main():
    pp=Path('metadata/conditional_indel_fit_audit_v2_plan_20260927.json')
    plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():
        assert sha(p)==h,p
    producer=json.loads(Path(plan['producer_plan']).read_text())
    identity=json.loads(Path(plan['producer_launch']).read_text())
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    done={}
    terminal_seen=False
    with threadpool_limits(limits=1):
        analytic_check()
        while True:
            for job in producer['jobs']:
                if job['job_id'] in done:
                    continue
                path=Path(job['output'])/'receipt.json'
                if not path.exists():
                    continue
                # Production creates receipts after each job; skip any incomplete write.
                try:
                    json.loads(path.read_text())
                except json.JSONDecodeError:
                    continue
                result=audit(job,path)
                dest=out/(job['job_id']+'.json')
                dest.write_text(json.dumps(result,indent=2)+'\n')
                done[job['job_id']]=result
                print(job['job_id'],result['status'],flush=True)
            state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',identity['unit'],'--property=ActiveState,MainPID,Result,ExecMainStatus'],text=True).splitlines())
            if state['ActiveState']=='inactive':
                assert state['Result']=='success' and state['ExecMainStatus']=='0',state
                if len(done)<len(producer['jobs']) and not terminal_seen:
                    terminal_seen=True
                    continue
                assert len(done)==len(producer['jobs'])
                break
            assert state['ActiveState']=='active',state
            p=psutil.Process(identity['pid'])
            assert p.create_time()==identity['created'] and p.cmdline()==identity['cmdline']
            assert int(state['MainPID'])==p.pid
            print('waiting_for_verified_live_producer',len(done),flush=True)
            time.sleep(30)
    production_path=Path(producer['output'])/'receipt.json'
    production=json.loads(production_path.read_text())
    assert production['plan_sha256']==sha(plan['producer_plan'])
    for r in production['jobs']:
        assert r['receipt_sha256']==done[r['job_id']]['source_receipt_sha256']
    nonempty=[r for r in done.values() if 'starts' in r]
    failures=sum(not s['within_tolerance'] for r in nonempty for s in r['starts'])
    result=dict(status=('all_312_input_dispositions_and_918_fit_likelihoods_verified' if failures==0 else 'complete_replay_with_unresolved_likelihood_discrepancies'),
        starts_outside_tolerance=failures,
        high_precision_starts=sum(s['precision']=='70_digit_recursive_pruning' for r in nonempty for s in r['starts']),
        input_dispositions=len(done),nonempty_models=len(nonempty),
        starts_checked=sum(len(r['starts']) for r in nonempty),
        maximum_likelihood_error=max(s['absolute_error'] for r in nonempty for s in r['starts']),
        best_optimizer_failures=sum(not r['best_optimizer_success'] for r in nonempty),
        best_bound_models=sum(any(r['best_near_bound']) for r in nonempty),
        models_with_start_spread_above_001=sum(r['start_likelihood_spread']>.001 for r in nonempty),
        producer_terminal_state=state,producer_receipt_sha256=sha(production_path),
        audit_plan_sha256=sha(pp),artifacts={p.name:sha(p) for p in out.glob('*.json')},
        scope='Independent arithmetic verification at all fitted parameters. Bounds, optimizer flags and start sensitivity retained; not proof of convergence, model adequacy or final ancestral sequence validity.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
