#!/usr/bin/env python3
"""Run full fixed-tree multistart refinement and retain numerical dispositions."""
import argparse
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
import csv
import fcntl
import json
from pathlib import Path
import re
import shutil
import subprocess
import traceback
import numpy as np
from Bio import Phylo,SeqIO
from ancestral_chain_attempt import run_attempt,sha,write_json
from readback_fastml_variant_roundoff import readback
from replay_fastml_indel_probabilities import infer,analytic_check


def effective_settings(stdout,log,item):
    values=dict(re.findall(r'^(_\w+)\s+\((?:Float|Int|Str)\)\s+(\S+)',stdout,re.M))
    expected={**item['expected_effective_options'],'_performOptimizationsBBL':'0','_performOptimizationsROOT':'0',
              '_isRootFreqEQstationary':'1','_performOptimizationsManyStarts':'0'}
    for name,value in expected.items():
        assert name in values,name
        if name=='_optimizationLevel':assert values[name]==value,(name,values[name],value)
        else:np.testing.assert_allclose(float(values[name]),float(value),rtol=1e-10,atol=0,err_msg=name)
    starts=re.findall(r'optimization starting- epsilonOptParam=([^ ]+) epsilonOptIter=\s*([^,]+), MaxNumIterations=(\d+)',log)
    assert starts and all(int(row[2])==100 for row in starts)
    initial=re.search(r' L=\s*\S+ gainLossRatio=\s*\S+ gain=\s*(\S+) loss=\s*(\S+) Alpha=\s*(\S+)',log)
    assert initial,'Missing initial model parameter diagnostic'
    observed=dict(zip(['gain','loss','alpha'],map(float,initial.groups())))
    for key,value in item['initial_parameters'].items():
        np.testing.assert_allclose(observed[key],value,rtol=6e-5,atol=1e-6,err_msg='Printed initial '+key)
    return dict(effective_options={k:values[k] for k in expected},initial_parameters_printed=observed,
                model_optimizer_calls=len(starts),model_iteration_limit_messages=log.count('Too many iterations in optimizeGainLossModel'))


def tree_branches(path):
    tree=Phylo.read(path,'newick');return {tuple(sorted(t.name for t in node.get_terminals())):float(node.branch_length or 0.) for node in tree.find_clades()}


def execute(item,plan):
    root=Path(plan['output'])/item['id'];root.mkdir(parents=True,exist_ok=True);final=root/'readback.json'
    if item['config'] is None:
        assert item['job']['character_count']==0
        write_json(final,dict(status='no_coded_characters',id=item['id'],job=item['job'],start=item['start']));return final
    rp=run_attempt(root/'attempts',item['config']);attempt=json.loads(rp.read_text())
    result=dict(status='native_refinement_failed',exit_code=attempt['exit_code'])
    if attempt['exit_code']==0:
        try:
            folder=rp.parent/'RESULTS'
            settings=effective_settings((rp.parent/'stdout.log').read_text(),(folder/'log.txt').read_text(),item)
            options=dict(line.split(None,1) for line in Path(item['config']['command'][-1]).read_text().splitlines() if line.strip() and not line.startswith('#'))
            before=tree_branches(options['_treeFile']);after=tree_branches(folder/'TheTree.INodes.ph')
            assert before.keys()==after.keys(),'Fixed tree topology changed'
            np.testing.assert_allclose([before[k] for k in before],[after[k] for k in before],rtol=1e-10,atol=1e-12)
            result=readback(rp.parent,item['job']);result['scipy_rate_replay']=dict(likelihood_difference=result['likelihood_difference'],maximum_probability_difference=result['maximum_probability_difference'])
            pars=result['fitted_parameters'];native=subprocess.check_output([plan['rate_exporter']],input='fit '+repr(pars['alpha'])+'\n',text=True)
            name,*values=native.strip().split('\t');assert name=='fit';rates=np.array(values,dtype=float)
            assert rates.shape==(4,) and np.isfinite(rates).all() and (rates>0).all()
            tree=Phylo.read(folder/'TheTree.INodes.ph','newick');tree.root.name=tree.root.name or tree.root.comment
            seq={r.id:str(r.seq)+'0' for r in SeqIO.parse(item['job']['characters'],'fasta')}
            post,ll=infer(tree,seq,pars['gain'],pars['loss'],rates);assert np.isfinite(ll).all() and ll[-1]<0
            corrected=float((ll[:-1]-np.log(-np.expm1(ll[-1]))).sum());maximum=0.;count=0
            with (folder/'AncestralReconstructPosterior.txt').open() as handle:
                for row in csv.DictReader(handle,delimiter='\t'):
                    value=abs(float(row['Prob'])-float(post[row['Node']][int(row['POS'])-1]));assert np.isfinite(value)
                    maximum=max(maximum,value);count+=1
            assert count==result['probability_rows']
            result.update(status='refined_native_settings_and_replay_checked_requires_optimization_review',optimizer_settings=settings,
                fixed_tree_branch_count=len(before),native_rates=rates.tolist(),native_rate_replay_log_likelihood=corrected,
                native_rate_likelihood_difference=corrected-result['reported_log_likelihood'],native_rate_maximum_probability_difference=maximum)
        except Exception as error:
            result=dict(status='refinement_validation_failed_requires_review',error=repr(error),traceback=traceback.format_exc())
    result.update(id=item['id'],job=item['job'],start=item['start'],initial_parameters=item['initial_parameters'],
        attempt_receipt=str(rp),attempt_receipt_sha256=sha(rp))
    write_json(final,result);return final


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();prep=Path(plan['preparation']);receipt=json.loads((prep/'receipt.json').read_text())
    for path,h in receipt['pins'].items():assert sha(path)==h,path
    for name,h in receipt['artifacts'].items():assert sha(prep/name)==h,name
    jobs=json.loads((prep/'jobs.json').read_text());assert len(jobs)==len({j['id'] for j in jobs})==780
    assert sum(j['config'] is not None for j in jobs)==765
    assert shutil.disk_usage('.').free>plan['minimum_free_disk_bytes']
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'controller.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    binding=out/'run_plan.json'
    if binding.exists():assert binding.read_bytes()==args.plan.read_bytes()
    else:binding.write_bytes(args.plan.read_bytes())
    analytic_check();completed={};pending={};iterator=iter(jobs);exhausted=False
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        while not exhausted or pending:
            while not exhausted and len(pending)<plan['workers']:
                try:item=next(iterator)
                except StopIteration:exhausted=True;break
                assert shutil.disk_usage('.').free>plan['minimum_free_disk_bytes']
                pending[pool.submit(execute,item,plan)]=item['id']
            if not pending:break
            done,_=wait(pending,return_when=FIRST_COMPLETED)
            for future in done:
                identifier=pending.pop(future);path=future.result();r=json.loads(path.read_text())
                completed[identifier]=dict(path=str(path),sha256=sha(path),status=r['status'])
                print('Refined FastML dispositions',len(completed),'/780',identifier,r['status'],flush=True)
                if r['status']=='refinement_validation_failed_requires_review':
                    raise RuntimeError('Refinement contract failed; stop further submissions and preserve outstanding outputs: '+identifier)
    verify();assert len(completed)==780
    write_json(out/'receipt.json',dict(status='full_multistart_refinement_dispositions_complete_requires_review',plan_sha256=ph,
        dispositions=completed,scope='All five starts and empty cases retained. Native settings, fixed tree and two rate discretization replays checked. Requires independent output audit, between-start comparison and stationarity/model adequacy before using ancestors.'))


if __name__=='__main__':main()
