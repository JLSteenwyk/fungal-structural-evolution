#!/usr/bin/env python3
"""Full source-bound refinement readback and five-start likelihood comparison."""
import argparse
from collections import Counter,defaultdict
import json
from pathlib import Path
import re
import subprocess
import time
import numpy as np
import pandas as pd
import psutil
from Bio import Phylo
from ancestral_chain_attempt import sha,write_json


def branches(path):
    tree=Phylo.read(path,'newick')
    return {tuple(sorted(t.name for t in n.get_terminals())):float(n.branch_length or 0.) for n in tree.find_clades()}


def summarize_group(job_id,rows,baseline):
    assert len(rows)==5 and {r['start'] for r in rows}=={'precision_parameters','refreshed_parameters','low_ratio','balanced','high_ratio'}
    statuses=Counter(r['status'] for r in rows)
    if rows[0]['job']['character_count']==0:
        assert statuses=={'no_coded_characters':5}
        return dict(job_id=job_id,status='no_coded_characters',successful_starts=0)
    valid=[r for r in rows if r['status']=='refined_native_settings_and_replay_checked_requires_optimization_review']
    row=dict(job_id=job_id,status='incomplete_start_set_requires_review',successful_starts=len(valid))
    if not valid:return row
    best=max(valid,key=lambda r:r['native_rate_replay_log_likelihood']);likelihoods=[r['native_rate_replay_log_likelihood'] for r in valid]
    lower=max(baseline)
    row.update(best_start=best['start'],best_log_likelihood=max(likelihoods),previous_best_log_likelihood=lower,
        improvement_over_previous_best=max(likelihoods)-lower,start_log_likelihood_range=max(likelihoods)-min(likelihoods),
        starts_with_model_iteration_limit=sum(r['optimizer_settings']['model_iteration_limit_messages']>0 for r in valid),
        maximum_native_likelihood_discrepancy=max(abs(r['native_rate_likelihood_difference']) for r in valid),
        maximum_native_probability_discrepancy=max(r['native_rate_maximum_probability_difference'] for r in valid),
        best_alpha=best['fitted_parameters']['alpha'],best_gain=best['fitted_parameters']['gain'],best_loss=best['fitted_parameters']['loss'])
    if len(valid)==5:
        row['status']='all_starts_accounted_requires_scientific_review'
        row['all_start_likelihoods_within_1e_minus5']=bool(np.ptp(likelihoods)<=1e-5)
        row['best_not_worse_than_previous_within_1e_minus7']=bool(max(likelihoods)>=lower-1e-7)
    return row


def check_record(item,r):
    assert r['id']==item['id'] and r['start']==item['start'] and r['job']==item['job']
    if item['config'] is None:
        assert r['status']=='no_coded_characters';return
    rp=Path(r['attempt_receipt']);assert sha(rp)==r['attempt_receipt_sha256'];ar=json.loads(rp.read_text())
    assert json.loads((rp.parent.parent/'configuration.json').read_text())==item['config']
    for name,h in ar['artifacts'].items():assert sha(rp.parent/name)==h
    for name,h in item['config']['pins'].items():assert sha(name)==h
    if r['status']!='refined_native_settings_and_replay_checked_requires_optimization_review':return
    assert ar['exit_code']==0
    stdout=(rp.parent/'stdout.log').read_text()
    native=dict(re.findall(r'^(_\w+)\s+\((?:Float|Int|Str)\)\s+(\S+)',stdout,re.M))
    for key,value in item['expected_effective_options'].items():
        assert key in native
        if key=='_optimizationLevel':assert native[key]==value=='mid'
        else:assert float(native[key])==float(value),(key,native[key],value)
    assert float(native['_performOptimizationsBBL'])==0 and float(native['_performOptimizationsROOT'])==0
    parameter_file=Path(item['config']['command'][-1]);options=dict(l.split(None,1) for l in parameter_file.read_text().splitlines() if l.strip() and not l.startswith('#'))
    before=branches(options['_treeFile']);after=branches(rp.parent/'RESULTS/TheTree.INodes.ph')
    assert before.keys()==after.keys();np.testing.assert_allclose(list(before.values()),[after[k] for k in before],rtol=1e-10,atol=1e-12)
    text=(rp.parent/'RESULTS/EstimatedParameters.txt').read_text()
    reported=float(re.search(r'Log-likelihood=\s*([\deE.+-]+)',text).group(1));assert reported==r['reported_log_likelihood']
    pars=dict(re.findall(r'^(_\w+)\s+([\deE.+-]+)',text,re.M))
    for key,name in [('alpha','_userAlphaRate'),('gain','_userGain'),('loss','_userLoss')]:assert r['fitted_parameters'][key]==float(pars[name])
    assert r['native_rate_likelihood_difference']==r['native_rate_replay_log_likelihood']-reported
    assert np.isfinite([r['native_rate_replay_log_likelihood'],r['native_rate_likelihood_difference'],r['native_rate_maximum_probability_difference']]).all()
    assert not r['raw_values_clipped']
    assert r['probability_rows']==len(after)*item['job']['character_count']
    log=(rp.parent/'RESULTS/log.txt').read_text()
    assert r['optimizer_settings']['model_iteration_limit_messages']==log.count('Too many iterations in optimizeGainLossModel')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    config=json.loads(args.plan.read_text());ch=sha(args.plan)
    def verify():
        assert sha(args.plan)==ch
        for p,h in config['pins'].items():assert sha(p)==h,p
    verify();launch=json.loads(Path(config['producer_launch']).read_text())
    assert launch['plan_sha256']==sha(config['producer_plan'])
    while True:
        state=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','MainPID','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:break
        assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
        try:
            p=psutil.Process(launch['pid']);assert p.create_time()==launch['created'] and p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:time.sleep(1);continue
        time.sleep(30)
    assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
    verify();plan=json.loads(Path(config['producer_plan']).read_text());root=Path(plan['output']);prep=Path(plan['preparation'])
    for p,h in plan['pins'].items():assert sha(p)==h,p
    preparation=json.loads((prep/'receipt.json').read_text())
    for name,h in preparation['artifacts'].items():assert sha(prep/name)==h
    jobs={j['id']:j for j in json.loads((prep/'jobs.json').read_text())};assert len(jobs)==780
    receipt=json.loads((root/'receipt.json').read_text());assert receipt['plan_sha256']==sha(config['producer_plan']) and set(receipt['dispositions'])==set(jobs)
    baseline_root=Path(config['baseline']);baseline_receipt=json.loads((baseline_root/'receipt.json').read_text());baseline=defaultdict(list)
    for name,h in baseline_receipt['artifacts'].items():
        assert sha(baseline_root/name)==h
        r=json.loads((baseline_root/name).read_text());baseline[r['job_id']].append(r['native_rate_replay_log_likelihood'])
    assert len(baseline)==153 and all(len(v)==2 for v in baseline.values())
    groups=defaultdict(list);counts=Counter()
    for key,item in jobs.items():
        entry=receipt['dispositions'][key];path=Path(entry['path']);assert path.resolve().is_relative_to(root.resolve()) and sha(path)==entry['sha256']
        r=json.loads(path.read_text());assert r['status']==entry['status'];check_record(item,r)
        groups[item['job']['job_id']].append(r);counts[r['status']]+=1
    assert len(groups)==156
    rows=[summarize_group(k,v,baseline[k]) for k,v in sorted(groups.items())]
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    table=pd.DataFrame(rows);table.to_csv(out/'five_start_comparison.tsv',sep='\t',index=False)
    verify();write_json(out/'receipt.json',dict(status='full_refinement_output_readback_and_start_comparison_complete',plan_sha256=ch,
        source_receipt_sha256=sha(root/'receipt.json'),dispositions=len(jobs),input_groups=len(rows),status_counts=dict(counts),
        group_status_counts=dict(Counter(r['status'] for r in rows)),
        all_start_objectives_within_descriptive_tolerance=sum(r.get('all_start_likelihoods_within_1e_minus5',False) for r in rows),
        best_worse_than_previous=sum(not r['best_not_worse_than_previous_within_1e_minus7'] for r in rows if 'best_not_worse_than_previous_within_1e_minus7' in r),
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='All780 expected starts/empty cases, source attempts, effective controls, fixed trees and reported parameters checked; native-rate replay scalars retained, not recomputed here. Between-start agreement is descriptive, not a global optimum or posterior/model-adequacy certificate.'))


if __name__=='__main__':main()
