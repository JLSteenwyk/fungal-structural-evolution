#!/usr/bin/env python3
"""Reconstruct every balance/coverage cell directly from selected identities and node records."""
import argparse,csv,json,math,time
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from run_ortholog_pair_guide_comparison import sha

FEATURES=['sequence_distance','log_positive_sequence_distance','mean_log_length','log_length_asymmetry','mean_plddt','minimum_plddt','mean_lowconf_fraction','maximum_lowconf_fraction']
def values(n):
    d=n['sequence_distance'];a=n['length_a'];b=n['length_b'];p=n['mean_ca_plddt_a'];q=n['mean_ca_plddt_b'];u=n['fraction_ca_plddt_below50_a'];v=n['fraction_ca_plddt_below50_b']
    return [d,math.log(d) if d>0 else float('nan'),(math.log(a)+math.log(b))/2,abs(math.log(a/b)),(p+q)/2,min(p,q),(u+v)/2,max(u,v)]
def reconstruct(x,y,base):
    valid=np.isfinite(x)&np.isfinite(y);x=pd.Series(x[valid]);y=pd.Series(y[valid]);base=pd.Series(base[np.isfinite(base)]);n=len(x);nb=len(base)
    r=dict(pairs=n,baseline_targets=nb,target_mean='',control_mean='',target_sd='',control_sd='',mean_difference='',standardized_mean_difference='',smd_status='no_pairs',mean_absolute_difference='',p95_absolute_difference='',maximum_absolute_difference='',baseline_target_mean=float(base.mean()) if nb else '',selection_mean_shift='',selection_shift_in_baseline_sd='',selection_shift_status='no_pairs')
    if n==0:return r
    delta=x-y;absolute=delta.abs();r.update(target_mean=float(x.mean()),control_mean=float(y.mean()),mean_difference=float(x.mean()-y.mean()),mean_absolute_difference=float(absolute.mean()),p95_absolute_difference=float(absolute.quantile(.95)),maximum_absolute_difference=float(absolute.max()))
    if n==1:r['smd_status']='insufficient_pairs'
    else:
        sx=float(x.std(ddof=1));sy=float(y.std(ddof=1));scale=math.hypot(sx,sy)/math.sqrt(2);r.update(target_sd=sx,control_sd=sy)
        if scale==0:r['smd_status']='zero_pooled_variance'
        else:r.update(standardized_mean_difference=r['mean_difference']/scale,smd_status='estimable')
    if nb:
        shift=float(x.mean()-base.mean());r['selection_mean_shift']=shift
        if nb==1:r['selection_shift_status']='insufficient_baseline'
        else:
            sd=float(base.std(ddof=1))
            if sd==0:r['selection_shift_status']='zero_baseline_variance'
            else:r.update(selection_shift_in_baseline_sd=shift/sd,selection_shift_status='estimable')
    return r

def compare(row,expected):
    assert set(row)==set(expected)
    for field,v in expected.items():
        if isinstance(v,(float,np.floating)):assert math.isclose(float(row[field]),v,rel_tol=1e-9,abs_tol=1e-10),(field,row[field],v)
        else:assert row[field]==str(v),(field,row[field],v)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify()
    if 'producer' in plan:
        dep=plan['producer']
        while psutil.pid_exists(dep['pid']):
            try:
                proc=psutil.Process(dep['pid'])
                if abs(proc.create_time()-dep['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
                assert proc.cmdline()==dep['cmdline']
            except psutil.NoSuchProcess:break
            time.sleep(30)
    verify();sp=json.loads(Path(plan['source_plan']).read_text());root=Path(sp['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp);sel=Path(sp['selection']);sr=json.loads((sel/'receipt.json').read_text())
    assert r['status']=='complete_background_control_balance_pending_readback' and r['plan_sha256']==sha(plan['source_plan']) and r['selection_receipt_sha256']==sha(sel/'receipt.json') and r['selection_readback_sha256']==sha(sp['readback'])
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    for name,h in sr['artifacts'].items():assert sha(sel/name)==h
    graph=Path(sp['graph']);gr=json.loads((graph/'receipt.json').read_text())
    def nodes(name):
        assert sha(graph/name)==gr['artifacts'][name]
        return {n['node_id']:n for line in open(graph/name) for n in [json.loads(line)]}
    targets=nodes('target_nodes.jsonl');backgrounds=nodes('background_nodes.jsonl');tf=pd.DataFrame.from_dict({k:values(v) for k,v in targets.items()},orient='index');bf=pd.DataFrame.from_dict({k:values(v) for k,v in backgrounds.items()},orient='index')
    baseline={guide:tf.loc[[k for k,n in targets.items() if n['guide']==guide]].to_numpy() for guide in ['profile','mafft']}
    selected=pd.read_csv(sel/'selections.tsv.gz',sep='\t',usecols=['target_id','background_id','policy','scenario_id']);selected['guide']=selected.target_id.map({k:n['guide'] for k,n in targets.items()});groups=selected.groupby(['guide','policy','scenario_id']).indices
    scenarios=json.loads((sel/'scenarios.json').read_text());expectedkeys={(g,p,s['scenario_id']) for g in ['profile','mafft'] for p in ['alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore'] for s in scenarios};assert set(groups)<=expectedkeys
    coverage={}
    for row in csv.DictReader(open(root/'selection_coverage.tsv'),delimiter='\t'):
        key=tuple(row[k] for k in ['guide','policy','scenario_id']);assert key not in coverage;coverage[key]=row
    table={}
    for row in csv.DictReader(open(root/'covariate_balance.tsv'),delimiter='\t'):
        key=tuple(row[k] for k in ['guide','policy','scenario_id','feature']);assert key not in table;table[key]=row
    assert set(coverage)==expectedkeys and set(table)=={(*k,f) for k in expectedkeys for f in FEATURES}
    for key in sorted(expectedkeys):
        guide,policy,sid=key;s=next(x for x in scenarios if x['scenario_id']==sid);chunk=selected.iloc[groups.get(key,[])];tn=[targets[k] for k in chunk.target_id];bn=[backgrounds[k] for k in chunk.background_id];x=tf.loc[chunk.target_id].to_numpy();y=bf.loc[chunk.background_id].to_numpy();reuse=Counter(chunk.background_id);weights=sorted(reuse.values(),reverse=True);n=len(chunk);total=len(baseline[guide])
        c=dict(guide=guide,policy=policy,scenario_id=sid,**{k:v for k,v in s.items() if k!='scenario_id'},all_target_records=total,matched_target_records=n,unmatched_target_records=total-n,matched_taxa=len({v['taxon_id'] for v in tn}),matched_families=len({v['family'] for v in tn}),unique_background_nodes=len(reuse),maximum_background_reuse=max(weights,default=0),top_five_background_fraction=sum(weights[:5])/n if n else '',identical_model_targets=sum(v['same_model'] for v in tn),identical_model_controls=sum(v['same_model'] for v in bn),zero_sequence_distance_targets=sum(v['sequence_distance']==0 for v in tn),zero_sequence_distance_controls=sum(v['sequence_distance']==0 for v in bn))
        compare(coverage[key],c)
        for j,f in enumerate(FEATURES):compare(table[(*key,f)],dict(guide=guide,policy=policy,scenario_id=sid,feature=f,**reconstruct(x[:,j],y[:,j],baseline[guide][:,j])))
    assert len(selected)==r['selected_records']==sr['selected_records'] and len(coverage)==r['strata'] and len(table)==r['balance_rows'] and r['features']==FEATURES
    verify();assert sha(rp)==rh
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    result=dict(status='passed_full_background_control_balance_readback',plan_sha256=ph,producer_receipt_sha256=rh,selected_records=len(selected),strata=len(coverage),balance_rows=len(table),scope='All selected identities directly rejoined; feature definitions independently reconstructed, all moments/quantiles/standardized differences/nonestimable statuses and baseline shifts checked, every coverage/taxon/family/reuse/zero/identical-model count verified. No producer feature or balance helpers imported. Descriptive balance, not independent-sample tests or biological effects.')
    Path(plan['output']).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
