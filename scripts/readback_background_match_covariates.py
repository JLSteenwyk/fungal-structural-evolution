#!/usr/bin/env python3
"""Recompute all edge covariates using vector arrays and independently count overlaps."""
import argparse,json,time
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from run_ortholog_pair_guide_comparison import sha

def matrix_differences(target,background):
    # Arrays: edge, endpoint, [length, mean confidence, low-confidence fraction].
    assert target.shape==background.shape and target.shape[1:]==(2,3)
    assert np.isfinite(target).all() and np.isfinite(background).all()
    result={}
    for order in [0,1]:
        b=background if order==0 else background[:,::-1,:]
        result['length_ratio_'+str(order)]=np.maximum(target[:,:,0]/b[:,:,0],b[:,:,0]/target[:,:,0]).max(axis=1)
        result['plddt_difference_'+str(order)]=np.abs(target[:,:,1]-b[:,:,1]).max(axis=1)
        result['lowconf_difference_'+str(order)]=np.abs(target[:,:,2]-b[:,:,2]).max(axis=1)
    for name,limits in [('tight',(1.1,5,.05)),('moderate',(1.25,10,.1)),('wide',(1.5,15,.2))]:
        passes=np.zeros(len(target),dtype=bool)
        for order in [0,1]:
            passes|=np.logical_and.reduce([result[field+'_'+str(order)]<=limit for field,limit in zip(['length_ratio','plddt_difference','lowconf_difference'],limits)])
        result['within_'+name]=passes.astype(int)
    return result

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
    verify();sp=json.loads(Path(plan['source_plan']).read_text());root=Path(sp['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp);graph=Path(sp['source']);gr=json.loads((graph/'receipt.json').read_text())
    assert r['status']=='complete_background_match_covariates_pending_readback' and r['plan_sha256']==sha(plan['source_plan'])
    assert r['graph_receipt_sha256']==sha(graph/'receipt.json') and r['graph_readback_sha256']==sha(sp['readback'])
    for name,h in gr['artifacts'].items():assert sha(graph/name)==h
    table=root/'edge_covariates.tsv';assert sha(table)==r['artifacts'][table.name]
    def nodes(name):return pd.DataFrame([json.loads(line) for line in open(graph/name)]).set_index('node_id')
    targets=nodes('target_nodes.jsonl');backgrounds=nodes('background_nodes.jsonl');counts=Counter();n=0
    fields=['length','mean_ca_plddt','fraction_ca_plddt_below50'];columns=[field+'_'+side for side in ['a','b'] for field in fields]
    source=pd.read_csv(graph/'edges.tsv',sep='\t',chunksize=50000,keep_default_na=False)
    actual=pd.read_csv(table,sep='\t',chunksize=50000,keep_default_na=False,float_precision='round_trip')
    for edge in source:
        row=next(actual,None);assert row is not None
        pd.testing.assert_frame_equal(edge.reset_index(drop=True),row[edge.columns].reset_index(drop=True),check_dtype=False)
        t=targets.loc[edge.target_id];b=backgrounds.loc[edge.background_id];expected=matrix_differences(t[columns].to_numpy(float).reshape(-1,2,3),b[columns].to_numpy(float).reshape(-1,2,3))
        for name,values in expected.items():
            if name.startswith('within_'):assert np.array_equal(row[name].to_numpy(),values)
            else:assert np.allclose(row[name].to_numpy(),values,rtol=1e-13,atol=1e-13)
        for kind,cols in [('genes',['gene_a','gene_b']),('models',['model_id_a','version_a','model_id_b','version_b']),('sequences',['sequence_sha256_a','sequence_sha256_b'])]:
            def convert(v):return set(zip(v[::2],v[1::2])) if kind=='models' else set(v)
            overlaps=np.array([len(convert(x)&convert(y)) for x,y in zip(t[cols].itertuples(index=False,name=None),b[cols].itertuples(index=False,name=None))])
            assert np.array_equal(overlaps,row['shared_'+kind].to_numpy());expected['shared_'+kind]=(overlaps>0).astype(int)
        summaries=pd.DataFrame({'guide':t.guide.to_numpy(),'policy':edge.policy.to_numpy(),'edges':np.ones(len(edge),dtype=int),**{k:v for k,v in expected.items() if k.startswith(('within_','shared_'))}}).groupby(['guide','policy']).sum()
        for (guide,policy),summary in summaries.iterrows():
            for field,value in summary.items():counts[guide+'|'+policy+'|'+field]+=int(value)
        n+=len(edge);print('Verified covariate edges',n,flush=True)
    assert next(actual,None) is None and n==r['edges']==gr['edges'] and dict(counts)==r['counts']
    assert r['bands']=={'tight':[1.1,5.,.05],'moderate':[1.25,10.,.1],'wide':[1.5,15.,.2]}
    verify();assert sha(rp)==rh and sha(table)==r['artifacts'][table.name]
    result=dict(status='passed_full_background_match_covariate_readback',plan_sha256=ph,producer_receipt_sha256=rh,edges=n,scope='Every original graph field/order, both endpoint-mapping numeric differences, three joint caliper flags, shared identity counts and all summary counts independently reconstructed with array operations. No selected matching, covariate balance, causal interpretation or reliability guarantee.')
    Path(plan['output']).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
