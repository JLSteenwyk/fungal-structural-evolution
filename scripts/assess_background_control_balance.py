#!/usr/bin/env python3
"""Assess all selected-control strata, target selection shifts and background reuse."""
import argparse,csv,json,time
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from background_control_balance import FEATURES,feature_vector,balance
from run_ortholog_pair_guide_comparison import sha

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
    verify();source=Path(plan['selection']);rp=source/'receipt.json';r=json.loads(rp.read_text());proof=Path(plan['readback']);audit=json.loads(proof.read_text());assert audit['status']=='passed_full_background_control_selection_readback' and audit['producer_receipt_sha256']==sha(rp)
    for name,h in r['artifacts'].items():assert sha(source/name)==h
    graph=Path(plan['graph']);gr=json.loads((graph/'receipt.json').read_text())
    def nodes(name):
        assert sha(graph/name)==gr['artifacts'][name]
        return {n['node_id']:n for line in open(graph/name) for n in [json.loads(line)]}
    targets=nodes('target_nodes.jsonl');backgrounds=nodes('background_nodes.jsonl')
    tf=pd.DataFrame.from_dict({k:feature_vector(n) for k,n in targets.items()},orient='index',columns=FEATURES);bf=pd.DataFrame.from_dict({k:feature_vector(n) for k,n in backgrounds.items()},orient='index',columns=FEATURES)
    baseline={g:tf.loc[[k for k,n in targets.items() if n['guide']==g]].to_numpy() for g in ['profile','mafft']}
    pieces=defaultdict(list);taxa=defaultdict(set);families=defaultdict(set);same=Counter();n=0
    for chunk in pd.read_csv(source/'selections.tsv.gz',sep='\t',chunksize=50000):
        chunk['guide']=[targets[k]['guide'] for k in chunk.target_id];x=tf.loc[chunk.target_id].to_numpy();y=bf.loc[chunk.background_id].to_numpy()
        for key,indices in chunk.groupby(['guide','policy','scenario_id']).indices.items():
            pieces[key].append((x[indices],y[indices]));rows=chunk.iloc[indices]
            for tid,bid in zip(rows.target_id,rows.background_id):
                t=targets[tid];b=backgrounds[bid];taxa[key].add(t['taxon_id']);families[key].add(t['family']);same[key+('target',)]+=t['same_model'];same[key+('control',)]+=b['same_model']
        n+=len(chunk)
    assert n==r['selected_records'];reuse=defaultdict(list)
    for row in csv.DictReader(open(source/'background_reuse.tsv'),delimiter='\t'):reuse[row['guide'],row['policy'],row['scenario_id']].append(int(row['selected_targets']))
    scenarios=json.loads((source/'scenarios.json').read_text());out=Path(plan['output']);out.mkdir(exist_ok=False);balance_rows=0;strata=0
    with (out/'covariate_balance.tsv').open('w') as f,(out/'selection_coverage.tsv').open('w') as cf:
        writer=None;cw=None
        for guide in ['profile','mafft']:
            for policy in ['alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore']:
                for scenario in scenarios:
                    sid=scenario['scenario_id'];key=guide,policy,sid;parts=pieces.get(key,[])
                    x=np.concatenate([p[0] for p in parts]) if parts else np.empty((0,len(FEATURES)));y=np.concatenate([p[1] for p in parts]) if parts else x.copy();total=len(baseline[guide]);weights=reuse.get(key,[]);prefix='|'.join(key)
                    assert len(x)==r['counts'].get(prefix+'|matched',0) and total==r['counts'][prefix+'|targets'] and sum(weights)==len(x)
                    coverage=dict(guide=guide,policy=policy,scenario_id=sid,**{k:v for k,v in scenario.items() if k!='scenario_id'},all_target_records=total,matched_target_records=len(x),unmatched_target_records=total-len(x),matched_taxa=len(taxa[key]),matched_families=len(families[key]),unique_background_nodes=len(weights),maximum_background_reuse=max(weights,default=0),top_five_background_fraction=sum(sorted(weights,reverse=True)[:5])/len(x) if len(x) else '',identical_model_targets=same[key+('target',)],identical_model_controls=same[key+('control',)],zero_sequence_distance_targets=int(np.count_nonzero(x[:,0]==0)),zero_sequence_distance_controls=int(np.count_nonzero(y[:,0]==0)))
                    if cw is None:cw=csv.DictWriter(cf,list(coverage),delimiter='\t',lineterminator='\n');cw.writeheader()
                    cw.writerow(coverage);strata+=1
                    for j,feature in enumerate(FEATURES):
                        row=dict(guide=guide,policy=policy,scenario_id=sid,feature=feature,**balance(x[:,j],y[:,j],baseline[guide][:,j]))
                        if writer is None:writer=csv.DictWriter(f,list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
                        writer.writerow(row);balance_rows+=1
    verify()
    for name,h in r['artifacts'].items():assert sha(source/name)==h
    result=dict(status='complete_background_control_balance_pending_readback',plan_sha256=ph,selection_receipt_sha256=sha(rp),selection_readback_sha256=sha(proof),selected_records=n,strata=strata,balance_rows=balance_rows,features=FEATURES,artifacts={p.name:sha(p) for p in out.iterdir()},scope='Descriptive selected target/control means, sample SDs, pooled-SD standardized mean differences, paired absolute differences and target-selection shifts relative to all modeled target records in each guide. Eight permutation-invariant pair features. Positive-log-distance rows exclude zero-distance pairs without epsilon; raw-distance and zero counts retained. Zero/insufficient variance explicitly nonestimable. Reuse, taxa/families and identical-model pairs reported. Event-weighted summaries are not independent-sample tests, phylogenetic correction, reliability guarantees or biological effects; no automatic balance-pass threshold imposed.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','selected_records','strata','balance_rows']}),flush=True)
if __name__=='__main__':main()
