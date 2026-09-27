#!/usr/bin/env python3
"""Independently enumerate eligible choices for every target, policy and scenario."""
import argparse,csv,gzip,itertools,json,math,time
from collections import Counter
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha

BANDS={'tight':(1.1,5.,.05),'moderate':(1.25,10.,.1),'wide':(1.5,15.,.2)}
def scenarios():
    rows=[]
    for setting in ['guide_native_ortholog','both_guides_native_ortholog','both_guides_unreported_parents']:
        for factor in [1.25,1.5,2.]:
            for band in ['tight','moderate','wide']:
                for focal in [False,True]:rows.append(dict(scenario_id=f'S{len(rows)+1:02d}',background_set=setting,distance_factor=factor,band=band,focal_only=focal))
    return rows

def reconstruct(edges,backgrounds,d,s):
    candidates={};limit=BANDS[s['band']]
    for edge in edges:
        bid=edge['background_id'];b=backgrounds[bid];db=b['sequence_distance'];factor=s['distance_factor']
        if any(int(edge['shared_'+kind])!=0 for kind in ['genes','models','sequences']):continue
        if s['focal_only'] and int(edge['focal_taxon'])!=1:continue
        if s['background_set']=='both_guides_native_ortholog' and b['both_guides']!=1:continue
        if s['background_set']=='both_guides_unreported_parents' and b['both_unreported_parents']!=1:continue
        if db<d/factor or db>d*factor:continue
        seq=0. if d==0 and db==0 else (math.log(db/d)/math.log(1.5))**2
        for order in [0,1]:
            lr=float(edge['length_ratio_'+str(order)]);pc=float(edge['plddt_difference_'+str(order)]);lc=float(edge['lowconf_difference_'+str(order)])
            if lr>limit[0] or pc>limit[1] or lc>limit[2]:continue
            score=seq+(math.log(lr)/math.log(1.25))**2+(pc/10)**2+(lc/.1)**2
            old=candidates.get(bid)
            if old is None or (score,order)<old:candidates[bid]=(score,order)
    ranking=sorted((score,bid,order) for bid,(score,order) in candidates.items())
    if not ranking:return None
    score,bid,order=ranking[0]
    return dict(background_id=bid,endpoint_order=order,score=score,eligible_candidates=len(ranking),equal_score_candidates=sum(score==r[0] for r in ranking),score_gap_to_second=ranking[1][0]-score if len(ranking)>1 else '')

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
    verify();sp=json.loads(Path(plan['source_plan']).read_text());out=Path(sp['output']);rp=out/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    assert r['status']=='complete_metadata_background_control_selection_pending_readback' and r['plan_sha256']==sha(plan['source_plan'])
    cr=Path(sp['covariates'])/'receipt.json';assert r['source_covariate_receipt_sha256']==sha(cr)
    for name,h in r['artifacts'].items():assert sha(out/name)==h
    grid=scenarios();assert json.loads((out/'scenarios.json').read_text())==grid
    graph=Path(sp['graph']);gr=json.loads((graph/'receipt.json').read_text())
    def nodes(name):return {r['node_id']:r for line in open(graph/name) for r in [json.loads(line)]}
    targets=nodes('target_nodes.jsonl');backgrounds=nodes('background_nodes.jsonl');summary=Counter();reuse=Counter();total=selected=source_edges=0
    with open(Path(sp['covariates'])/'edge_covariates.tsv') as ef,open(graph/'target_policy_dispositions.tsv') as df,gzip.open(out/'selections.tsv.gz','rt') as sf,open(out/'target_policy_selection_status.tsv') as uf:
        groups=iter(itertools.groupby(csv.DictReader(ef,delimiter='\t'),key=lambda x:(x['target_id'],x['policy'])));group=next(groups,None);selections=iter(csv.DictReader(sf,delimiter='\t'));states=iter(csv.DictReader(uf,delimiter='\t'))
        for original in csv.DictReader(df,delimiter='\t'):
            tid=original['target_id'];policy=original['policy'];key=tid,policy;t=targets[tid];edges=[]
            if group is not None and group[0]==key:edges=list(group[1]);group=next(groups,None)
            assert len(edges)==int(original['eligible_edges_within_factor_2']);source_edges+=len(edges);matched=[];unmatched=[]
            for s in grid:
                prefix=t['guide']+'|'+policy+'|'+s['scenario_id'];summary[prefix+'|targets']+=1;expected=reconstruct(edges,backgrounds,t['sequence_distance'],s)
                if expected is None:unmatched.append(s['scenario_id']);summary[prefix+'|unmatched']+=1;continue
                actual=next(selections,None);assert actual is not None and (actual['target_id'],actual['policy'],actual['scenario_id'])==(tid,policy,s['scenario_id'])
                for field,value in expected.items():
                    if isinstance(value,float):assert math.isclose(float(actual[field]),value,rel_tol=1e-13,abs_tol=1e-13),(key,s,field)
                    else:assert actual[field]==str(value),(key,s,field)
                matched.append(s['scenario_id']);selected+=1;summary[prefix+'|matched']+=1;reuse[t['guide'],policy,s['scenario_id'],expected['background_id']]+=1
            state=next(states,None);expected=dict(target_id=tid,policy=policy,architecture_status=original['architecture_status'],candidate_edges=str(len(edges)),shared_identity_edges=str(sum(any(int(e['shared_'+k]) for k in ['genes','models','sequences']) for e in edges)),matched_scenarios=','.join(matched),unmatched_scenarios=','.join(unmatched));assert state==expected;total+=1
            if total%100000==0:print('Verified target-policy choices',total,flush=True)
        assert group is None and next(selections,None) is None and next(states,None) is None
    actual={}
    for row in csv.DictReader(open(out/'background_reuse.tsv'),delimiter='\t'):
        key=tuple(row[k] for k in ['guide','policy','scenario_id','background_id']);assert key not in actual;actual[key]=int(row['selected_targets'])
    assert actual==dict(reuse) and dict(summary)==r['counts'] and total==r['target_policy_records']==gr['target_policy_dispositions'] and selected==r['selected_records'] and source_edges==r['source_edges'] and r['scenarios']==len(grid)
    verify();assert sha(rp)==rh
    for name,h in r['artifacts'].items():assert sha(out/name)==h
    result=dict(status='passed_full_background_control_selection_readback',plan_sha256=ph,producer_receipt_sha256=rh,target_policy_records=total,scenarios=len(grid),scenario_decisions=total*len(grid),selected_records=selected,source_edges=source_edges,scope='Every scenario candidate independently enumerated; selected identity/order/score, eligible/tied count and runner-up gap reconstructed. All unmatched statuses, full target universe, guide/policy/scenario summaries and background reuse counts verified. Same frozen covariates and design, not independent biological evidence. Balance and outcome eligibility remain pending.')
    Path(plan['output']).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
