#!/usr/bin/env python3
"""Choose metadata-nearest controls across the complete prespecified sensitivity grid."""
import argparse,csv,gzip,itertools,json,math
from collections import Counter
from pathlib import Path
from background_match_covariates import BANDS
from run_ortholog_pair_guide_comparison import sha

SETS=['guide_native_ortholog','both_guides_native_ortholog','both_guides_unreported_parents']
SCENARIOS=[dict(scenario_id=f'S{i:02d}',background_set=s,distance_factor=f,band=b,focal_only=local) for i,(s,f,b,local) in enumerate(itertools.product(SETS,[1.25,1.5,2.],BANDS,[False,True]),1)]

def candidate(target,background,edge):
    d=target['sequence_distance'];db=background['sequence_distance']
    if d==0:
        if db!=0:raise ValueError('Zero distance requires exact zero background')
        seq=0.;ratio=1.
    else:
        if db<=0:raise ValueError('Positive distance requires positive background')
        seq=(math.log(db/d)/math.log(1.5))**2;ratio=max(db/d,d/db)
    ranks={}
    for band,(length,confidence,fraction) in BANDS.items():
        choices=[]
        for o in [0,1]:
            lr=float(edge['length_ratio_'+str(o)]);pc=float(edge['plddt_difference_'+str(o)]);lc=float(edge['lowconf_difference_'+str(o)])
            if lr<=length and pc<=confidence and lc<=fraction:
                score=seq+(math.log(lr)/math.log(1.25))**2+(pc/10)**2+(lc/.1)**2
                choices.append((score,background['node_id'],o))
        ranks[band]=min(choices) if choices else None
    return dict(ranks=ranks,background=background,focal=int(edge['focal_taxon']),distance_ratio=ratio,sequence_distance=db,shared=any(int(edge['shared_'+k]) for k in ['genes','models','sequences']))

def choose(candidates,scenario,target_distance):
    acceptable=[]
    for c in candidates:
        b=c['background'];s=scenario['background_set'];f=scenario['distance_factor'];rank=c['ranks'][scenario['band']]
        if c['shared'] or rank is None or (scenario['focal_only'] and not c['focal']):continue
        if s=='both_guides_native_ortholog' and not b['both_guides']:continue
        if s=='both_guides_unreported_parents' and not b['both_unreported_parents']:continue
        # Match the original multiplicative bounds exactly, avoiding reciprocal rounding.
        if not target_distance/f<=c['sequence_distance']<=target_distance*f:continue
        acceptable.append(rank)
    acceptable.sort()
    if not acceptable:return None
    score,bid,order=acceptable[0]
    return dict(background_id=bid,endpoint_order=order,score=score,eligible_candidates=len(acceptable),equal_score_candidates=sum(x[0]==score for x in acceptable),score_gap_to_second=acceptable[1][0]-score if len(acceptable)>1 else '')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();root=Path(plan['covariates']);rp=root/'receipt.json';r=json.loads(rp.read_text());proof=json.loads(Path(plan['readback']).read_text());assert proof['status']=='passed_full_background_match_covariate_readback' and proof['producer_receipt_sha256']==sha(rp)
    graph=Path(plan['graph']);gr=json.loads((graph/'receipt.json').read_text());assert r['graph_receipt_sha256']==sha(graph/'receipt.json')
    for name,h in gr['artifacts'].items():assert sha(graph/name)==h
    assert sha(root/'edge_covariates.tsv')==r['artifacts']['edge_covariates.tsv']
    def nodes(name):return {r['node_id']:r for line in open(graph/name) for r in [json.loads(line)]}
    targets=nodes('target_nodes.jsonl');backgrounds=nodes('background_nodes.jsonl');out=Path(plan['output']);out.mkdir(exist_ok=False)
    (out/'scenarios.json').write_text(json.dumps(SCENARIOS,indent=2)+'\n');summary=Counter();reuse=Counter();groups_seen=0;selected=0;edge_count=0
    with open(root/'edge_covariates.tsv') as ef,open(graph/'target_policy_dispositions.tsv') as df,gzip.open(out/'selections.tsv.gz','wt',compresslevel=3) as sf,(out/'target_policy_selection_status.tsv').open('w') as uf:
        groups=iter(itertools.groupby(csv.DictReader(ef,delimiter='\t'),key=lambda x:(x['target_id'],x['policy'])));current=next(groups,None)
        sw=csv.DictWriter(sf,fieldnames=['target_id','policy','scenario_id','background_id','endpoint_order','score','eligible_candidates','equal_score_candidates','score_gap_to_second'],delimiter='\t',lineterminator='\n');sw.writeheader()
        uw=csv.writer(uf,delimiter='\t',lineterminator='\n');uw.writerow(['target_id','policy','architecture_status','candidate_edges','shared_identity_edges','matched_scenarios','unmatched_scenarios'])
        for disposition in csv.DictReader(df,delimiter='\t'):
            tid=disposition['target_id'];policy=disposition['policy'];t=targets[tid];key=tid,policy;candidates=[]
            if current is not None and current[0]==key:
                for edge in current[1]:candidates.append(candidate(t,backgrounds[edge['background_id']],edge));edge_count+=1
                current=next(groups,None)
            assert len(candidates)==int(disposition['eligible_edges_within_factor_2'])
            matched=[];unmatched=[]
            for scenario in SCENARIOS:
                sid=scenario['scenario_id'];prefix=t['guide']+'|'+policy+'|'+sid;summary[prefix+'|targets']+=1;match=choose(candidates,scenario,t['sequence_distance'])
                if match is None:unmatched.append(sid);summary[prefix+'|unmatched']+=1
                else:
                    sw.writerow(dict(target_id=tid,policy=policy,scenario_id=sid,**match));matched.append(sid);selected+=1;summary[prefix+'|matched']+=1;reuse[t['guide'],policy,sid,match['background_id']]+=1
            uw.writerow([tid,policy,disposition['architecture_status'],len(candidates),sum(c['shared'] for c in candidates),','.join(matched),','.join(unmatched)]);groups_seen+=1
            if groups_seen%100000==0:print('Target-policy selections',groups_seen,flush=True)
        assert current is None
    assert edge_count==r['edges'] and groups_seen==gr['target_policy_dispositions']
    with (out/'background_reuse.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['guide','policy','scenario_id','background_id','selected_targets']);w.writerows((*k,n) for k,n in sorted(reuse.items()))
    verify();result=dict(status='complete_metadata_background_control_selection_pending_readback',plan_sha256=ph,source_covariate_receipt_sha256=sha(rp),target_policy_records=groups_seen,scenarios=len(SCENARIOS),selected_records=selected,source_edges=edge_count,counts=dict(summary),artifacts={p.name:sha(p) for p in out.iterdir()},scope='One deterministic nearest metadata control per eligible target/policy/scenario with replacement. Score sums squared log sequence-distance ratio/log1.5, endpoint-maximum log length ratio/log1.25, mean-pLDDT difference/10 and low-confidence-fraction difference/0.1. Endpoint mapping must meet every caliper jointly. Exact zero sequence distances match zero only. Ties use background node hash then order. Shared gene/model/sequence edges excluded explicitly; identical-model within-target/background pairs remain. All unmatched scenario IDs and control reuse retained. No structural response used. Independent selection readback, actual covariate balance, phylogenetic dependence, outcome eligibility and effects remain pending. Confidence adjustment can change the estimand; this is not a causal claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','target_policy_records','scenarios','selected_records']}),flush=True)
if __name__=='__main__':main()
