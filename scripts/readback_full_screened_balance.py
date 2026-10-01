#!/usr/bin/env python3
"""Independently reconstruct all full-design balance, representation and reuse rows."""
import argparse
import csv
import gzip
import itertools
import json
import math
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from readback_background_control_balance import FEATURES, values, compare
from full_screened_balance_sources import load, KEY, MASKS, REUSE_FIELDS, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def statistics(x, y, baseline):
    finite = np.isfinite(x) & np.isfinite(y)
    a, b, c = [pd.Series(v, dtype=float) for v in [x[finite], y[finite], baseline[np.isfinite(baseline)]]]
    n, nb = len(a), len(c)
    mean = lambda v: float(v.iloc[0]) if v.nunique() == 1 else float(v.mean())
    std = lambda v: 0.0 if v.nunique() == 1 else float(v.std(ddof=1))
    result = dict(pairs=n,baseline_targets=nb,target_mean='',control_mean='',target_sd='',control_sd='',mean_difference='',standardized_mean_difference='',smd_status='no_pairs',mean_absolute_difference='',p95_absolute_difference='',maximum_absolute_difference='',baseline_target_mean=mean(c) if nb else '',selection_mean_shift='',selection_shift_in_baseline_sd='',selection_shift_status='no_pairs')
    if not n: return result
    ma, mb = mean(a), mean(b); absolute = (a-b).abs()
    result.update(target_mean=ma,control_mean=mb,mean_difference=ma-mb,mean_absolute_difference=float(absolute.mean()),p95_absolute_difference=float(absolute.quantile(.95)),maximum_absolute_difference=float(absolute.max()))
    if n == 1: result['smd_status'] = 'insufficient_pairs'
    else:
        sa, sb = std(a), std(b); pooled = math.sqrt((sa*sa+sb*sb)/2); result.update(target_sd=sa,control_sd=sb)
        if pooled == 0: result['smd_status'] = 'zero_pooled_variance'
        else: result.update(standardized_mean_difference=(ma-mb)/pooled,smd_status='estimable')
    if not nb: result['selection_shift_status'] = 'no_baseline'
    else:
        shift = ma-mean(c); result['selection_mean_shift'] = shift
        if nb == 1: result['selection_shift_status'] = 'insufficient_baseline'
        elif std(c) == 0: result['selection_shift_status'] = 'zero_baseline_variance'
        else: result.update(selection_shift_in_baseline_sd=shift/std(c),selection_shift_status='estimable')
    return result


def keyed(path, fields):
    result = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key = tuple(row[f] for f in fields); assert key not in result; result[key] = row
    return result


def run(plan_path, output):
    plan_path=Path(plan_path); plan=json.loads(plan_path.read_text()); source,bindings=load(plan,plan_path)
    out=Path(plan['output']);rp=out/'receipt.json';r=json.loads(rp.read_text());assert r['status']=='complete_full_screened_balance_pending_independent_readback'
    assert r['plan_sha256']==sha(plan_path) and r['matching_receipt_sha256']==sha(source['prior_receipt']) and r['source_hashes']==bindings and r['scientific_eligibility'] is False
    bind(bindings,rp)
    for name,digest in r['artifacts'].items():bind(bindings,out/name,digest)
    verify(bindings);assert plan['features']==FEATURES
    nodes=source['nodes'];vectors={k:pd.DataFrame.from_dict({ident:values(n) for ident,n in index.items()},orient='index',columns=FEATURES) for k,index in nodes.items()}
    flags={};seen=set()
    with gzip.open(source['root']/'target_policy_coverage_status.tsv.gz','rt') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=row['target_id'],row['policy'];assert key not in seen and row['policy'] in plan['policies'];seen.add(key)
            v=tuple(int(row['target_'+m+'_pass_bits']) for m in MASKS);n=nodes['target'][key[0]]
            assert all(0<=b<64 for b in v) and v[2]==v[0]&v[1] and row['guide']==n['guide'] and row['target_pair_key']==n['pair_key']
            assert key[0] not in flags or flags[key[0]]==v;flags[key[0]]=v
            if n['same_model']:assert v==(0,0,0)
    assert seen=={(tid,p) for tid in nodes['target'] for p in plan['policies']}
    baseline={g:vectors['target'].loc[[tid for tid,n in nodes['target'].items() if n['guide']==g]].to_numpy() for g in plan['guides']}
    eligible_baseline={(g,m,s['id']):vectors['target'].loc[[tid for tid,n in nodes['target'].items() if n['guide']==g and flags[tid][mi]//2**bit%2]].to_numpy() for g in plan['guides'] for mi,m in enumerate(MASKS) for bit,s in enumerate(plan['screens'])}
    cols=['target_id','background_id','guide','policy','scenario_id']+[side+'_'+m+'_pass_bits' for side in ['target','control'] for m in MASKS]
    selected=pd.read_csv(source['root']/'selected_pair_coverage.tsv.gz',sep='\t',usecols=cols,dtype={c:str for c in cols[:5]})
    assert len(selected)==plan['expected']['selected_records'] and not selected.duplicated(['target_id','policy','scenario_id']).any()
    for side,kind in [('target','target'),('control','background')]:
        ident='target_id' if kind=='target' else 'background_id'
        assert selected[ident].map({k:n['guide'] for k,n in nodes[kind].items()}).equals(selected.guide)
        assert ((selected[side+'_both_pass_bits']==(selected[side+'_full_pass_bits']&selected[side+'_plddt70_pass_bits'])) & selected[[side+'_'+m+'_pass_bits' for m in MASKS]].ge(0).all(axis=1) & selected[[side+'_'+m+'_pass_bits' for m in MASKS]].lt(64).all(axis=1)).all()
    for mi,m in enumerate(MASKS):assert selected.target_id.map({tid:v[mi] for tid,v in flags.items()}).equals(selected['target_'+m+'_pass_bits'])
    grouped=selected.groupby(KEY[:3]).indices;coverage=keyed(out/'coverage.tsv',KEY);balance=keyed(out/'balance.tsv',KEY+['feature'])
    expected=set(source['attrition']);assert set(coverage)==expected and set(balance)=={(*key,f) for key in expected for f in FEATURES}
    reuse_rows=retained_cells=strata=0
    with gzip.open(out/'retained_control_reuse.tsv.gz','rt') as handle:
        reuse_reader=csv.DictReader(handle,delimiter='\t');assert reuse_reader.fieldnames==REUSE_FIELDS
        for group in itertools.product(plan['guides'],plan['policies'],[s['scenario_id'] for s in source['scenarios']]):
            chunk=selected.iloc[grouped.get(group,[])];x=vectors['target'].loc[chunk.target_id].to_numpy();y=vectors['background'].loc[chunk.background_id].to_numpy();strata+=1
            for mi,mask in enumerate(MASKS):
                for bit,spec in enumerate(plan['screens']):
                    key=(*group,mask,spec['id']);attr=source['attrition'][key]
                    keep=(chunk['target_'+mask+'_pass_bits'].to_numpy()//2**bit%2==1)&(chunk['control_'+mask+'_pass_bits'].to_numpy()//2**bit%2==1)
                    retained=chunk.loc[keep];xx,yy=x[keep],y[keep];n=len(retained);retained_cells+=n;total=len(baseline[group[0]]);eb=eligible_baseline[group[0],mask,spec['id']]
                    assert n==int(attr['joint_pass_matched_records']) and len(chunk)==int(attr['matched_records']) and total==int(attr['all_target_records']) and len(eb)==int(attr['target_pass_all_records'])
                    tn=[nodes['target'][tid] for tid in retained.target_id];bn=[nodes['background'][bid] for bid in retained.background_id];assert not any(v['same_model'] for v in tn+bn)
                    counts=Counter(retained.background_id);pairs=Counter(v['pair_key'] for v in bn)
                    # Each group of k records has k*(1/k)=1 mass. Compute concentration from group counts, independently of per-record weight arrays.
                    nc,pc=len(counts),len(pairs);row=dict(zip(KEY,key));row.update(all_target_records=total,metadata_matched_records=len(chunk),metadata_unmatched_records=int(attr['unmatched_records']),target_eligible_all_records=len(eb),retained_matches=n,screen_excluded_matches=len(chunk)-n,retained_target_taxa=len({v['taxon_id'] for v in tn}),retained_target_families=len({v['family'] for v in tn}),distinct_target_genes=len({v['gene_'+s] for v in tn for s in ['a','b']}),distinct_target_physical_pairs=len({v['pair_key'] for v in tn}),distinct_control_nodes=nc,distinct_control_physical_pairs=pc,distinct_control_genes=len({v['gene_'+s] for v in bn for s in ['a','b']}),distinct_control_taxa=len({v['taxon_'+s] for v in bn for s in ['a','b']}),maximum_control_node_reuse=max(counts.values(),default=0),maximum_control_physical_pair_reuse=max(pairs.values(),default=0),top_five_control_node_fraction=sum(sorted(counts.values(),reverse=True)[:5])/n if n else '',top_five_control_physical_pair_fraction=sum(sorted(pairs.values(),reverse=True)[:5])/n if n else '',reciprocal_node_weight_sum=float(nc),reciprocal_physical_pair_weight_sum=float(pc),node_weight_kish_concentration=nc*nc/sum(1/k for k in counts.values()) if n else '',physical_pair_weight_kish_concentration=pc*pc/sum(1/k for k in pairs.values()) if n else '',zero_sequence_distance_targets=int(sum(xx[:,0]==0)),zero_sequence_distance_controls=int(sum(yy[:,0]==0)))
                    compare(coverage[key],row)
                    for bid in sorted(counts):
                        pair=nodes['background'][bid]['pair_key'];wanted=dict(zip(KEY,key),background_id=bid,background_pair_key=pair,retained_target_records=counts[bid],reciprocal_node_reuse_weight=1/counts[bid],retained_records_using_physical_pair=pairs[pair],reciprocal_physical_pair_reuse_weight=1/pairs[pair]);actual=next(reuse_reader,None);assert actual is not None;compare(actual,wanted);reuse_rows+=1
                    for col,feature in enumerate(FEATURES):
                        stats=statistics(xx[:,col],yy[:,col],baseline[group[0]][:,col]);extra={}
                        mapping={'baseline_targets':'baseline_targets','baseline_target_mean':'target_mean','selection_mean_shift':'mean_shift','selection_shift_in_baseline_sd':'shift_in_baseline_sd','selection_shift_status':'shift_status'}
                        for name,base in [('metadata_matched',x),('target_eligible_all',eb)]:
                            ss=statistics(xx[:,col],yy[:,col],base[:,col]);extra.update({name+'_'+new:ss[old] for old,new in mapping.items()})
                        wanted=dict(zip(KEY,key),feature=feature,retained_feature_excluded_pairs=n-stats['pairs'],**stats,**extra);compare(balance[(*key,feature)],wanted)
            if strata%12==0:print('Independent complete screened balance strata',strata,'/',len(plan['guides'])*len(plan['policies'])*plan['expected']['scenarios'],flush=True)
        assert next(reuse_reader,None) is None
    summary=dict(target_nodes=len(nodes['target']),background_nodes=len(nodes['background']),selected_records=len(selected),scenarios=plan['expected']['scenarios'],strata=strata,coverage_rows=len(coverage),balance_rows=len(balance),reuse_rows=reuse_rows,retained_selection_screen_cells=retained_cells,features=FEATURES,screens=plan['screens'],masks=MASKS,guides=plan['guides'],policies=plan['policies'])
    assert all(r[k]==summary[k] for k in SUMMARY_FIELDS);verify(bindings)
    result=dict(status='passed_full_screened_balance_independent_readback',plan_sha256=sha(plan_path),producer_receipt_sha256=sha(rp),**summary,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2));return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();run(args.plan,args.output)
