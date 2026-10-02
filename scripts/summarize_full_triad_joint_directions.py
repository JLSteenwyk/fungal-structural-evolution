#!/usr/bin/env python3
"""Retain joint contrast envelopes across both correspondence sources and all alternatives."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import itertools
import json
import math
from pathlib import Path
import shutil
from full_triad_joint_direction_sources import load_sources,selectors,MASKS,METHODS,METHOD_SCENARIOS,DEFINITIONS,CORE_SCENARIOS,PERMUTATIONS,DIRECTIONS,QUALIFIED,FULL_BITS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def classify(n,total,unique,low,high,epsilon):
    if n!=total:return 'unavailable'
    if not unique:return 'nonunique_fit'
    if low>epsilon:return 'positive'
    if high < -epsilon:return 'negative'
    if low>=-epsilon and high<=epsilon:return 'within_numerical_tolerance'
    return 'sign_uncertain'


def envelope(rows,orders,epsilon):
    ranges=[r['metric_ranges']['rmsd_ar_minus_br'] for r in rows];available=sum(r['available_orders'] for r in ranges);total=orders*len(rows)
    bounds=[(r['minimum'],r['maximum']) for r in ranges if r['available_orders']]
    low=min((a for a,b in bounds),default=None);high=max((b for a,b in bounds),default=None);unique=all(r['all_orders_unique'] for r in rows)
    return dict(expected_orders=total,available_orders=available,all_selected_orders_unique=unique,minimum_contrast_angstrom=low,
        maximum_contrast_angstrom=high,contrast_span_angstrom=high-low if bounds else None,
        strict_zero_direction=classify(available,total,unique,low,high,0.),numerical_tolerance_direction=classify(available,total,unique,low,high,epsilon))


def scenario(triad,sequence,structural,pairs,mask,method,core,plan):
    masks,methods,cores=selectors(mask,method,core);epsilon=plan['contrast_numerical_tolerance_angstrom']
    skeys=list(itertools.product(masks,methods));tkeys=list(itertools.product(masks,cores));pkeys=list(itertools.product(masks,methods,cores))
    a=envelope([sequence[k] for k in skeys],6,epsilon);b=envelope([structural[k] for k in tkeys],8,epsilon)
    n=a['available_orders']+b['available_orders'];total=a['expected_orders']+b['expected_orders'];unique=a['all_selected_orders_unique'] and b['all_selected_orders_unique']
    bounds=[r for r in [a,b] if r['available_orders']];low=min((r['minimum_contrast_angstrom'] for r in bounds),default=None);high=max((r['maximum_contrast_angstrom'] for r in bounds),default=None)
    joint=dict(expected_orders=total,available_orders=n,all_selected_orders_unique=unique,minimum_contrast_angstrom=low,maximum_contrast_angstrom=high,
        contrast_span_angstrom=high-low if bounds else None,strict_zero_direction=classify(n,total,unique,low,high,0.),numerical_tolerance_direction=classify(n,total,unique,low,high,epsilon))
    da=a['numerical_tolerance_direction'];db=b['numerical_tolerance_direction']
    relation='unavailable_or_nonunique' if da in DIRECTIONS[4:] or db in DIRECTIONS[4:] else 'matching_categories' if da==db else 'different_categories'
    screens={}
    for s in plan['screens']:
        sid=s['id'];bits=FULL_BITS
        for key in pkeys:bits &=pairs[key]['screens'][sid]['pair_pass_bits']
        passed=bits==FULL_BITS
        if passed:assert n==total and unique
        screens[sid]=dict(joint_pair_pass_bits=bits,all_selected_pairs_pass=passed,qualified_joint_direction=joint['numerical_tolerance_direction'] if passed else 'excluded_by_quality')
    return dict(triad_id=triad['triad_id'],models=triad['models'],role_order=['a','b','reference'],mask_scenario=mask,method_scenario=method,core_scenario=core,
        sequence_group_keys=[[triad['triad_id'],m,q] for m,q in skeys],structural_group_keys=[[triad['triad_id'],m,c] for m,c in tkeys],
        comparison_group_keys=[[triad['triad_id'],m,q,c] for m,q,c in pkeys],sequence_envelope=a,structural_envelope=b,joint_envelope=joint,
        source_direction_category_relation=relation,screens=screens)


def checked_group(row,triad,orders,plan):
    assert row['triad_id']==triad['triad_id'] and row['models']==triad['models'] and type(row['all_orders_unique']) is bool
    r=row['metric_ranges']['rmsd_ar_minus_br'];assert type(r['available_orders']) is int and 0<=r['available_orders']<=orders
    assert (r['available_orders']==0)==(r['minimum'] is None)==(r['maximum'] is None)
    if r['available_orders']:assert math.isfinite(r['minimum']) and math.isfinite(r['maximum']) and r['minimum']<=r['maximum']
    for s in plan['screens']:
        b=row['screens'][s['id']]['order_pass_bits'];assert type(b) is int and 0<=b<(1<<orders)


def execute(plan,path,paths,triads,bindings,out):
    config=dict(plan_sha256=sha(path),source_hashes=bindings);configuration=out/'configuration.json'
    if configuration.exists():assert json.loads(configuration.read_text())==config
    else:
        temp=out/'configuration.tmp';temp.write_text(json.dumps(config,indent=2)+'\n');temp.replace(configuration)
    chunks=out/'checkpoints';chunks.mkdir(exist_ok=True);expected_chunks=set();artifacts={};chunk=[];number=0;reused=0
    def save():
        nonlocal number,reused
        if not chunk:return
        p=chunks/f'{number:06d}.jsonl.gz';expected_chunks.add(p.name);data=gzip.compress(('\n'.join(chunk)+'\n').encode(),compresslevel=1,mtime=0)
        if p.exists():assert p.read_bytes()==data;reused+=1
        else:
            temp=p.with_suffix('.tmp');temp.write_bytes(data);temp.replace(p)
        artifacts[str(p.relative_to(out))]=sha(p);chunk.clear();number+=1
    counts=Counter();qualified=Counter();relations=Counter();maximum=0.;records=0
    with gzip.open(paths['sequence_groups'],'rt') as sf,gzip.open(paths['structural_groups'],'rt') as tf,gzip.open(paths['comparison_groups'],'rt') as pf,gzip.open(out/'joint_directions.jsonl.gz','wt',compresslevel=1) as dest:
        for i,triad in enumerate(triads,1):
            sequence={};structural={};pairs={}
            for m,q in itertools.product(MASKS[:2],METHODS):
                r=json.loads(next(sf));assert (r['mask'],r['sequence_method'])==(m,q) and r['order_layout']==PERMUTATIONS;checked_group(r,triad,6,plan);sequence[m,q]=r
            for m,c in itertools.product(MASKS[:2],DEFINITIONS):
                r=json.loads(next(tf));assert (r['mask'],r['mapping_definition'])==(m,c) and r['order_bits_layout']==[list(o) for o in itertools.product([0,1],repeat=3)];checked_group(r,triad,8,plan);structural[m,c]=r
            for m,q,c in itertools.product(MASKS[:2],METHODS,DEFINITIONS):
                r=json.loads(next(pf));assert r['triad_id']==triad['triad_id'] and r['models']==triad['models'] and (r['mask'],r['sequence_method'],r['mapping_definition'])==(m,q,c)
                assert r['pair_order_layout']==[[p,f'{i:03b}'] for p in PERMUTATIONS for i in range(8)] and r['pair_states']==48
                for s in plan['screens']:
                    b=r['screens'][s['id']]['pair_pass_bits'];assert type(b) is int and 0<=b<=FULL_BITS
                pairs[m,q,c]=r
            for m,q,c in itertools.product(MASKS,METHOD_SCENARIOS,CORE_SCENARIOS):
                r=scenario(triad,sequence,structural,pairs,m,q,c,plan);encoded=json.dumps(r,separators=(',',':'));dest.write(encoded+'\n');chunk.append(encoded);records+=1
                key=f'{m}|{q}|{c}';counts[key+'|'+r['joint_envelope']['numerical_tolerance_direction']]+=1;relations[key+'|'+r['source_direction_category_relation']]+=1
                if r['joint_envelope']['contrast_span_angstrom'] is not None:maximum=max(maximum,r['joint_envelope']['contrast_span_angstrom'])
                for sid,s in r['screens'].items():qualified[key+'|'+sid+'|'+s['qualified_joint_direction']]+=1
            if i%plan['checkpoint_triads']==0:
                save();temp=out/'state.tmp';temp.write_text(json.dumps(dict(status='running_full_joint_direction_sensitivity',completed_triads=i,total_triads=len(triads)))+'\n');temp.replace(out/'state.json')
                assert sum(p.stat().st_size for p in chunks.glob('*.jsonl.gz'))+dest.tell()<=plan['resources']['maximum_output_gib']*2**30 and shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
            if i%1000==0:print('Full joint correspondence direction',i,'/',len(triads),flush=True)
        assert next(sf,None) is next(tf,None) is next(pf,None) is None;save()
    assert {p.name for p in chunks.glob('*.jsonl.gz')}==expected_chunks
    fields=['mask_scenario','method_scenario','core_scenario','screen','qualified_direction','physical_triads','total_physical_triads'];rows=0
    with (out/'joint_direction_counts.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for m,q,c,s in itertools.product(MASKS,METHOD_SCENARIOS,CORE_SCENARIOS,[s['id'] for s in plan['screens']]):
            key=f'{m}|{q}|{c}|{s}';assert sum(qualified[key+'|'+d] for d in QUALIFIED)==len(triads)
            for d in QUALIFIED:w.writerow(dict(zip(fields,[m,q,c,s,d,qualified[key+'|'+d],len(triads)])));rows+=1
    summary=dict(measured_triads=len(triads),sequence_fit_rows=24*len(triads),structural_fit_rows=32*len(triads),source_sequence_groups=4*len(triads),source_structural_groups=4*len(triads),
        source_comparison_groups=8*len(triads),joint_direction_groups=records,screen_decisions=records*len(plan['screens']),summary_rows=rows,joint_direction_counts=dict(counts),
        qualified_direction_counts=dict(qualified),source_category_relation_counts=dict(relations),maximum_joint_contrast_span=maximum)
    assert all(summary[k]==v for k,v in plan['expected'].items());bind(bindings,configuration);verify(bindings)
    for name in ['joint_directions.jsonl.gz','joint_direction_counts.tsv']:artifacts[name]=sha(out/name)
    result=dict(status='complete_full_triad_joint_directions_pending_independent_readback',**summary,plan_sha256=sha(path),checked_reused_chunks=reused,
        contrast_numerical_tolerance_angstrom=plan['contrast_numerical_tolerance_angstrom'],source_hashes=bindings,artifacts=artifacts,scientific_eligibility=False,scope=plan['scope'])
    assert not (out/'receipt.json').exists();temp=out/'receipt.tmp';temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(out/'receipt.json')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts','joint_direction_counts','qualified_direction_counts','source_category_relation_counts']}),flush=True);return result


def run(path):
    plan=json.loads(Path(path).read_text());assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    paths,triads,bindings=load_sources(plan,path);out=Path(plan['output']);out.mkdir(exist_ok=True,parents=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (out/'receipt.json').exists();return execute(plan,path,paths,triads,bindings,out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
