#!/usr/bin/env python3
"""Reconstruct joint direction sensitivity from complete original raw-fit grids."""
import argparse
from collections import Counter
import csv
import gzip
import itertools
from itertools import zip_longest
import json
import math
from pathlib import Path
import sqlite3
from full_triad_joint_direction_sources import load_sources,MASKS,METHODS,METHOD_SCENARIOS,DEFINITIONS,CORE_SCENARIOS,PERMUTATIONS,DIRECTIONS,QUALIFIED,FULL_BITS,SUMMARY_FIELDS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def direction(n,total,u,low,high,epsilon):
    if n<total:return 'unavailable'
    assert n==total
    if u<total:return 'nonunique_fit'
    assert u==total
    if min(low,high)>epsilon:return 'positive'
    if max(low,high)<-epsilon:return 'negative'
    if max(abs(low),abs(high))<=epsilon:return 'within_numerical_tolerance'
    return 'sign_uncertain'


def run(path,output):
    plan=json.loads(Path(path).read_text());paths,triads,bindings=load_sources(plan,path);root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    assert r['status']=='complete_full_triad_joint_directions_pending_independent_readback' and r['plan_sha256']==sha(path)
    assert r['scientific_eligibility'] is False and r['scope']==plan['scope'] and r['contrast_numerical_tolerance_angstrom']==plan['contrast_numerical_tolerance_angstrom']
    assert json.loads((root/'configuration.json').read_text())==dict(plan_sha256=sha(path),source_hashes=bindings)
    screens=[s['id'] for s in plan['screens']];models={t['triad_id']:t['models'] for t in triads};assert len(models)==len(triads)
    db=sqlite3.connect(':memory:');db.execute('CREATE TABLE fits(t TEXT,k INT,m TEXT,x TEXT,o INT,v REAL,u INT,'+','.join('p'+str(i)+' INT' for i in range(len(screens)))+',PRIMARY KEY(t,k,m,x,o))')
    row_counts=[];counts=Counter();qualified=Counter();relations=Counter();records=0;maximum=0.
    try:
        for kind,name,axes in [(0,'raw_sequence',METHODS),(1,'raw_structural',DEFINITIONS)]:
            n=0;batch=[]
            with gzip.open(paths[name],'rt') as f:
                reader=csv.DictReader(f,delimiter='\t');assert len(reader.fieldnames)==len(set(reader.fieldnames))
                for row in reader:
                    tid=row['triad_id'];m=row['mask'];axis=row['method'] if kind==0 else row['mapping_definition'];assert tid in models and m in MASKS[:2] and axis in axes
                    assert [[row['model_'+a],int(row['version_'+a])] for a in ['a','b','reference']]==models[tid]
                    if kind==0:order=PERMUTATIONS.index(row['permutation'])
                    else:
                        digits=[int(row['order_'+e]) for e in ['ab','ar','br']];assert all(v in [0,1] for v in digits);order=4*digits[0]+2*digits[1]+digits[2]
                    v=float(row['rmsd_ar_minus_br']) if row['rmsd_ar_minus_br']!='' else None
                    if v is not None:assert math.isfinite(v) and abs(v-(float(row['rmsd_ar'])-float(row['rmsd_br'])))<=plan['contrast_numerical_tolerance_angstrom']
                    u=int(row['fit_status']=='computed_unique_at_numeric_tolerance');flags=[int(row[s+'_pass']) for s in screens];assert all(v in [0,1] for v in flags)
                    assert not any(flags) or u and v is not None
                    batch.append((tid,kind,m,axis,order,v,u,*flags));n+=1
                    if len(batch)==5000:db.executemany('INSERT INTO fits VALUES('+','.join('?' for _ in batch[0])+')',batch);batch.clear()
                if batch:db.executemany('INSERT INTO fits VALUES('+','.join('?' for _ in batch[0])+')',batch)
            row_counts.append(n);assert n==len(triads)*(24 if kind==0 else 32)
        assert db.execute('SELECT COUNT(*) FROM (SELECT t,k,m,x FROM fits GROUP BY t,k,m,x HAVING COUNT(*)!=CASE WHEN k=0 THEN 6 ELSE 8 END OR MIN(o)!=0 OR MAX(o)!=CASE WHEN k=0 THEN 5 ELSE 7 END)').fetchone()[0]==0
        assert db.execute('SELECT COUNT(*) FROM (SELECT t FROM fits GROUP BY t HAVING COUNT(*)!=56)').fetchone()[0]==0
        nchunks=(len(triads)+plan['checkpoint_triads']-1)//plan['checkpoint_triads'];chunks=[root/'checkpoints'/f'{i:06d}.jsonl.gz' for i in range(nchunks)]
        assert set(chunks)==set((root/'checkpoints').glob('*.jsonl.gz'))
        def archived():
            for i,p in enumerate(chunks):
                n=0
                with gzip.open(p,'rt') as f:
                    for line in f:n+=1;yield json.loads(line)
                assert n==27*min(plan['checkpoint_triads'],len(triads)-i*plan['checkpoint_triads'])
        epsilon=plan['contrast_numerical_tolerance_angstrom']
        def envelope(tid,kind,ms,xs):
            query='SELECT COUNT(*),COUNT(v),SUM(u),MIN(v),MAX(v) FROM fits WHERE t=? AND k=? AND m IN ('+','.join('?' for _ in ms)+') AND x IN ('+','.join('?' for _ in xs)+')'
            total,n,u,low,high=db.execute(query,(tid,kind,*ms,*xs)).fetchone();assert total==(6 if kind==0 else 8)*len(ms)*len(xs)
            return dict(expected_orders=total,available_orders=n,all_selected_orders_unique=u==total,minimum_contrast_angstrom=low,maximum_contrast_angstrom=high,
                contrast_span_angstrom=high-low if n else None,strict_zero_direction=direction(n,total,u,low,high,0.),numerical_tolerance_direction=direction(n,total,u,low,high,epsilon))
        def expected():
            for triad in triads:
                tid=triad['triad_id'];leaf_bits={}
                for kind,axes in [(0,METHODS),(1,DEFINITIONS)]:
                    for m,x in itertools.product(MASKS[:2],axes):
                        leaf_bits[kind,m,x]=db.execute('SELECT '+','.join('SUM(p'+str(i)+'*(1<<o))' for i in range(len(screens)))+' FROM fits WHERE t=? AND k=? AND m=? AND x=?',(tid,kind,m,x)).fetchone()
                seq_envelopes={};struct_envelopes={}
                for m,q in itertools.product(MASKS,METHOD_SCENARIOS):seq_envelopes[m,q]=envelope(tid,0,MASKS[:2] if m=='both_masks' else [m],METHODS if q=='both_methods' else [q])
                for m,c in itertools.product(MASKS,CORE_SCENARIOS):struct_envelopes[m,c]=envelope(tid,1,MASKS[:2] if m=='both_masks' else [m],DEFINITIONS if c=='both_cores' else [c])
                for m,q,c in itertools.product(MASKS,METHOD_SCENARIOS,CORE_SCENARIOS):
                    ms=MASKS[:2] if m=='both_masks' else [m];qs=METHODS if q=='both_methods' else [q];cs=DEFINITIONS if c=='both_cores' else [c]
                    a=seq_envelopes[m,q];b=struct_envelopes[m,c];total=a['expected_orders']+b['expected_orders'];n=a['available_orders']+b['available_orders']
                    unique=a['all_selected_orders_unique'] and b['all_selected_orders_unique'];u=total if unique else 0
                    query='SELECT MIN(v),MAX(v) FROM fits WHERE t=? AND m IN ('+','.join('?' for _ in ms)+') AND ((k=0 AND x IN ('+','.join('?' for _ in qs)+')) OR (k=1 AND x IN ('+','.join('?' for _ in cs)+')))'
                    low,high=db.execute(query,(tid,*ms,*qs,*cs)).fetchone()
                    joint=dict(expected_orders=total,available_orders=n,all_selected_orders_unique=unique,minimum_contrast_angstrom=low,maximum_contrast_angstrom=high,
                        contrast_span_angstrom=high-low if n else None,strict_zero_direction=direction(n,total,u,low,high,0.),numerical_tolerance_direction=direction(n,total,u,low,high,epsilon))
                    da=a['numerical_tolerance_direction'];dbdir=b['numerical_tolerance_direction']
                    relation='unavailable_or_nonunique' if da not in DIRECTIONS[:4] or dbdir not in DIRECTIONS[:4] else 'matching_categories' if da==dbdir else 'different_categories'
                    checks={}
                    for i,sid in enumerate(screens):
                        absent=0
                        for mask,method,core in itertools.product(ms,qs,cs):
                            sbits=leaf_bits[0,mask,method][i];tbits=leaf_bits[1,mask,core][i]
                            bits=sum(1<<(s*8+t) for s in range(6) for t in range(8) if sbits&(1<<s) and tbits&(1<<t))
                            absent |=FULL_BITS^bits
                        bits=FULL_BITS^absent;passed=bits==FULL_BITS
                        if passed:assert n==total and unique
                        checks[sid]=dict(joint_pair_pass_bits=bits,all_selected_pairs_pass=passed,qualified_joint_direction=joint['numerical_tolerance_direction'] if passed else 'excluded_by_quality')
                    yield dict(triad_id=tid,models=triad['models'],role_order=['a','b','reference'],mask_scenario=m,method_scenario=q,core_scenario=c,
                        sequence_group_keys=[[tid,mask,method] for mask,method in itertools.product(ms,qs)],structural_group_keys=[[tid,mask,core] for mask,core in itertools.product(ms,cs)],
                        comparison_group_keys=[[tid,mask,method,core] for mask,method,core in itertools.product(ms,qs,cs)],sequence_envelope=a,structural_envelope=b,joint_envelope=joint,
                        source_direction_category_relation=relation,screens=checks)
        canonical=lambda v:json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False)
        with gzip.open(root/'joint_directions.jsonl.gz','rt') as f:
            for wanted,line,checkpoint in zip_longest(expected(),f,archived()):
                assert wanted is not None and line is not None and checkpoint is not None;actual=json.loads(line)
                assert canonical(actual)==canonical(wanted)==canonical(checkpoint),'False joint contrast sensitivity group'
                key='|'.join(actual[k] for k in ['mask_scenario','method_scenario','core_scenario']);counts[key+'|'+actual['joint_envelope']['numerical_tolerance_direction']]+=1
                relations[key+'|'+actual['source_direction_category_relation']]+=1
                if actual['joint_envelope']['contrast_span_angstrom'] is not None:maximum=max(maximum,actual['joint_envelope']['contrast_span_angstrom'])
                for sid,s in actual['screens'].items():qualified[key+'|'+sid+'|'+s['qualified_joint_direction']]+=1
                records+=1
                if records%27000==0:print('Independent raw-fit joint direction SQL',records,'/',plan['expected']['joint_direction_groups'],flush=True)
        fields=['mask_scenario','method_scenario','core_scenario','screen','qualified_direction','physical_triads','total_physical_triads'];rows=[]
        for m,q,c,s in itertools.product(MASKS,METHOD_SCENARIOS,CORE_SCENARIOS,screens):
            key=f'{m}|{q}|{c}|{s}';assert sum(qualified[key+'|'+d] for d in QUALIFIED)==len(triads)
            for d in QUALIFIED:rows.append(dict(zip(fields,map(str,[m,q,c,s,d,qualified[key+'|'+d],len(triads)]))))
        with (root/'joint_direction_counts.tsv').open() as f:
            reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==fields and list(reader)==rows
    finally:db.close()
    summary=dict(measured_triads=len(triads),sequence_fit_rows=row_counts[0],structural_fit_rows=row_counts[1],source_sequence_groups=row_counts[0]//6,source_structural_groups=row_counts[1]//8,
        source_comparison_groups=len(triads)*8,joint_direction_groups=records,screen_decisions=records*len(screens),summary_rows=len(rows),joint_direction_counts=dict(counts),
        qualified_direction_counts=dict(qualified),source_category_relation_counts=dict(relations),maximum_joint_contrast_span=maximum)
    assert set(summary)==set(SUMMARY_FIELDS) and all(r[k]==v for k,v in summary.items()) and all(summary[k]==v for k,v in plan['expected'].items())
    bind(bindings,root/'configuration.json');assert r['source_hashes']==bindings
    assert set(r['artifacts'])=={'joint_directions.jsonl.gz','joint_direction_counts.tsv',*[str(p.relative_to(root)) for p in chunks]}
    for name,h in r['artifacts'].items():bind(bindings,root/name,h)
    bind(bindings,rp,rh);verify(bindings)
    result=dict(status='passed_full_triad_joint_directions_raw_fit_sql_readback',**summary,plan_sha256=sha(path),producer_receipt_sha256=rh,
        source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','joint_direction_counts','qualified_direction_counts','source_category_relation_counts']}),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.output)
