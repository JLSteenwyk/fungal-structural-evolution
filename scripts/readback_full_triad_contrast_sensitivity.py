#!/usr/bin/env python3
"""Reconstruct every contrast envelope and screen from original raw fits with SQL."""
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
from full_triad_contrast_sensitivity_sources import load_sources, MASKS, CORES, QUALIFIED, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def classify(count, total, unique, low, high, epsilon):
    if count < total: return 'unavailable'
    assert count == total
    if unique != total: return 'nonunique_fit'
    if low > epsilon and high > epsilon: return 'positive'
    if low < -epsilon and high < -epsilon: return 'negative'
    if max(abs(low),abs(high)) <= epsilon: return 'within_numerical_tolerance'
    return 'sign_uncertain'


def run(path, output):
    plan=json.loads(Path(path).read_text());raw,groups,triads,bindings=load_sources(plan,path)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text());rh=sha(rp)
    assert receipt['status']=='complete_full_triad_contrast_sensitivity_pending_independent_readback'
    assert receipt['plan_sha256']==sha(path) and receipt['scientific_eligibility'] is False and receipt['scope']==plan['scope']
    assert receipt['contrast_numerical_tolerance_angstrom']==plan['contrast_numerical_tolerance_angstrom']
    assert json.loads((root/'configuration.json').read_text())==dict(plan_sha256=sha(path),source_hashes=bindings)
    screens=[s['id'] for s in plan['screens']];models={t['triad_id']:t['models'] for t in triads};assert len(models)==len(triads)
    db=sqlite3.connect(':memory:');cols=','.join('s'+str(i)+' INT' for i in range(len(screens)))
    db.execute('CREATE TABLE fits(t TEXT,m TEXT,c TEXT,o INT,v REAL,u INT,'+cols+',PRIMARY KEY(t,m,c,o))')
    fits=0;records=0;maximum=0.;counts=Counter();qualified=Counter()
    try:
        with gzip.open(raw,'rt') as f:
            reader=csv.DictReader(f,delimiter='\t');assert len(reader.fieldnames)==len(set(reader.fieldnames));batch=[]
            for r in reader:
                t=r['triad_id'];m=r['mask'];c=r['mapping_definition'];assert t in models and m in MASKS[:2] and c in CORES[:2]
                assert [[r['model_'+a],int(r['version_'+a])] for a in ['a','b','reference']]==models[t]
                order=[int(r['order_'+e]) for e in ['ab','ar','br']];assert all(v in [0,1] for v in order)
                number=order[0]*4+order[1]*2+order[2];v=float(r['rmsd_ar_minus_br']) if r['rmsd_ar_minus_br']!='' else None
                if v is not None:
                    assert math.isfinite(v)
                    assert abs(v-(float(r['rmsd_ar'])-float(r['rmsd_br'])))<=plan['contrast_numerical_tolerance_angstrom']
                u=int(r['fit_status']=='computed_unique_at_numeric_tolerance')
                flags=[int(r[s+'_pass']) for s in screens];assert all(b in [0,1] for b in flags)
                assert not any(flags) or u and v is not None
                batch.append((t,m,c,number,v,u,*flags));fits+=1
                if len(batch)==5000:db.executemany('INSERT INTO fits VALUES('+','.join('?' for _ in batch[0])+')',batch);batch.clear()
            if batch:db.executemany('INSERT INTO fits VALUES('+','.join('?' for _ in batch[0])+')',batch)
        assert fits==plan['expected']['source_fit_rows']==len(triads)*32
        assert db.execute('SELECT COUNT(*) FROM (SELECT t,m,c FROM fits GROUP BY t,m,c HAVING COUNT(*)!=8 OR MIN(o)!=0 OR MAX(o)!=7)').fetchone()[0]==0
        assert db.execute('SELECT COUNT(*) FROM (SELECT t FROM fits GROUP BY t HAVING COUNT(*)!=32)').fetchone()[0]==0
        number_chunks=(len(triads)+plan['checkpoint_triads']-1)//plan['checkpoint_triads']
        paths=[root/'checkpoints'/f'{i:06d}.jsonl.gz' for i in range(number_chunks)]
        assert set(paths)==set((root/'checkpoints').glob('*.jsonl.gz'))
        def checkpoints():
            for i,p in enumerate(paths):
                n=0
                with gzip.open(p,'rt') as f:
                    for line in f:n+=1;yield json.loads(line)
                assert n==9*min(plan['checkpoint_triads'],len(triads)-i*plan['checkpoint_triads'])
        def expected_records():
            for triad in triads:
                t=triad['triad_id'];leaf_bits={}
                for m,c in itertools.product(MASKS[:2],CORES[:2]):
                    bits=db.execute('SELECT '+','.join('SUM(s'+str(i)+'*(1<<o))' for i in range(len(screens)))+' FROM fits WHERE t=? AND m=? AND c=?',(t,m,c)).fetchone()
                    leaf_bits[m,c]=list(bits)
                for mask,core in itertools.product(MASKS,CORES):
                    ms=MASKS[:2] if mask=='both_masks' else [mask];cs=CORES[:2] if core=='both_cores' else [core]
                    selected=[(m,c) for m in ms for c in cs]
                    query='SELECT COUNT(*),COUNT(v),SUM(u),MIN(v),MAX(v) FROM fits WHERE t=? AND m IN ('+','.join('?' for _ in ms)+') AND c IN ('+','.join('?' for _ in cs)+')'
                    total,n,u,low,high=db.execute(query,(t,*ms,*cs)).fetchone();assert total==8*len(selected)
                    numerical=classify(n,total,u,low,high,plan['contrast_numerical_tolerance_angstrom']);strict=classify(n,total,u,low,high,0.)
                    screen_values={}
                    for i,sid in enumerate(screens):
                        bits=[leaf_bits[m,c][i] for m,c in selected];passed=sum(b==255 for b in bits)==len(selected)
                        if passed:assert n==u==total
                        screen_values[sid]=dict(leaf_order_pass_bits=bits,all_selected_orders_pass=passed,qualified_direction=numerical if passed else 'excluded_by_quality')
                    yield dict(triad_id=t,models=triad['models'],role_order=['a','b','reference'],mask_scenario=mask,core_scenario=core,
                        source_order_group_keys=[[t,m,c] for m,c in selected],expected_orders=total,available_orders=n,all_selected_orders_unique=u==total,
                        minimum_contrast_angstrom=low,maximum_contrast_angstrom=high,contrast_span_angstrom=high-low if n else None,
                        strict_zero_direction=strict,numerical_tolerance_direction=numerical,screens=screen_values)
        with gzip.open(root/'contrast_sensitivity.jsonl.gz','rt') as f:
            for wanted,exported,checkpoint in zip_longest(expected_records(),f,checkpoints()):
                assert wanted is not None and exported is not None and checkpoint is not None
                actual=json.loads(exported)
                canonical=lambda x:json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)
                assert canonical(actual)==canonical(wanted)==canonical(checkpoint),'False full contrast group'
                m=actual['mask_scenario'];c=actual['core_scenario'];counts[f'{m}|{c}|'+actual['numerical_tolerance_direction']]+=1
                for sid,s in actual['screens'].items():qualified[f'{m}|{c}|{sid}|'+s['qualified_direction']]+=1
                if actual['contrast_span_angstrom'] is not None:maximum=max(maximum,actual['contrast_span_angstrom'])
                records+=1
                if records%9000==0:print('Independent raw-fit contrast SQL groups',records,'/',plan['expected']['sensitivity_groups'],flush=True)
        fields=['mask_scenario','core_scenario','screen','qualified_direction','physical_triads','total_physical_triads'];wanted_rows=[]
        for m,c,s in itertools.product(MASKS,CORES,screens):
            assert sum(qualified[f'{m}|{c}|{s}|{d}'] for d in QUALIFIED)==len(triads)
            for d in QUALIFIED:wanted_rows.append(dict(zip(fields,map(str,[m,c,s,d,qualified[f'{m}|{c}|{s}|{d}'],len(triads)]))))
        with (root/'contrast_sensitivity_counts.tsv').open() as f:
            reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==fields and list(reader)==wanted_rows
    finally:db.close()
    summary=dict(measured_triads=len(triads),source_fit_rows=fits,source_order_groups=fits//8,sensitivity_groups=records,
        screen_decisions=records*len(screens),summary_rows=len(wanted_rows),direction_counts=dict(counts),qualified_direction_counts=dict(qualified),maximum_contrast_span=maximum)
    assert set(summary)==set(SUMMARY_FIELDS) and all(receipt[k]==v for k,v in summary.items())
    assert all(summary[k]==v for k,v in plan['expected'].items() if k in summary)
    bind(bindings,root/'configuration.json');assert receipt['source_hashes']==bindings
    expected_artifacts={'contrast_sensitivity.jsonl.gz','contrast_sensitivity_counts.tsv',*[str(p.relative_to(root)) for p in paths]}
    assert set(receipt['artifacts'])==expected_artifacts
    for name,digest in receipt['artifacts'].items():bind(bindings,root/name,digest)
    bind(bindings,rp,rh);verify(bindings)
    result=dict(status='passed_full_triad_contrast_sensitivity_raw_fit_sql_readback',**summary,plan_sha256=sha(path),producer_receipt_sha256=rh,
                source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','direction_counts','qualified_direction_counts']}),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.output)
