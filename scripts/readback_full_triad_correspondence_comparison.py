#!/usr/bin/env python3
"""Independently rebuild complete source/order pair comparisons with SQL."""
import argparse
from collections import Counter
import csv
import gzip
import json
import math
from pathlib import Path
import sqlite3

from full_triad_correspondence_comparison_sources import load_sources,blocks,NUMERIC_FIELDS,SEQUENCE_NUMERIC_FIELDS,METHODS,PERMUTATIONS,DEFINITIONS,pair_fields,SUMMARY_FIELDS
from readback_whole_protein_common_fits import compare_row
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def compare(actual,wanted,tolerance=1e-9):
    maximum=0.
    if isinstance(wanted,dict):
        assert isinstance(actual,dict) and set(actual)==set(wanted),'Changed summary keys'
        for k in wanted:maximum=max(maximum,compare(actual[k],wanted[k],tolerance))
    elif isinstance(wanted,list):
        assert isinstance(actual,list) and len(actual)==len(wanted)
        for a,w in zip(actual,wanted):maximum=max(maximum,compare(a,w,tolerance))
    elif isinstance(wanted,float):
        assert type(actual) in [float,int] and math.isfinite(actual)
        maximum=abs(actual-wanted);assert maximum<=tolerance,(actual,wanted)
    else:assert type(actual) is type(wanted) and actual==wanted,(actual,wanted)
    return maximum


def overlap(left,right):
    a=sorted(left);b=sorted(right);i=j=common=0
    while i<len(a) and j<len(b):
        if a[i]==b[j]:common+=1;i+=1;j+=1
        elif a[i]<b[j]:i+=1
        else:j+=1
    return common,len(a)+len(b)-common


def aggregate(db,columns,condition,params,availability,table='source'):
    sql=','.join(f'COUNT({c}),MIN({c}),MAX({c}),MAX({c})-MIN({c}),AVG({c})' for c in columns)
    values=db.execute('SELECT '+sql+' FROM '+table+' WHERE '+condition,params).fetchone()
    return {field:dict(zip([availability,'minimum','maximum','span','mean'],values[5*i:5*(i+1)])) for i,field in enumerate(columns)}


def run(plan_path,output):
    plan=json.loads(Path(plan_path).read_text());source,bindings=load_sources(plan,plan_path)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text());rh=sha(rp)
    assert receipt['status']=='complete_full_triad_correspondence_comparison_pending_independent_readback'
    assert receipt['plan_sha256']==sha(plan_path) and receipt['scientific_eligibility'] is False and receipt['scope']==plan['scope']
    assert json.loads((root/'configuration.json').read_text())==dict(plan_sha256=sha(plan_path),source_hashes=bindings)
    counts=Counter();seq_counts=Counter();joint_counts=Counter();passes=Counter();all_passes=Counter()
    maxima={f:0. for f in NUMERIC_FIELDS};groups=sgroups=states=seqrows=structrows=0;difference=0.;expected_artifacts=set()
    names=['sequence_order_robustness.jsonl.gz','correspondence_pair_states.tsv.gz','correspondence_comparison_groups.jsonl.gz']
    db=sqlite3.connect(':memory:')
    numeric=','.join(f'n{i} REAL' for i in range(len(SEQUENCE_NUMERIC_FIELDS)))
    screens=','.join(f's{i}_{name} '+('TEXT' if name=='ex' else 'INT') for i in range(len(plan['screens'])) for name in ['pass','core','inherited','ex'])
    db.execute('CREATE TABLE source(kind TEXT,mask TEXT,k TEXT,o INT,status TEXT,h TEXT,native TEXT,payload TEXT,'+numeric+','+screens+',PRIMARY KEY(kind,mask,k,o))')
    db.execute('CREATE TABLE pairs(o INT PRIMARY KEY,identical INT,unique_fit INT,agreement TEXT,'+
        ','.join(f'n{i} REAL' for i in range(len(NUMERIC_FIELDS)))+','+','.join(f's{i} INT' for i in range(len(plan['screens'])))+')')
    try:
        with gzip.open(root/names[0],'rt') as sf,gzip.open(root/names[1],'rt') as pf,gzip.open(root/names[2],'rt') as jf:
            pair_reader=csv.DictReader(pf,delimiter='\t');assert pair_reader.fieldnames==pair_fields(plan)
            for number,(triad,sequence,structural,st,tt) in enumerate(blocks(source,raw_msa=True),1):
                db.execute('DELETE FROM source');checkpoint_sequences=[];checkpoint_pairs=[];checkpoint_groups=[]
                for kind,records in [('sequence',sequence),('structural',structural)]:
                    for (mask,key,order),row in records.items():
                        index=PERMUTATIONS.index(order) if kind=='sequence' else order
                        values=[kind,mask,key,index,row['fit_status'],row['triples_sha256'],row.get('sequence_native_status'),json.dumps(row),
                            *[None if row.get(f,'')=='' else float(row[f]) for f in SEQUENCE_NUMERIC_FIELDS]]
                        assert all(v is None or math.isfinite(v) for v in values[8:])
                        for i,s in enumerate(plan['screens']):
                            sid=s['id'];flags=[int(row[sid+suffix]) for suffix in ['_pass','_core_pass','_three_pair_pass']]
                            assert all(v in [0,1] for v in flags) and flags[0]==flags[1]*flags[2]
                            values.extend([*flags,row[sid+'_exclusions']])
                        db.execute('INSERT INTO source VALUES('+','.join('?' for _ in values)+')',values)
                for mask in ['full','plddt70']:
                    for method in METHODS:
                        condition='kind=? AND mask=? AND k=?';params=('sequence',mask,method)
                        raw=list(db.execute('SELECT status,h,native,payload FROM source WHERE '+condition+' ORDER BY o',params));assert len(raw)==6
                        rows=[json.loads(r[3]) for r in raw];statuses=[r[0] for r in raw];unique=all(s=='computed_unique_at_numeric_tolerance' for s in statuses)
                        values=aggregate(db,[f'n{i}' for i in range(len(SEQUENCE_NUMERIC_FIELDS))],condition,params,'available_orders')
                        metrics={f:values['n'+str(i)] for i,f in enumerate(SEQUENCE_NUMERIC_FIELDS)};c=metrics['rmsd_ar_minus_br']
                        direction='uncomputed_or_excluded_order' if c['available_orders']<6 else 'nonunique_order_fit' if not unique else 'positive' if c['minimum']>0 else 'negative' if c['maximum']<0 else 'includes_zero'
                        wanted=dict(triad_id=triad['triad_id'],models=triad['models'],mask=mask,sequence_method=method,order_layout=PERMUTATIONS,
                            order_fit_statuses=statuses,order_native_statuses=[r[2] for r in raw],order_native_payload_hashes=[r['native_payload_sha256'] for r in rows],
                            order_triplet_hashes=[r[1] for r in raw],distinct_order_triplet_sets=db.execute('SELECT COUNT(DISTINCT h) FROM source WHERE '+condition,params).fetchone()[0],
                            all_orders_unique=unique,strict_all_orders_contrast_direction=direction,metric_ranges=metrics,screens={})
                        for i,s in enumerate(plan['screens']):
                            sid=s['id'];bits,core,inherited,n,all_pass,any_pass=db.execute(f'SELECT SUM(s{i}_pass << o),SUM(s{i}_core << o),SUM(s{i}_inherited << o),SUM(s{i}_pass),MIN(s{i}_pass),MAX(s{i}_pass) FROM source WHERE '+condition,params).fetchone()
                            exclusions=[r[0].split(';') if r[0] else [] for r in db.execute(f'SELECT s{i}_ex FROM source WHERE '+condition+' ORDER BY o',params)]
                            wanted['screens'][sid]=dict(order_pass_bits=bits,core_order_pass_bits=core,inherited_order_pass_bits=inherited,passing_orders=n,all_orders_pass=bool(all_pass),any_order_pass=bool(any_pass),order_exclusions=exclusions)
                            if all_pass:assert unique
                        actual=json.loads(next(sf));difference=max(difference,compare(actual,wanted));checkpoint_sequences.append(actual)
                        seq_counts[mask+':'+method+':'+direction]+=1;sgroups+=1
                        for definition in DEFINITIONS:
                            db.execute('DELETE FROM pairs');paired=[]
                            columns=','.join(f'a.n{i}-b.n{i}' for i in range(len(NUMERIC_FIELDS)))
                            sql='SELECT a.o,b.o,a.status,b.status,a.h,b.h,a.payload,b.payload,'+columns+' FROM source a JOIN source b ON a.mask=b.mask WHERE a.kind="sequence" AND b.kind="structural" AND a.mask=? AND a.k=? AND b.k=? ORDER BY a.o,b.o'
                            grid=list(db.execute(sql,(mask,method,definition)));assert len(grid)==48
                            for seq_index,index,ss,ts,sh,th,sp,tp,*numbers in grid:
                                permutation=PERMUTATIONS[seq_index];s,t=json.loads(sp),json.loads(tp)
                                a,b=st[mask,method,permutation],tt[mask,definition,index];common,union=overlap(a,b)
                                is_unique=ss==ts=='computed_unique_at_numeric_tolerance'
                                if is_unique:
                                    sa=float(s['rmsd_ar_minus_br']);sb=float(t['rmsd_ar_minus_br']);agreement='same' if (sa>0)-(sa<0)==(sb>0)-(sb<0) else 'different'
                                else:agreement='unavailable'
                                wanted=dict(triad_id=triad['triad_id'],mask=mask,sequence_method=method,sequence_permutation=permutation,mapping_definition=definition,
                                    structural_order=f'{index:03b}',sequence_fit_status=ss,structural_fit_status=ts,sequence_triples_sha256=sh,structural_triples_sha256=th,
                                    sequence_common_residues=len(a),structural_common_residues=len(b),intersection_triples=common,union_triples=union,triple_jaccard=common/union if union else '',
                                    identical_correspondence=int(common==len(a)==len(b)),both_unique_numeric=int(is_unique),strict_numeric_contrast_sign_agreement=agreement)
                                wanted.update({'sequence_minus_structural_'+field:value if value is not None else '' for field,value in zip(NUMERIC_FIELDS,numbers)})
                                flags=[];key=mask+':'+method+':'+definition
                                for screen in plan['screens']:
                                    sid=screen['id'];ap=int(s[sid+'_pass']);bp=int(t[sid+'_pass']);flags.append(ap*bp)
                                    wanted.update({sid+'_sequence_pass':ap,sid+'_structural_pass':bp,sid+'_both_pass':ap*bp});passes[key+':'+sid]+=ap*bp
                                actual=next(pair_reader,None);assert actual is not None,'Missing correspondence order pair'
                                difference=max(difference,compare_row(actual,wanted));paired.append(actual);checkpoint_pairs.append(actual);states+=1
                                counts[key+':'+ss+':'+ts]+=1
                                for field,value in zip(NUMERIC_FIELDS,numbers):
                                    if value is not None:maxima[field]=max(maxima[field],abs(value))
                                values=[seq_index*8+index,wanted['identical_correspondence'],int(is_unique),agreement,*numbers,*flags]
                                db.execute('INSERT INTO pairs VALUES('+','.join('?' for _ in values)+')',values)
                            values=aggregate(db,[f'n{i}' for i in range(len(NUMERIC_FIELDS))],'1',[],'available_pairs','pairs')
                            metrics={'sequence_minus_structural_'+field:values['n'+str(i)] for i,field in enumerate(NUMERIC_FIELDS)};c=metrics['sequence_minus_structural_rmsd_ar_minus_br']
                            identical,unique=db.execute('SELECT SUM(identical),SUM(unique_fit) FROM pairs').fetchone()
                            direction='uncomputed_or_nonunique_pair' if unique<48 else 'positive' if c['minimum']>0 else 'negative' if c['maximum']<0 else 'includes_zero'
                            wanted=dict(triad_id=triad['triad_id'],models=triad['models'],mask=mask,sequence_method=method,mapping_definition=definition,
                                sequence_robustness_key=[triad['triad_id'],mask,method],structural_robustness_key=[triad['triad_id'],mask,definition],
                                pair_order_layout=[[p,f'{i:03b}'] for p in PERMUTATIONS for i in range(8)],pair_states=48,identical_correspondence_pairs=identical,both_unique_numeric_pairs=unique,
                                strict_numeric_contrast_sign_agreement_counts=dict(db.execute('SELECT agreement,COUNT(*) FROM pairs GROUP BY agreement')),
                                strict_all_pairs_contrast_difference_direction=direction,metric_difference_ranges=metrics,screens={})
                            joint_counts[key+':'+direction]+=1
                            for i,s in enumerate(plan['screens']):
                                sid=s['id'];bits,n,all_pass,any_pass=db.execute(f'SELECT SUM(s{i} << o),SUM(s{i}),MIN(s{i}),MAX(s{i}) FROM pairs').fetchone()
                                wanted['screens'][sid]=dict(pair_pass_bits=bits,passing_pairs=n,all_pairs_pass=bool(all_pass),any_pair_pass=bool(any_pass));all_passes[key+':'+sid]+=all_pass
                            actual=json.loads(next(jf));difference=max(difference,compare(actual,wanted));checkpoint_groups.append(actual);groups+=1
                checkpoint=root/'checkpoints'/(triad['triad_id']+'.json.gz');payload=json.loads(gzip.decompress(checkpoint.read_bytes()))
                assert payload['triad_id']==triad['triad_id'] and set(payload)=={'triad_id','sequence_groups','pairs','comparison_groups'}
                compare(payload['sequence_groups'],checkpoint_sequences);compare(payload['comparison_groups'],checkpoint_groups)
                assert [{k:str(v) for k,v in r.items()} for r in payload['pairs']]==checkpoint_pairs
                expected_artifacts.add(str(checkpoint.relative_to(root)));seqrows+=len(sequence);structrows+=len(structural)
                if number%1000==0:print('Independent raw-MSA/SQL correspondence comparison',number,'/',len(source['ready']),'pairs',states,flush=True)
            assert next(sf,None) is next(jf,None) is next(pair_reader,None) is None,'Extra correspondence output row'
    finally:db.close()
    summary=dict(ordered_model_triads=source['ordered_model_triads'],source_ready_triads=len(source['ready']),sequence_fit_rows=seqrows,structural_fit_rows=structrows,
        sequence_groups=sgroups,comparison_groups=groups,pair_states=states,pair_status_counts=dict(counts),sequence_direction_counts=dict(seq_counts),joint_direction_counts=dict(joint_counts),
        both_screen_pass_counts=dict(passes),all_pair_screen_pass_group_counts=dict(all_passes),maximum_absolute_metric_difference=maxima)
    assert all(summary[k]==v for k,v in plan['expected'].items())
    assert set(summary)==set(SUMMARY_FIELDS)
    for field in SUMMARY_FIELDS:difference=max(difference,compare(receipt[field],summary[field]))
    bind(bindings,root/'configuration.json');assert receipt['source_hashes']==bindings
    assert set(receipt['artifacts'])==expected_artifacts|set(names)
    bind(bindings,rp,rh)
    for name,digest in receipt['artifacts'].items():bind(bindings,root/name,digest)
    verify(bindings)
    result=dict(status='passed_full_triad_correspondence_raw_msa_sql_readback',**summary,plan_sha256=sha(plan_path),producer_receipt_sha256=rh,
        maximum_absolute_numeric_difference=difference,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.output)
