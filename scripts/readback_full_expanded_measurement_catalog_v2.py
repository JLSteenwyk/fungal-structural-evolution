#!/usr/bin/env python3
"""Reconstruct every normalized state directly from the complete closed sources."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from full_expanded_measurement_catalog_sources_v2 import load,FIELDS,NUMERIC_FIELDS,GEOMETRY_FIELDS,SUMMARY_FIELDS
from background_measurement_union_sources import NUMERIC_INTS,GEOMETRY_INTS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def run(path, output):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);original_bindings=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
    assert r['status']=='complete_full_expanded_measurement_catalog_pending_independent_readback' and r['plan_sha256']==sha(path) and r['scientific_eligibility'] is False and r['source_hashes']==original_bindings
    bind(bindings,rp);assert set(r['artifacts'])=={'stage_plan.json','directed_measurements.tsv.gz'}
    for n,d in r['artifacts'].items():bind(bindings,root/n,d)
    verify(bindings);assert json.loads((root/'stage_plan.json').read_text())==dict(plan_sha256=sha(path),schema='full-directed-measurement-catalog-v1')
    quality={};pairs={}
    for role,key,count in [('target','target_quality','target_pairs'),('background','background_quality','background_pairs')]:
        records={}
        with source[key].open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                k=row['pair_key'],row['mask'];assert k not in records and row['mask'] in ['full','plddt70'];records[k]=row
                endpoints=[(row['model_'+s],int(row['version_'+s])) for s in ['a','b']]
                assert endpoints[0]<endpoints[1] and hashlib.sha256(json.dumps(endpoints,separators=(',',':')).encode()).hexdigest()==row['pair_key']
                for spec in plan['screens']:assert row[spec['id']+'_pass'] in ['0','1'] and (row[spec['id']+'_pass']=='1')==(not row[spec['id']+'_exclusions'])
        pairset={p for p,m in records};assert len(pairset)==plan['expected'][count] and set(records)=={(p,m) for p in pairset for m in ['full','plddt70']};quality[role]=records;pairs[role]=pairset
    queue=set()
    with source['target_pairs'].open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            assert row['pair_key'] not in queue;queue.add(row['pair_key'])
            for mask in ['full','plddt70']:assert all(row[k]==quality['target'][row['pair_key'],mask][k] for k in ['model_a','version_a','model_b','version_b'])
    assert queue==pairs['target']
    numeric={};geometry={}
    for dest,key in [(numeric,'target_numeric'),(geometry,'target_geometry')]:
        with source[key].open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                k=row['pair_key'],row['mask'],int(row['order']);assert k not in dest;dest[k]=row
    assert set(numeric)==set(geometry) and len(numeric)==plan['expected']['target_numeric_states']
    total=raw_numeric=usable_total=target_usable=0;seen={'target':set(),'background':set()};statuscounts=Counter();sourcecounts=Counter();exclusions=Counter()
    with gzip.open(root/'directed_measurements.tsv.gz','rt') as handle:
        actual_rows=csv.DictReader(handle,delimiter='\t');assert actual_rows.fieldnames==FIELDS
        def check(role,pair,mask,order,status,usable,why,number,shape,kind,checkpoint,cp_hash,source_order):
            nonlocal total,raw_numeric,usable_total,target_usable
            key=pair,mask,order;assert key not in seen[role] and pair in pairs[role] and order in [0,1];seen[role].add(key)
            q=quality[role][pair,mask];other=quality[role][pair,'full' if mask=='plddt70' else 'plddt70']
            n=None if number is None else int(number['aligned_length'])
            assert (number is not None)==(shape is not None)==(status=='aligned') and bool(usable)==(not why)
            expected={'role':role,'pair_key':pair,'mask':mask,'order':str(order),**{k:q[k] for k in ['model_a','version_a','model_b','version_b']},'original_length_a':q['length_a'],'original_length_b':q['length_b'],'source_kind':kind,'source_checkpoint':str(checkpoint),'source_checkpoint_sha256':cp_hash,'source_order':str(source_order),'native_status':status,'numerical_usable':str(int(usable)),'numerical_exclusion_reasons':why}
            flag=intersection=0
            for i,spec in enumerate(plan['screens']):
                if int(q[spec['id']+'_pass']):flag|=2**i
                if int(q[spec['id']+'_pass']) and int(other[spec['id']+'_pass']):intersection|=2**i
            expected.update(pair_mask_pass_bits=str(flag),both_masks_pass_bits=str(intersection))
            if flag:assert usable
            for side in ['a','b']:expected['original_coverage_'+side]='' if n is None else str(n/int(q['length_'+side]))
            for fields,data in [(NUMERIC_FIELDS,number),(GEOMETRY_FIELDS,shape)]:
                for field in fields:
                    value=None if data is None else data.get(field)
                    if value is None:text=''
                    elif field.endswith('_status'):text=str(value)
                    elif field in NUMERIC_INTS+GEOMETRY_INTS:text=str(int(value))
                    else:
                        value=float(value);assert math.isfinite(value);text=str(value)
                    expected[field]=text
            if n is not None:
                assert n==int(shape['aligned_length']) and 0<n<=min(int(q['length_a']),int(q['length_b']))
            if usable:assert n>=3 and number['rmsd_status']=='within_printed_rounding' and shape['geometry_status']=='unique_at_numeric_tolerance'
            assert q[f'order{order}_native_status']==status
            assert float(q[f'order{order}_aligned_length'])==n if n is not None else q[f'order{order}_aligned_length']==''
            if role=='target':assert (q[f'order{order}_status']=='aligned')==bool(usable)
            else:assert int(q[f'order{order}_numerical_usable'])==int(usable)
            actual=next(actual_rows,None);assert actual is not None and actual==expected,(role,pair,mask,order)
            assert bindings[str(checkpoint)]==cp_hash
            total+=1;raw_numeric+=n is not None;usable_total+=bool(usable);target_usable+=bool(usable) and role=='target'
            statuscounts[role+':'+mask+':'+status]+=1;sourcecounts[role+':'+kind]+=1
            if why:exclusions[role+':'+why]+=1
            if total%100000==0:print('independent_full_expanded_measurement_states',total,flush=True)
        with source['target_checkpoints'].open() as native:
            for row in csv.DictReader(native,delimiter='\t'):
                pair,mask,os=Path(row['path']).stem.rsplit('-',2);order=int(os);key=pair,mask,order;status=row['status'];number,shape=numeric.get(key),geometry.get(key)
                reasons=[]
                if status!='aligned':reasons=['native_status:'+status]
                else:
                    assert number is not None and shape is not None
                    if number['rmsd_status']!='within_printed_rounding':reasons.append('rmsd_discrepancy')
                    if int(number['aligned_length'])<3:reasons.append('fewer_than_three_pairs')
                    if shape['geometry_status']!='unique_at_numeric_tolerance':reasons.append('nonunique_rotation')
                check('target',pair,mask,order,status,not reasons,';'.join(reasons),number,shape,'target_original',source['target_native_root']/row['path'],row['sha256'],order)
        assert set(numeric)<=seen['target'] and target_usable==source['target_summary_readback']['numerically_usable_directions']
        del numeric,geometry
        with gzip.open(source['background_states'],'rt') as native:
            for line in native:
                row=json.loads(line);pair,mask,order=row['pair_key'],row['mask'],row['order'];q=quality['background'][pair,mask]
                endpoints=[(q['model_'+s],int(q['version_'+s])) for s in ['a','b']]
                assert [(e['model_id'],e['version']) for e in row['directed_endpoints']]==(endpoints if order==0 else list(reversed(endpoints)))
                reasons=row['numerical_exclusion_reasons'];number,shape=row['numerical'],row['geometry']
                if row['source_native_status']=='aligned':
                    recomputed=[]
                    if number['rmsd_status']!='within_printed_rounding':recomputed.append('rmsd_discrepancy')
                    if number['aligned_length']<3:recomputed.append('fewer_than_three_pairs')
                    if shape['geometry_status']!='unique_at_numeric_tolerance':recomputed.append('nonunique_rotation')
                    assert reasons==recomputed
                else:assert reasons==[row['source_native_status']] and not row['numerical_usable']
                check('background',pair,mask,order,row['source_native_status'],row['numerical_usable'],';'.join(reasons),number,shape,row['selected_source'],row['source_checkpoint'],row['source_checkpoint_sha256'],row['source_order'])
        assert next(actual_rows,None) is None
    assert all(seen[role]=={(p,m,o) for p in pairs[role] for m in ['full','plddt70'] for o in [0,1]} for role in seen)
    summary=dict(target_pairs=len(pairs['target']),background_pairs=len(pairs['background']),pair_mask_rows=total//2,directed_states=total,numeric_states=raw_numeric,usable_states=usable_total,unavailable_states=total-raw_numeric,native_status_counts=dict(statuscounts),numerical_exclusion_counts=dict(exclusions),source_kind_counts=dict(sourcecounts),screens=plan['screens'],masks=['full','plddt70'])
    assert all(summary[k]==r[k] for k in SUMMARY_FIELDS);verify(bindings)
    result=dict(status='passed_full_expanded_measurement_catalog_source_and_numeric_readback',plan_sha256=sha(path),producer_receipt_sha256=sha(rp),**summary,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.plan,a.output)
