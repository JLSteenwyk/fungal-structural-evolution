#!/usr/bin/env python3
"""Normalize all target/background orders with quarantined raw numerical failures."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
from full_expanded_measurement_catalog_sources import load,FIELDS,NUMERIC_FIELDS,GEOMETRY_FIELDS
from background_measurement_union_sources import NUMERIC_INTS,GEOMETRY_INTS
from reference_measurement_union_sources import keyed_table,verify
from run_ortholog_pair_guide_comparison import sha


def qualities(path,screens,count):
    records={}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=row['pair_key'],row['mask'];assert key not in records and row['mask'] in ['full','plddt70'];records[key]=row
            ends=[(row['model_'+s],int(row['version_'+s])) for s in ['a','b']]
            assert ends==sorted(ends) and ends[0]!=ends[1] and row['pair_key']==hashlib.sha256(json.dumps(ends,separators=(',',':')).encode()).hexdigest()
            for spec in screens:assert row[spec['id']+'_pass'] in ['0','1'] and (row[spec['id']+'_pass']=='1')==(row[spec['id']+'_exclusions']=='')
    pairs={p for p,m in records};assert len(pairs)==count and set(records)=={(p,m) for p in pairs for m in ['full','plddt70']}
    return records


def normalized(role,q,order,status,usable,why,numeric,geometry,source_kind,checkpoint,checkpoint_hash,source_order,other,screens):
    n=int(numeric['aligned_length']) if numeric is not None else None
    assert (status=='aligned')==(n is not None) and isinstance(usable,bool)
    assert usable==(not why) and (not usable or status=='aligned')
    row=dict(role=role,pair_key=q['pair_key'],mask=q['mask'],order=order,**{k:q[k] for k in ['model_a','version_a','model_b','version_b']},original_length_a=int(q['length_a']),original_length_b=int(q['length_b']),source_kind=source_kind,source_checkpoint=checkpoint,source_checkpoint_sha256=checkpoint_hash,source_order=source_order,native_status=status,numerical_usable=int(usable),numerical_exclusion_reasons=why,pair_mask_pass_bits=sum(1<<i for i,s in enumerate(screens) if q[s['id']+'_pass']=='1'),both_masks_pass_bits=sum(1<<i for i,s in enumerate(screens) if q[s['id']+'_pass']==other[s['id']+'_pass']=='1'))
    for side in ['a','b']:row['original_coverage_'+side]=n/int(q['length_'+side]) if n is not None else ''
    for names,data,ints in [(NUMERIC_FIELDS,numeric,NUMERIC_INTS),(GEOMETRY_FIELDS,geometry,GEOMETRY_INTS)]:
        for field in names:
            raw=data.get(field) if data is not None else None
            value='' if raw is None else str(raw) if field.endswith('_status') else int(raw) if field in ints else float(raw)
            if isinstance(value,(int,float)):assert math.isfinite(value)
            row[field]=value
    if n is not None:
        assert 0<n<=min(row['original_length_a'],row['original_length_b'])
        assert numeric['aligned_length']==geometry['aligned_length'] or int(numeric['aligned_length'])==int(geometry['aligned_length'])
    if usable:
        assert row['aligned_length']>=3 and row['rmsd_status']=='within_printed_rounding' and row['geometry_status']=='unique_at_numeric_tolerance'
    if row['pair_mask_pass_bits']:assert usable
    prefix=f'order{order}_'
    assert q[prefix+'native_status']==status
    assert float(q[prefix+'aligned_length'])==n if n is not None else q[prefix+'aligned_length']==''
    if role=='background':assert int(q[prefix+'numerical_usable'])==int(usable)
    else:assert (q[prefix+'status']=='aligned')==usable
    return row


def run(path):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path)
    out=Path(plan['output']);assert shutil.disk_usage(out.parent).free>=plan['resources']['minimum_free_disk_gib']*2**30;out.mkdir(exist_ok=True)
    lock=(out/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (out/'receipt.json').exists(),'Completed stage cannot be restarted'
    marker=out/'stage_plan.json';state=dict(plan_sha256=sha(path),schema='full-directed-measurement-catalog-v1')
    if marker.exists():assert json.loads(marker.read_text())==state
    else:marker.write_text(json.dumps(state,indent=2)+'\n')
    numeric=keyed_table(source['target_numeric']);geometry=keyed_table(source['target_geometry']);assert set(numeric)==set(geometry)
    summaries={}
    with source['target_summary'].open() as handle:
        for r in csv.DictReader(handle,delimiter='\t'):
            key=r['pair_key'],r['mask'];assert key not in summaries;summaries[key]=r
    tq=qualities(source['target_quality'],plan['screens'],plan['expected']['target_pairs']);assert set(tq)==set(summaries)
    native_pairs=set()
    with source['target_pairs'].open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            pair=row['pair_key'];assert pair not in native_pairs;native_pairs.add(pair)
            for mask in ['full','plddt70']:
                assert all(row[k]==tq[pair,mask][k] for k in ['model_a','version_a','model_b','version_b'])
    assert native_pairs=={p for p,m in tq}
    counts=Counter();sources=Counter();exclusions=Counter();total=raw_numeric=usable_total=0
    temp=out/'directed_measurements.tsv.gz.tmp'
    with gzip.open(temp,'wt',compresslevel=1) as handle:
        writer=csv.DictWriter(handle,FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader()
        def emit(row):
            nonlocal total,raw_numeric,usable_total
            writer.writerow(row);total+=1;raw_numeric+=row['aligned_length']!='';usable_total+=row['numerical_usable']
            counts[row['role']+':'+row['mask']+':'+row['native_status']]+=1;sources[row['role']+':'+row['source_kind']]+=1
            if row['numerical_exclusion_reasons']:exclusions[row['role']+':'+row['numerical_exclusion_reasons']]+=1
        seen=set()
        with source['target_checkpoints'].open() as native:
            for raw in csv.DictReader(native,delimiter='\t'):
                pair,mask,os=Path(raw['path']).stem.rsplit('-',2);order=int(os);key=pair,mask,order
                assert key not in seen and (pair,mask) in tq and order in [0,1];seen.add(key)
                checkpoint=source['target_native_root']/raw['path'];assert bindings[str(checkpoint)]==raw['sha256']
                q=tq[pair,mask];s=summaries[pair,mask];status=raw['status'];assert s[f'order{order}_native_status']==status
                usable=s[f'order{order}_status']=='aligned';why=s[f'order{order}_numerical_exclusion_reasons'] if status=='aligned' else 'native_status:'+status
                n=numeric.get(key);g=geometry.get(key);assert (status=='aligned')==(n is not None)==(g is not None)
                emit(normalized('target',q,order,status,usable,why,n,g,'target_original',str(checkpoint),raw['sha256'],order,tq[pair,'plddt70' if mask=='full' else 'full'],plan['screens']))
        assert len(seen)==4*plan['expected']['target_pairs'] and set(numeric)=={k for k in seen if k in numeric}
        assert len(numeric)==plan['expected']['target_numeric_states'] and usable_total==source['target_summary_readback']['numerically_usable_directions']
        del numeric,geometry,summaries,tq,seen
        bq=qualities(source['background_quality'],plan['screens'],plan['expected']['background_pairs']);seen=set()
        with gzip.open(source['background_states'],'rt') as handle:
            for line in handle:
                raw=json.loads(line);pair,mask,order=raw['pair_key'],raw['mask'],raw['order'];key=pair,mask,order
                assert key not in seen and (pair,mask) in bq and order in [0,1];seen.add(key);q=bq[pair,mask]
                ends=[(q['model_'+s],int(q['version_'+s])) for s in ['a','b']]
                assert [(e['model_id'],e['version']) for e in raw['directed_endpoints']]==(ends if order==0 else ends[::-1])
                assert bindings[raw['source_checkpoint']]==raw['source_checkpoint_sha256']
                emit(normalized('background',q,order,raw['source_native_status'],raw['numerical_usable'],';'.join(raw['numerical_exclusion_reasons']),raw['numerical'],raw['geometry'],raw['selected_source'],raw['source_checkpoint'],raw['source_checkpoint_sha256'],raw['source_order'],bq[pair,'plddt70' if mask=='full' else 'full'],plan['screens']))
        assert len(seen)==4*plan['expected']['background_pairs']
    assert total==4*(plan['expected']['target_pairs']+plan['expected']['background_pairs'])
    verify(bindings);temp.replace(out/'directed_measurements.tsv.gz')
    summary=dict(target_pairs=plan['expected']['target_pairs'],background_pairs=plan['expected']['background_pairs'],pair_mask_rows=total//2,directed_states=total,numeric_states=raw_numeric,usable_states=usable_total,unavailable_states=total-raw_numeric,native_status_counts=dict(counts),numerical_exclusion_counts=dict(exclusions),source_kind_counts=dict(sources),screens=plan['screens'],masks=['full','plddt70'])
    result=dict(status='complete_full_expanded_measurement_catalog_pending_independent_readback',plan_sha256=sha(path),**summary,source_hashes=bindings,artifacts={n:sha(out/n) for n in ['stage_plan.json','directed_measurements.tsv.gz']},scientific_eligibility=False,scope=plan['scope'])
    with (out/'receipt.json').open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
