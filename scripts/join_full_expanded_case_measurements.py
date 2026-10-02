#!/usr/bin/env python3
"""Join all logical cases and both masks; retain every target/control order pair."""
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
from full_expanded_case_measurement_sources import load,OUTCOMES,STATES,STATE_FIELDS,BASE_FIELDS,FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def metrics(state):
    if state is None:return None
    if state['numerical_usable']!='1':return None
    assert state['native_status']=='aligned' and not state['numerical_exclusion_reasons']
    values=dict(rmsd=float(state['rmsd_recomputed']),native_tm_dissimilarity=1-(float(state['tm_left_native'])+float(state['tm_right_native']))/2)
    assert all(math.isfinite(v) for v in values.values()) and values['rmsd']>=0 and 0<=values['native_tm_dissimilarity']<=1
    return values


def joined(case,mask,catalog):
    row={k:case[k] for k in BASE_FIELDS if k in case};row['mask']=mask
    row['row_identity']=hashlib.sha256(json.dumps(['expanded-matched-measurement-row-v1',case['case_id'],mask],separators=(',',':')).encode()).hexdigest()
    for role in ['target','control','joint']:row[f'{role}_mask_pass_bits']=case[f'{role}_{mask}_pass_bits']
    states={};values={}
    for role,order in STATES:
        key=role,case['target_pair_key' if role=='target' else 'background_pair_key'],mask,order
        same=case[role+'_same_model']=='1';state=None if same else catalog[key];states[role,order]=state;values[role,order]=metrics(state)
        if same:assert key not in catalog
        source_values=dict(state_key='' if same else '|'.join(map(str,key)),native_status='identical_model_no_alignment' if same else state['native_status'],numerical_usable='0' if same else state['numerical_usable'],numerical_exclusion_reasons='identical_model_no_alignment' if same else state['numerical_exclusion_reasons'],aligned_length='' if same else state['aligned_length'])
        source_values.update({o:'' if values[role,order] is None else values[role,order][o] for o in OUTCOMES})
        row.update({f'{role}_order{order}_{f}':source_values[f] for f in STATE_FIELDS})
        bits=int(case[f"{'target' if role=='target' else 'control'}_{mask}_pass_bits"])
        both=int(case[f"{'target' if role=='target' else 'control'}_both_pass_bits"])
        if same:assert bits==both==0
        else:assert int(state['pair_mask_pass_bits'])==bits and int(state['both_masks_pass_bits'])==both
        if bits:assert values[role,order] is not None
    for outcome in OUTCOMES:
        means={r:sum(values[r,o][outcome] for o in [0,1])/2 if all(values[r,o] is not None for o in [0,1]) else None for r in ['target','background']}
        row[outcome+'_target_both_orders_mean']=means['target'] if means['target'] is not None else ''
        row[outcome+'_background_both_orders_mean']=means['background'] if means['background'] is not None else ''
        row[outcome+'_both_orders_mean_delta']=means['target']-means['background'] if None not in means.values() else ''
        deltas=[];bits=0
        for i,(a,b) in enumerate([(0,0),(0,1),(1,0),(1,1)]):
            value=values['target',a][outcome]-values['background',b][outcome] if values['target',a] is not None and values['background',b] is not None else None
            row[f'{outcome}_order_pair_{a}{b}_delta']='' if value is None else value
            if value is not None:bits|=1<<i;deltas.append(value)
        row[outcome+'_usable_order_pair_bits']=bits
        for f,value in [('min',min(deltas) if bits==15 else ''),('max',max(deltas) if bits==15 else ''),('span',max(deltas)-min(deltas) if bits==15 else '')]:row[outcome+'_complete_order_envelope_'+f]=value
        row[outcome+'_disposition']='complete_four_order_pairs' if bits==15 else 'identical_model_no_alignment' if any(case[r+'_same_model']=='1' for r in ['target','background']) else 'incomplete_numerical_orders'
    return row


def run(path):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);root=Path(plan['output'])
    assert shutil.disk_usage(root.parent).free>=plan['resources']['minimum_free_disk_gib']*2**30;root.mkdir(exist_ok=True)
    lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (root/'receipt.json').exists(),'Completed join cannot be restarted'
    marker=root/'stage_plan.json';s=dict(plan_sha256=sha(path),schema='expanded-matched-measurement-row-v1')
    if marker.exists():assert json.loads(marker.read_text())==s
    else:marker.write_text(json.dumps(s,indent=2)+'\n')
    catalog={}
    with gzip.open(source['catalog']['root']/'directed_measurements.tsv.gz','rt') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=row['role'],row['pair_key'],row['mask'],int(row['order']);assert key not in catalog;catalog[key]=row
    assert len(catalog)==plan['expected']['directed_states']
    seen=set();rows=0;complete=Counter();valid=Counter();same=0;temp=root/'case_mask_contrasts.tsv.gz.tmp'
    with gzip.open(source['case_index']['root']/'case_index.tsv.gz','rt') as handle,gzip.open(temp,'wt',compresslevel=1) as output:
        writer=csv.DictWriter(output,FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader()
        for case in csv.DictReader(handle,delimiter='\t'):
            assert case['case_id'] not in seen;seen.add(case['case_id'])
            for mask in ['full','plddt70']:
                row=joined(case,mask,catalog);writer.writerow(row);rows+=1;same+=any(case[r+'_same_model']=='1' for r in ['target','background'])
                for outcome in OUTCOMES:
                    bits=row[outcome+'_usable_order_pair_bits'];valid[outcome+':'+mask]+=bits.bit_count();complete[outcome+':'+mask]+=bits==15
    assert len(seen)==plan['expected']['logical_cases'] and rows==2*len(seen)
    verify(bindings);temp.replace(root/'case_mask_contrasts.tsv.gz')
    summary=dict(logical_cases=len(seen),physical_cases=source['case_index']['receipt']['physical_cases'],selected_records=plan['expected']['selected_records'],unmatched_decisions=plan['expected']['unmatched_decisions'],case_mask_rows=rows,order_pair_cells=rows*4*len(OUTCOMES),valid_order_pair_cells=dict(valid),complete_envelopes=dict(complete),incomplete_envelopes={k:len(seen)-v for k,v in complete.items()},same_model_case_mask_rows=same,outcomes=OUTCOMES,masks=['full','plddt70'])
    r=dict(status='complete_full_expanded_case_measurement_join_pending_independent_readback',plan_sha256=sha(path),**summary,source_hashes=bindings,artifacts={n:sha(root/n) for n in ['stage_plan.json','case_mask_contrasts.tsv.gz']},scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as handle:handle.write(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:v for k,v in r.items() if k!='source_hashes'},indent=2),flush=True);return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
