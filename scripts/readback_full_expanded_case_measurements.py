#!/usr/bin/env python3
"""Independently check all identities, exclusions and Decimal order contrasts."""
import argparse
from collections import Counter
import csv
from decimal import Decimal,localcontext
import gzip
import hashlib
import json
from pathlib import Path
from full_expanded_case_measurement_sources import load,OUTCOMES,FIELDS,SUMMARY_FIELDS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def run(path,output):
    path=Path(path);plan=json.loads(path.read_text());sources,bindings=load(plan,path);original=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
    assert r['status']=='complete_full_expanded_case_measurement_join_pending_independent_readback' and r['plan_sha256']==sha(path) and r['scientific_eligibility'] is False and r['source_hashes']==original
    bind(bindings,rp);assert set(r['artifacts'])=={'stage_plan.json','case_mask_contrasts.tsv.gz'}
    for n,d in r['artifacts'].items():bind(bindings,root/n,d)
    verify(bindings);assert json.loads((root/'stage_plan.json').read_text())==dict(plan_sha256=sha(path),schema='expanded-matched-measurement-row-v1')
    measurements={}
    with gzip.open(sources['catalog']['root']/'directed_measurements.tsv.gz','rt') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=row['role'],row['pair_key'],row['mask'],int(row['order']);assert key not in measurements;measurements[key]=row
    assert len(measurements)==plan['expected']['directed_states']
    seen=set();count=same=0;valid=Counter();complete=Counter()
    with localcontext() as context,gzip.open(sources['case_index']['root']/'case_index.tsv.gz','rt') as cases,gzip.open(root/'case_mask_contrasts.tsv.gz','rt') as output_rows:
        context.prec=50;actual_rows=csv.DictReader(output_rows,delimiter='\t');assert actual_rows.fieldnames==FIELDS
        for case in csv.DictReader(cases,delimiter='\t'):
            assert case['case_id'] not in seen;seen.add(case['case_id'])
            for mask in ['full','plddt70']:
                actual=next(actual_rows,None);assert actual is not None
                expected={k:case[k] for k in ['case_id','physical_case_id','target_id','background_id','guide','target_sequence_distance','background_sequence_distance','target_same_model','background_same_model','target_pair_key','background_pair_key','target_both_pass_bits','control_both_pass_bits','joint_both_pass_bits']}
                expected['mask']=mask;expected['row_identity']=hashlib.sha256(json.dumps(['expanded-matched-measurement-row-v1',case['case_id'],mask],separators=(',',':')).encode()).hexdigest()
                for role in ['target','control','joint']:expected[f'{role}_mask_pass_bits']=case[f'{role}_{mask}_pass_bits']
                values={};numeric={};source_scale=Decimal(1)
                for role in ['target','background']:
                    bit_role='target' if role=='target' else 'control';is_same=case[role+'_same_model']=='1'
                    for order in [0,1]:
                        key=role,case[role+'_pair_key'],mask,order;prefix=f'{role}_order{order}_'
                        state=measurements.get(key);assert (state is None)==is_same
                        if is_same:
                            assert case[f'{bit_role}_{mask}_pass_bits']==case[f'{bit_role}_both_pass_bits']=='0'
                            expected.update({prefix+'state_key':'',prefix+'native_status':'identical_model_no_alignment',prefix+'numerical_usable':'0',prefix+'numerical_exclusion_reasons':'identical_model_no_alignment',prefix+'aligned_length':''});values[role,order]=None
                        else:
                            assert state['pair_mask_pass_bits']==case[f'{bit_role}_{mask}_pass_bits'] and state['both_masks_pass_bits']==case[f'{bit_role}_both_pass_bits']
                            expected[prefix+'state_key']='|'.join(str(x) for x in key)
                            for field in ['native_status','numerical_usable','numerical_exclusion_reasons','aligned_length']:expected[prefix+field]=state[field]
                            if state['numerical_usable']=='1':
                                assert state['native_status']=='aligned' and not state['numerical_exclusion_reasons']
                                a=Decimal(state['rmsd_recomputed']);left=Decimal(state['tm_left_native']);right=Decimal(state['tm_right_native'])
                                assert a.is_finite() and a>=0 and 0<=left<=1 and 0<=right<=1
                                values[role,order]={'rmsd':a,'native_tm_dissimilarity':Decimal(1)-(left+right)/Decimal(2)};source_scale=max(source_scale,a)
                            else:values[role,order]=None
                            if int(state['pair_mask_pass_bits']):assert values[role,order] is not None
                        for outcome in OUTCOMES:
                            field=prefix+outcome
                            if values[role,order] is None:expected[field]=''
                            else:numeric[field]=values[role,order][outcome]
                for outcome in OUTCOMES:
                    averages={}
                    for role in ['target','background']:
                        field=f'{outcome}_{role}_both_orders_mean'
                        if values[role,0] is not None and values[role,1] is not None:averages[role]=(values[role,0][outcome]+values[role,1][outcome])/Decimal(2);numeric[field]=averages[role]
                        else:expected[field]=''
                    field=outcome+'_both_orders_mean_delta'
                    if len(averages)==2:numeric[field]=averages['target']-averages['background']
                    else:expected[field]=''
                    matrix=[];bitmap=0
                    for index in range(4):
                        a,b=divmod(index,2);field=f'{outcome}_order_pair_{a}{b}_delta'
                        if values['target',a] is None or values['background',b] is None:expected[field]=''
                        else:
                            difference=values['target',a][outcome]-values['background',b][outcome];numeric[field]=difference;matrix.append(difference);bitmap+=2**index
                    expected[outcome+'_usable_order_pair_bits']=str(bitmap)
                    for name in ['min','max','span']:
                        field=outcome+'_complete_order_envelope_'+name
                        if bitmap!=15:expected[field]=''
                        else:numeric[field]=min(matrix) if name=='min' else max(matrix) if name=='max' else max(matrix)-min(matrix)
                    expected[outcome+'_disposition']='complete_four_order_pairs' if bitmap==15 else 'identical_model_no_alignment' if any(case[x+'_same_model']=='1' for x in ['target','background']) else 'incomplete_numerical_orders'
                    valid[outcome+':'+mask]+=bitmap.bit_count();complete[outcome+':'+mask]+=bitmap==15
                assert set(expected)|set(numeric)==set(FIELDS) and not set(expected)&set(numeric)
                for field,value in expected.items():assert actual[field]==value,(case['case_id'],mask,field)
                for field,value in numeric.items():
                    observed=Decimal(actual[field]);assert observed.is_finite()
                    scale=Decimal(1) if 'native_tm_dissimilarity' in field else source_scale
                    assert abs(observed-value)<=Decimal('2e-14')*max(scale,abs(value)),(case['case_id'],mask,field,observed,value)
                count+=1;same+=any(case[x+'_same_model']=='1' for x in ['target','background'])
                if count%50000==0:print('independent_full_case_mask_contrasts',count,flush=True)
        assert next(actual_rows,None) is None
    assert len(seen)==plan['expected']['logical_cases'] and count==2*len(seen)
    summary=dict(logical_cases=len(seen),physical_cases=sources['case_index']['receipt']['physical_cases'],selected_records=plan['expected']['selected_records'],unmatched_decisions=plan['expected']['unmatched_decisions'],case_mask_rows=count,order_pair_cells=count*4*len(OUTCOMES),valid_order_pair_cells=dict(valid),complete_envelopes=dict(complete),incomplete_envelopes={k:len(seen)-v for k,v in complete.items()},same_model_case_mask_rows=same,outcomes=OUTCOMES,masks=['full','plddt70'])
    assert all(summary[k]==r[k] for k in SUMMARY_FIELDS);verify(bindings)
    result=dict(status='passed_full_expanded_case_measurement_join_decimal_readback',plan_sha256=sha(path),producer_receipt_sha256=sha(rp),**summary,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.plan,a.output)
