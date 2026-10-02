#!/usr/bin/env python3
"""Build every case/mask/order design row and every original matching stratum."""
import argparse
from collections import Counter, defaultdict
import csv
import fcntl
import gzip
import itertools
import json
import math
from pathlib import Path
import shutil
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from full_expanded_model_input_sources import load, schema, identity, FIELDS, NUMERIC, ORDERS, MASKS, GATES, COUNT_FIELDS, SUMMARY_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def static_features(case, nodes):
    a, b = [float(case[r + '_sequence_distance']) for r in ['target', 'background']]
    assert math.isfinite(a) and math.isfinite(b) and min(a, b) >= 0
    result = dict(zip(['gene_distance_delta', 'gene_distance_square_delta', 'gene_distance_cube_delta'], [a**k-b**k for k in [1,2,3]]))
    sides = {}
    for role in ['target', 'background']:
        n = nodes[role][case[role + '_id']]
        assert all(n['length_' + e] > 0 for e in ['a','b'])
        sides[role] = [sum(math.log(n['length_' + e]) for e in ['a','b'])/2,
            sum(n['mean_ca_plddt_' + e] for e in ['a','b'])/2,
            sum(n['fraction_ca_plddt_below50_' + e] for e in ['a','b'])/2]
        assert all(math.isfinite(v) for v in sides[role]) and 0<=sides[role][1]<=100 and 0<=sides[role][2]<=1
    result.update(zip(['log_original_length_ratio','mean_ca_plddt_delta','fraction_ca_plddt_below50_delta'],
        [a-b for a,b in zip(sides['target'],sides['background'])]))
    assert all(math.isfinite(v) for v in result.values())
    return result


def model_row(case, case_row, contrast, cov, order, source):
    mask = contrast['mask']; result = {k:case[k] for k in FIELDS if k in case}
    result.update({k:cov[k] for k in ['family_component','species_pattern_id']})
    result.update(input_id=identity(contrast['row_identity'],order),row_identity=contrast['row_identity'],case_row=case_row,
        species_pattern_row=int(cov['species_pattern_row']),mask=mask,order_contrast=order,
        joint_mask_pass_bits=int(contrast['joint_mask_pass_bits']),joint_both_pass_bits=int(contrast['joint_both_pass_bits']),
        selection_records=int(case['selection_records']))
    result.update(static_features(case,source['nodes']))
    selected = {'target':[0,1] if order=='mean' else [int(order[0])], 'background':[0,1] if order=='mean' else [int(order[1])]}
    usable = True; side_values = {}
    for role, indices in selected.items():
        states = [None if case[role+'_same_model']=='1' else source['catalog'][role,case[role+'_pair_key'],mask,i] for i in indices]
        result[role+'_state_keys'] = json.dumps([contrast[f'{role}_order{i}_state_key'] for i in indices],separators=(',',':'))
        if any(s is None or s['numerical_usable']!='1' for s in states): usable=False; continue
        values = []
        node=source['nodes'][role][case[role+'_id']]
        for s in states:
            assert s['native_status']=='aligned'
            original_ends=sorted((node['model_id_'+e],int(node['version_'+e]),int(node['length_'+e])) for e in ['a','b'])
            measured_ends=sorted((s['model_'+e],int(s['version_'+e]),int(s['original_length_'+e])) for e in ['a','b'])
            assert original_ends==measured_ends
            i=float(s['sequence_identity_exact']);f=float(s['joint_plddt70_fraction']);n=int(s['aligned_length'])
            coverage=(float(s['original_coverage_a'])+float(s['original_coverage_b']))/2
            assert 0<=i<=1 and 0<=f<=1 and 0<=coverage<=1 and n>0
            values.append([i,i*i,i*i*i,coverage,math.log(n),f])
        side_values[role]=np.mean(values,axis=0)
    result['numerical_usable']=int(usable)
    result['disposition']='usable_order_contrast' if usable else 'identical_model_no_alignment' if any(case[r+'_same_model']=='1' for r in selected) else 'unavailable_or_quarantined_order_state'
    state_columns=['aligned_identity_delta','aligned_identity_square_delta','aligned_identity_cube_delta','original_coverage_delta','log_aligned_length_ratio','aligned_plddt70_fraction_delta']
    result.update(zip(state_columns,side_values['target']-side_values['background'] if usable else [None]*6))
    for outcome in ['rmsd','native_tm_dissimilarity']:
        key=outcome+('_both_orders_mean_delta' if order=='mean' else '_order_pair_'+order+'_delta')
        assert bool(contrast[key])==usable
        result[outcome+'_delta']=float(contrast[key]) if usable else None
    if not usable: assert result['joint_mask_pass_bits']==result['joint_both_pass_bits']==0
    assert set(result)==set(FIELDS)
    assert all(v is None or math.isfinite(v) for k,v in result.items() if k in NUMERIC)
    return result


def grouped_selections(source, plan):
    index={c['case_id']:i for i,c in enumerate(source['cases'])};groups=defaultdict(list);reuse=np.zeros(len(index),dtype=np.int64)
    p=source['stages']['cases']['root']/'selection_case_links.tsv.gz'
    with gzip.open(p,'rt') as f:
        for ordinal,row in enumerate(csv.DictReader(f,delimiter='\t'),1):
            assert int(row['source_row_ordinal'])==ordinal
            i=index[row['case_id']];case=source['cases'][i]
            for field in ['physical_case_id','target_id','background_id','guide']:assert row[field]==case[field]
            groups[row['guide'],row['policy'],row['scenario_id']].append(i);reuse[i]+=1
    assert ordinal==plan['expected']['selected_records'] and np.array_equal(reuse,[int(c['selection_records']) for c in source['cases']])
    return {k:np.asarray(v,dtype=np.int32) for k,v in groups.items()}


def stratum_counts(source,plan,valid):
    groups=grouped_selections(source,plan);cases=source['cases'];cov=source['cov']
    arrays={}
    for name,values in [('target',[c['target_id'] for c in cases]),('background',[c['background_id'] for c in cases]),
        ('background_pair',[c['background_pair_key'] for c in cases]),('family',[c['target_family'] for c in cases]),
        ('component',[cov[c['case_id']]['family_component'] for c in cases]),('pattern',[cov[c['case_id']]['species_pattern_id'] for c in cases]),('taxon',[c['focal_taxon'] for c in cases])]:
        arrays[name]=np.unique(values,return_inverse=True)[1]
    rows=[];same=np.asarray([any(c[r+'_same_model']=='1' for r in ['target','background']) for c in cases])
    for guide,policy,scenario in itertools.product(plan['guides'],plan['policies'],[s['scenario_id'] for s in source['scenarios']]):
        members=groups.get((guide,policy,scenario),np.empty(0,dtype=np.int32))
        assert len(np.unique(arrays['target'][members]))==len(members),'Original matching permits only one selected control per target/stratum'
        for mask,gate in itertools.product(MASKS,GATES):
            source_mask=mask if gate=='mask' else 'both'
            bits=np.asarray([int(c['joint_'+source_mask+'_pass_bits']) for c in cases],dtype=np.uint8)
            for bit,screen in enumerate(plan['screens']):
                retained=members[(bits[members] & (1<<bit))!=0]
                attrs=source['attrition'][guide,policy,scenario,source_mask,screen['id']]
                assert len(members)==int(attrs['matched_records']) and len(retained)==int(attrs['joint_pass_matched_records'])
                counts={};maxima={}
                for name in arrays:
                    ids,n=np.unique(arrays[name][retained],return_counts=True);counts[name]=len(ids);maxima[name]=int(n.max(initial=0))
                for order in ORDERS:
                    number=int(valid[mask,order][retained].sum());assert number==len(retained)
                    rows.append(dict(guide=guide,policy=policy,scenario_id=scenario,mask=mask,eligibility_gate=gate,screen=screen['id'],order_contrast=order,
                        all_target_records=int(attrs['all_target_records']),unmatched_records=int(attrs['unmatched_records']),selected_records=len(members),same_model_records=int(same[members].sum()),quality_retained_records=len(retained),quality_excluded_records=len(members)-len(retained),
                        quality_retained_numeric_records=number,unique_targets=counts['target'],unique_background_nodes=counts['background'],
                        unique_background_physical_pairs=counts['background_pair'],unique_target_families=counts['family'],unique_family_components=counts['component'],
                        unique_species_patterns=counts['pattern'],unique_focal_taxa=counts['taxon'],maximum_background_node_reuse=maxima['background'],maximum_background_physical_pair_reuse=maxima['background_pair'],
                        sum_inverse_background_node_reuse_weights=counts['background'],sum_inverse_background_physical_pair_reuse_weights=counts['background_pair'],sum_inverse_family_component_weights=counts['component']))
    return rows


def run(path, stop_after_cases=None):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);root=Path(plan['output'])
    assert shutil.disk_usage(root.parent).free>=plan['resources']['minimum_free_disk_gib']*2**30;root.mkdir(exist_ok=True)
    lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (root/'receipt.json').exists(),'Completed input stage cannot restart'
    marker=root/'stage_plan.json';state=dict(plan_sha256=sha(path),schema='expanded-model-input-v1')
    if marker.exists():assert json.loads(marker.read_text())==state
    else:marker.write_text(json.dumps(state,indent=2)+'\n')
    folder=root/'inputs';folder.mkdir(exist_ok=True);writers={};buffers={};valid={};parts=[]
    for mask,order in itertools.product(MASKS,ORDERS):
        key=mask,order;tmp=folder/(mask+'-'+order+'.parquet.tmp')
        writers[key]=pq.ParquetWriter(tmp,schema(),compression='zstd');buffers[key]=[];valid[key]=np.zeros(len(source['cases']),dtype=bool)
    def flush():
        for key,rows in buffers.items():
            if rows:writers[key].write_table(pa.Table.from_pylist(rows,schema=schema()));rows.clear()
    contrasts=source['stages']['measurements']['root']/'case_mask_contrasts.tsv.gz'
    try:
        with gzip.open(contrasts,'rt') as f:
            actual=csv.DictReader(f,delimiter='\t')
            for i,case in enumerate(source['cases']):
                for mask in MASKS:
                    row=next(actual,None);assert row and row['case_id']==case['case_id'] and row['mask']==mask
                    for order in ORDERS:
                        model=model_row(case,i,row,source['cov'][case['case_id']],order,source)
                        buffers[mask,order].append(model);valid[mask,order][i]=model['numerical_usable']
                if (i+1)%500==0:flush()
                if stop_after_cases is not None and i+1==stop_after_cases:raise InterruptedError('Software interruption contract')
            assert next(actual,None) is None
        flush()
    finally:
        for w in writers.values():w.close()
    counts=stratum_counts(source,plan,valid)
    expected_settings=len(plan['guides'])*len(plan['policies'])*plan['expected']['scenarios']*2*2*len(plan['screens'])*5
    assert len(counts)==expected_settings
    for mask,order in itertools.product(MASKS,ORDERS):
        tmp=folder/(mask+'-'+order+'.parquet.tmp');tmp.replace(folder/(mask+'-'+order+'.parquet'))
        p=folder/(mask+'-'+order+'.parquet');parts.append(dict(mask=mask,order_contrast=order,path=str(p.relative_to(root)),rows=len(source['cases']),sha256=sha(p)))
    mp=root/'partition_manifest.json';mp.write_text(json.dumps(parts,indent=2)+'\n')
    cp=root/'setting_counts.tsv'
    with cp.open('w') as f:w=csv.DictWriter(f,COUNT_FIELDS,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(counts)
    verify(bindings)
    rows=len(source['cases'])*10;valid_counts={'|'.join(k):int(v.sum()) for k,v in valid.items()}
    summary=dict(logical_cases=len(source['cases']),physical_cases=plan['expected']['physical_cases'],selected_records=plan['expected']['selected_records'],unmatched_decisions=plan['expected']['unmatched_decisions'],
        case_mask_rows=2*len(source['cases']),model_input_rows=rows,partitions=10,settings=len(counts),valid_input_rows=valid_counts,invalid_input_rows={k:len(source['cases'])-v for k,v in valid_counts.items()},
        quality_retained_setting_occurrences=sum(r['quality_retained_numeric_records'] for r in counts),model_specifications_per_setting=12,future_model_setting_rows=len(counts)*12,working_tree_alternatives=5)
    names=['stage_plan.json','partition_manifest.json','setting_counts.tsv']+[p['path'] for p in parts]
    receipt=dict(status='complete_full_expanded_model_inputs_pending_independent_readback',plan_sha256=sha(path),**summary,
        source_hashes=bindings,artifacts={n:sha(root/n) for n in names},scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as f:f.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(summary),flush=True);return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
