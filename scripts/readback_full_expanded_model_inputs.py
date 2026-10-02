#!/usr/bin/env python3
"""Independently reconstruct every expanded input cell and SQLite cohort count."""
import argparse
from collections import Counter
import csv
from decimal import Decimal, localcontext
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3
import numpy as np
import pyarrow.parquet as pq
from full_expanded_model_input_sources import load, schema, identity, FIELDS, NUMERIC, ORDERS, MASKS, GATES, COUNT_FIELDS, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

D=Decimal


def original_static(case,source):
    t,b=[D(case[r+'_sequence_distance']) for r in ['target','background']]
    numbers={};scales={}
    for k,name in zip([1,2,3],['gene_distance_delta','gene_distance_square_delta','gene_distance_cube_delta']):
        numbers[name]=t**k-b**k;scales[name]=max(D(1),abs(t**k),abs(b**k))
    pair={}
    for role in ['target','background']:
        n=source['nodes'][role][case[role+'_id']]
        pair[role]=[(D(n['length_a']).ln()+D(n['length_b']).ln())/2,
            (D(str(n['mean_ca_plddt_a']))+D(str(n['mean_ca_plddt_b'])))/2,
            (D(str(n['fraction_ca_plddt_below50_a']))+D(str(n['fraction_ca_plddt_below50_b'])))/2]
    for i,name in enumerate(['log_original_length_ratio','mean_ca_plddt_delta','fraction_ca_plddt_below50_delta']):
        numbers[name]=pair['target'][i]-pair['background'][i];scales[name]=max(D(1),abs(pair['target'][i]),abs(pair['background'][i]))
    return numbers,scales


def expected_state(state):
    if state is None or state['numerical_usable']!='1':return None
    i=D(state['sequence_identity_exact']);rmsd=D(state['rmsd_recomputed'])
    values=[rmsd,D(1)-(D(state['tm_left_native'])+D(state['tm_right_native']))/2,
        i,i**2,i**3,(D(state['original_coverage_a'])+D(state['original_coverage_b']))/2,
        D(state['aligned_length']).ln(),D(state['joint_plddt70_fraction'])]
    assert all(v.is_finite() for v in values)
    return values


def sqlite_counts(source,plan,database):
    assert not database.exists();db=sqlite3.connect(database);db.execute('PRAGMA foreign_keys=ON');db.execute('PRAGMA cache_size=-65536')
    db.execute('CREATE TABLE cases (case_id TEXT PRIMARY KEY, guide TEXT, target TEXT, background TEXT, pair TEXT, family TEXT, component TEXT, pattern TEXT, taxon TEXT, same INTEGER, full INTEGER, plddt70 INTEGER, both INTEGER)')
    cs=source['cases'];cov=source['cov'];index={c['case_id']:c for c in cs}
    db.executemany('INSERT INTO cases VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',[(c['case_id'],c['guide'],c['target_id'],c['background_id'],c['background_pair_key'],c['target_family'],cov[c['case_id']]['family_component'],cov[c['case_id']]['species_pattern_id'],c['focal_taxon'],int(c['target_same_model']=='1' or c['background_same_model']=='1'),int(c['joint_full_pass_bits']),int(c['joint_plddt70_pass_bits']),int(c['joint_both_pass_bits'])) for c in cs])
    db.execute('CREATE TABLE selections (ordinal INTEGER PRIMARY KEY,case_id TEXT REFERENCES cases(case_id),guide TEXT,policy TEXT,scenario_id TEXT,target TEXT,UNIQUE(guide,policy,scenario_id,target))')
    batch=[];seen=Counter()
    with gzip.open(source['stages']['cases']['root']/'selection_case_links.tsv.gz','rt') as f:
        for number,r in enumerate(csv.DictReader(f,delimiter='\t'),1):
            assert int(r['source_row_ordinal'])==number;c=index[r['case_id']]
            for name in ['guide','physical_case_id','target_id','background_id']:assert r[name]==c[name]
            seen[r['case_id']]+=1;batch.append((number,r['case_id'],r['guide'],r['policy'],r['scenario_id'],r['target_id']))
            if len(batch)==10000:db.executemany('INSERT INTO selections VALUES (?,?,?,?,?,?)',batch);batch=[]
    if batch:db.executemany('INSERT INTO selections VALUES (?,?,?,?,?,?)',batch)
    assert number==plan['expected']['selected_records'] and seen=={c['case_id']:int(c['selection_records']) for c in cs}
    db.execute('CREATE INDEX selection_case ON selections(case_id)');db.commit()
    common='FROM selections s JOIN cases c ON c.case_id=s.case_id'
    keycols='s.guide,s.policy,s.scenario_id'
    baseline={tuple(r[:3]):(r[3],r[4]) for r in db.execute('SELECT '+keycols+',COUNT(*),SUM(c.same) '+common+' GROUP BY '+keycols)}
    expected={}
    for mask in ['full','plddt70','both']:
        for bit,screen in enumerate(plan['screens']):
            where=' WHERE (c.'+mask+' & ?)<>0'
            values={tuple(r[:3]):list(r[3:]) for r in db.execute('SELECT '+keycols+',COUNT(*),COUNT(DISTINCT c.target),COUNT(DISTINCT c.background),COUNT(DISTINCT c.pair),COUNT(DISTINCT c.family),COUNT(DISTINCT c.component),COUNT(DISTINCT c.pattern),COUNT(DISTINCT c.taxon) '+common+where+' GROUP BY '+keycols,(1<<bit,))}
            maxima={}
            for field in ['background','pair']:
                query='WITH reuse AS (SELECT '+keycols+',c.'+field+' entity,COUNT(*) n '+common+where+' GROUP BY '+keycols+',c.'+field+') SELECT guide,policy,scenario_id,MAX(n) FROM reuse GROUP BY guide,policy,scenario_id'
                maxima[field]={tuple(r[:3]):r[3] for r in db.execute(query,(1<<bit,))}
            for group in itertools.product(plan['guides'],plan['policies'],[s['scenario_id'] for s in source['scenarios']]):
                selected,same=baseline.get(group,(0,0));counts=values.get(group,[0]*8);attrs=source['attrition'][group+(mask,screen['id'])]
                assert selected==int(attrs['matched_records']) and counts[0]==int(attrs['joint_pass_matched_records'])
                expected[group+(mask,screen['id'])]=dict(all_target_records=int(attrs['all_target_records']),unmatched_records=int(attrs['unmatched_records']),selected_records=selected,same_model_records=same,
                    quality_retained_records=counts[0],quality_excluded_records=selected-counts[0],quality_retained_numeric_records=counts[0],unique_targets=counts[1],unique_background_nodes=counts[2],unique_background_physical_pairs=counts[3],
                    unique_target_families=counts[4],unique_family_components=counts[5],unique_species_patterns=counts[6],unique_focal_taxa=counts[7],maximum_background_node_reuse=maxima['background'].get(group,0),maximum_background_physical_pair_reuse=maxima['pair'].get(group,0),
                    sum_inverse_background_node_reuse_weights=counts[2],sum_inverse_background_physical_pair_reuse_weights=counts[3],sum_inverse_family_component_weights=counts[5])
    db.close();database.unlink();return expected


def run(path,output):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);original=dict(bindings);root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
    assert r['status']=='complete_full_expanded_model_inputs_pending_independent_readback' and r['source_hashes']==original and r['plan_sha256']==sha(path) and r['scientific_eligibility'] is False
    bind(bindings,rp)
    for n,d in r['artifacts'].items():bind(bindings,root/n,d)
    verify(bindings);assert json.loads((root/'stage_plan.json').read_text())==dict(plan_sha256=sha(path),schema='expanded-model-input-v1')
    parts=json.loads((root/'partition_manifest.json').read_text());assert [(p['mask'],p['order_contrast']) for p in parts]==list(itertools.product(MASKS,ORDERS))
    assert set(r['artifacts'])=={'stage_plan.json','partition_manifest.json','setting_counts.tsv'}|{p['path'] for p in parts}
    valid_counts={};total=0;cells=0
    with localcontext() as ctx:
        ctx.prec=50;static=[original_static(case,source) for case in source['cases']]
        states={key:expected_state(state) for key,state in source['catalog'].items()}
        for part in parts:
            mask,order=part['mask'],part['order_contrast'];fp=root/part['path'];assert part['path']==f'inputs/{mask}-{order}.parquet' and part['sha256']==r['artifacts'][part['path']] and part['rows']==len(source['cases'])
            reader=pq.ParquetFile(fp);assert reader.schema_arrow==schema();index=valid=0
            for batch in reader.iter_batches(batch_size=256):
                for actual in batch.to_pylist():
                    case=source['cases'][index];cov=source['cov'][case['case_id']]
                    rowid=hashlib.sha256(json.dumps(['expanded-matched-measurement-row-v1',case['case_id'],mask],separators=(',',':')).encode()).hexdigest()
                    labels={k:case[k] for k in FIELDS if k in case}
                    labels.update({k:cov[k] for k in ['family_component','species_pattern_id']})
                    labels.update(input_id=identity(rowid,order),row_identity=rowid,case_row=index,species_pattern_row=int(cov['species_pattern_row']),mask=mask,order_contrast=order,
                        joint_mask_pass_bits=int(case['joint_'+mask+'_pass_bits']),joint_both_pass_bits=int(case['joint_both_pass_bits']),selection_records=int(case['selection_records']))
                    expected,scale=[dict(v) for v in static[index]]
                    selected={};chosen={'target':[0,1] if order=='mean' else [int(order[0])],'background':[0,1] if order=='mean' else [int(order[1])]}
                    for role,indices in chosen.items():
                        keys=[(role,case[role+'_pair_key'],mask,i) for i in indices]
                        aliases=case[role+'_same_model']=='1'
                        labels[role+'_state_keys']=json.dumps(['' if aliases else '|'.join(map(str,k)) for k in keys],separators=(',',':'))
                        selected[role]=[None if aliases else states[k] for k in keys]
                    usable=all(v is not None for values in selected.values() for v in values)
                    labels['numerical_usable']=int(usable);labels['disposition']='usable_order_contrast' if usable else 'identical_model_no_alignment' if case['target_same_model']=='1' or case['background_same_model']=='1' else 'unavailable_or_quarantined_order_state'
                    if usable:
                        averages={role:[sum(v[j] for v in values)/len(values) for j in range(8)] for role,values in selected.items()}
                        names=['rmsd_delta','native_tm_dissimilarity_delta','aligned_identity_delta','aligned_identity_square_delta','aligned_identity_cube_delta','original_coverage_delta','log_aligned_length_ratio','aligned_plddt70_fraction_delta']
                        for j,name in enumerate(names):expected[name]=averages['target'][j]-averages['background'][j];scale[name]=max(D(1),*(abs(v[j]) for vv in selected.values() for v in vv))
                        valid+=1
                    else:
                        for name in ['rmsd_delta','native_tm_dissimilarity_delta','aligned_identity_delta','aligned_identity_square_delta','aligned_identity_cube_delta','original_coverage_delta','log_aligned_length_ratio','aligned_plddt70_fraction_delta']:expected[name]=None
                        assert labels['joint_mask_pass_bits']==labels['joint_both_pass_bits']==0
                    assert set(labels)|set(expected)==set(FIELDS) and not set(labels)&set(expected)
                    for name,value in labels.items():assert actual[name]==value,(mask,order,index,name)
                    for name,value in expected.items():
                        if value is None:assert actual[name] is None,(mask,order,index,name)
                        else:
                            observed=D(str(actual[name]));assert observed.is_finite() and abs(observed-value)<=D('2e-14')*max(scale[name],abs(value)),(mask,order,index,name,observed,value)
                        cells+=1
                    index+=1
            assert index==len(source['cases']);total+=index;valid_counts[mask+'|'+order]=valid
            print('Independent expanded model input cells',mask,order,index,flush=True)
    expected=sqlite_counts(source,plan,root/'independent_selection_readback.sqlite');count=0;seen=set();retained=0
    with (root/'setting_counts.tsv').open() as f:
        table=csv.DictReader(f,delimiter='\t');assert table.fieldnames==COUNT_FIELDS
        for row in table:
            key=tuple(row[k] for k in ['guide','policy','scenario_id','mask','eligibility_gate','screen','order_contrast']);assert key not in seen;seen.add(key)
            assert key[3] in MASKS and key[4] in GATES and key[6] in ORDERS
            e=expected[key[:3]+(key[3] if key[4]=='mask' else 'both',key[5])]
            for field,value in e.items():assert int(row[field])==value,(key,field,row[field],value)
            count+=1;retained+=e['quality_retained_numeric_records']
    assert seen==set(itertools.product(plan['guides'],plan['policies'],[s['scenario_id'] for s in source['scenarios']],MASKS,GATES,[s['id'] for s in plan['screens']],ORDERS))
    summary=dict(logical_cases=len(source['cases']),physical_cases=plan['expected']['physical_cases'],selected_records=plan['expected']['selected_records'],unmatched_decisions=plan['expected']['unmatched_decisions'],case_mask_rows=2*len(source['cases']),model_input_rows=total,partitions=len(parts),settings=count,
        valid_input_rows=valid_counts,invalid_input_rows={k:len(source['cases'])-v for k,v in valid_counts.items()},quality_retained_setting_occurrences=retained,model_specifications_per_setting=12,future_model_setting_rows=count*12,working_tree_alternatives=5)
    assert all(r[k]==summary[k] for k in SUMMARY_FIELDS);verify(bindings)
    result=dict(status='passed_full_expanded_model_input_decimal_and_sql_readback',plan_sha256=sha(path),producer_receipt_sha256=sha(rp),**summary,numeric_cells_checked=cells,
        source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.plan,a.output)
