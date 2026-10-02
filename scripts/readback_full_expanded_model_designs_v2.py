#!/usr/bin/env python3
"""Verify full cohort/recipe identities with SQLite and independent QR/SVD ranks."""
import argparse
from collections import Counter
import csv
import gzip
import itertools
import json
from pathlib import Path
import sqlite3
import numpy as np
from scipy import linalg
from full_expanded_model_design_sources import load, array_digest, cohort_id, design_id, fit_id, AXES, OUTCOMES, DEGREES, SETTING_FIELDS, SUMMARY_FIELDS
from full_expanded_model_input_sources import ORDERS, NUISANCE, STRATUM
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def original_cohorts(source, database):
    assert not database.exists(); db=sqlite3.connect(database); db.execute('PRAGMA cache_size=-65536')
    db.execute('PRAGMA foreign_keys=ON')
    db.execute('CREATE TABLE cases (row INTEGER PRIMARY KEY,id TEXT UNIQUE,guide TEXT,target TEXT,full INTEGER,plddt70 INTEGER,both INTEGER)')
    db.executemany('INSERT INTO cases VALUES (?,?,?,?,?,?,?)',[(i,c['case_id'],c['guide'],c['target_id'],int(c['joint_full_pass_bits']),int(c['joint_plddt70_pass_bits']),int(c['joint_both_pass_bits'])) for i,c in enumerate(source['cases'])])
    db.execute('CREATE TABLE selections (ordinal INTEGER PRIMARY KEY,case_row INTEGER REFERENCES cases(row),guide TEXT,policy TEXT,scenario TEXT,target TEXT,UNIQUE(guide,policy,scenario,target))')
    lookup={c['case_id']:i for i,c in enumerate(source['cases'])};batch=[];reuse=Counter()
    with gzip.open(source['selections'],'rt') as f:
        for number,r in enumerate(csv.DictReader(f,delimiter='\t'),1):
            assert number==int(r['source_row_ordinal']);i=lookup[r['case_id']];c=source['cases'][i]
            for key in ['guide','physical_case_id','target_id','background_id']:assert r[key]==c[key]
            batch.append((number,i,r['guide'],r['policy'],r['scenario_id'],r['target_id']));reuse[r['case_id']]+=1
            if len(batch)==10000:db.executemany('INSERT INTO selections VALUES (?,?,?,?,?,?)',batch);batch=[]
    if batch:db.executemany('INSERT INTO selections VALUES (?,?,?,?,?,?)',batch)
    assert number==source['config']['expected']['selected_records'] and reuse=={c['case_id']:int(c['selection_records']) for c in source['cases']}
    db.execute('CREATE INDEX setting_selection ON selections(guide,policy,scenario)');db.commit()
    cohorts={};settings={};screen_bits={s['id']:i for i,s in enumerate(source['config']['screens'])}
    for r in source['counts']:
        key=tuple(r[k] for k in STRATUM[:-1])
        if key in settings:continue
        g,p,s,m,gate,screen=key;column=m if gate=='mask' else 'both'
        rows=np.asarray([v[0] for v in db.execute('SELECT c.row FROM selections s JOIN cases c ON c.row=s.case_row WHERE s.guide=? AND s.policy=? AND s.scenario=? AND (c.'+column+' & ?)<>0 ORDER BY c.id',(g,p,s,1<<screen_bits[screen]))],dtype=np.int64)
        assert len(rows)==int(r['quality_retained_numeric_records'])
        ids_sha=array_digest([source['cases'][i]['case_id'] for i in rows],'S64');cid=cohort_id(g,m,ids_sha)
        if cid not in cohorts:cohorts[cid]=dict(guide=g,mask=m,members=rows,ordered_case_ids_sha256=ids_sha,membership_occurrences=0)
        else:assert np.array_equal(cohorts[cid]['members'],rows)
        cohorts[cid]['membership_occurrences']+=1;settings[key]=cid
    db.close();database.unlink();return cohorts,settings


def numeric_audit(matrix, degree, saved):
    n,p=matrix.shape;assert np.isfinite(matrix).all()
    maximum=np.max(abs(matrix),axis=0,initial=0);scaled=matrix.astype(np.longdouble)/np.where(maximum>0,maximum,1)
    norms=np.sqrt(np.sum(scaled*scaled,axis=0,dtype=np.longdouble))
    active=list(range(p)) if n==0 else [i for i,v in enumerate(maximum) if v!=0]
    assert saved['records']==n and saved['columns']==p and saved['active_columns']==len(active)
    assert saved['active_column_indices']==active and saved['exactly_zero_column_indices']==[i for i in range(p) if i not in active]
    np.testing.assert_array_equal(saved['column_maxabs'],maximum)
    np.testing.assert_allclose(saved['column_l2_after_maxabs'],norms.astype(float),rtol=2e-14,atol=1e-14)
    # Independently checked saved scales permit an exact byte-level export check.
    normalized=(matrix/np.where(maximum>0,maximum,1)/np.where(np.asarray(saved['column_l2_after_maxabs'])>0,saved['column_l2_after_maxabs'],1))[:,active]
    assert saved['scaled_design_sha256']==array_digest(normalized,'<f8')
    independently_normalized=np.asarray((scaled/np.where(norms>0,norms,1))[:,active],dtype=float)
    if n:
        _,r,pivots=linalg.qr(independently_normalized,mode='economic',pivoting=True,check_finite=True)
        singular=linalg.svd(r,full_matrices=False,compute_uv=False,lapack_driver='gesvd')
    else:singular=np.empty(0)
    np.testing.assert_allclose(saved['singular_values'],singular,rtol=2e-10,atol=2e-13)
    tolerance=max(n,p)*np.finfo(float).eps*(float(singular[0]) if len(singular) else 0)
    np.testing.assert_allclose(saved['rank_tolerance'],tolerance,rtol=2e-10,atol=0)
    ranks=[int((singular>tolerance*f).sum()) for f in [.1,1,10]]
    assert [saved[k] for k in ['rank_at_one_tenth_tolerance','rank','rank_at_ten_times_tolerance']]==ranks
    q=len(active)
    status='empty_setting' if n==0 else 'sequence_axis_uninformative_requires_review' if all(i not in active for i in range(1,degree+1)) else 'insufficient_residual_dimension' if n<=q else 'rank_boundary_requires_review' if ranks[0]!=ranks[2] else 'rank_deficient_requires_review' if ranks[1]<q else 'full_rank_design'
    assert saved['disposition']==status
    condition=float(singular[0]/singular[-1]) if len(singular)==q and ranks[1]==q and len(singular) else None
    if condition is None:assert saved['normalized_condition_number'] is None
    else:
        # Raw condition numbers amplify harmless last-singular-value roundoff.
        # Verify the exported derived metric exactly, then check its stable reciprocal
        # against the independently normalized QR/SVD within a backward-error bound.
        exported = saved['normalized_condition_number']
        derived = saved['singular_values'][0] / saved['singular_values'][-1]
        np.testing.assert_allclose(exported, derived, rtol=2e-14, atol=0)
        np.testing.assert_allclose(1/exported, 1/condition, rtol=2e-10,
            atol=4*np.finfo(float).eps*max(n,p))
    return status


def run(path,output):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);original=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_full_expanded_model_designs_pending_independent_readback'
    assert receipt['plan_sha256']==sha(path) and receipt['source_hashes']==original and receipt['scientific_eligibility'] is False
    bind(bindings,rp)
    for name,d in receipt['artifacts'].items():bind(bindings,root/name,d)
    verify(bindings);assert json.loads((root/'stage_plan.json').read_text())==dict(plan_sha256=sha(path),schema='expanded-model-design-v1')
    cohorts,mapping=original_cohorts(source,root/'selection_readback.sqlite');manifest=json.loads((root/'cohort_manifest.json').read_text())
    assert [r['cohort_id'] for r in manifest]==sorted(cohorts)
    expected_artifacts={'stage_plan.json','cohort_manifest.json','unique_designs.jsonl','unique_fit_inputs.jsonl','model_settings.tsv.gz'}|{r['path'] for r in manifest}
    assert set(receipt['artifacts'])==expected_artifacts
    ds=Counter();fs=Counter();ss=Counter();model_map={};design_count=fit_count=0
    with (root/'unique_designs.jsonl').open() as df,(root/'unique_fit_inputs.jsonl').open() as ff:
        for number,entry in enumerate(manifest,1):
            cid=entry['cohort_id'];c=cohorts[cid];members=c['members'];cases=source['cases']
            assert entry['guide']==c['guide'] and entry['mask']==c['mask'] and entry['records']==len(members)
            assert entry['membership_occurrences']==c['membership_occurrences'] and entry['ordered_case_ids_sha256']==c['ordered_case_ids_sha256']
            assert entry['path']=='cohorts/'+cid+'.npz' and entry['sha256']==receipt['artifacts'][entry['path']]
            with np.load(root/entry['path']) as a:assert a.files==['case_rows'];np.testing.assert_array_equal(a['case_rows'],members)
            assert entry['case_rows_sha256']==array_digest(members,'<i8')
            dependence={
                'species_pattern_rows_sha256':array_digest([int(source['cov'][cases[i]['case_id']]['species_pattern_row']) for i in members],'<i8'),
                'family_components_sha256':array_digest([source['cov'][cases[i]['case_id']]['family_component'] for i in members],'S64'),
                'background_nodes_sha256':array_digest([cases[i]['background_id'] for i in members],'S64'),
                'background_physical_pairs_sha256':array_digest([cases[i]['background_pair_key'] for i in members],'S64'),
                'unique_family_components':len({source['cov'][cases[i]['case_id']]['family_component'] for i in members}),
                'unique_species_patterns':len({source['cov'][cases[i]['case_id']]['species_pattern_id'] for i in members}),
                'largest_family_component':max(Counter(source['cov'][cases[i]['case_id']]['family_component'] for i in members).values(),default=0)}
            for order,axis,degree in itertools.product(ORDERS,AXES,DEGREES):
                values=source['arrays'][c['mask'],order];assert np.all(values['numerical_usable'][members]==1)
                d=json.loads(next(df));did=design_id(cid,order,axis,degree,source['contract'])
                labels=dict(design_id=did,cohort_id=cid,mask=c['mask'],order_contrast=order,sequence_axis=axis,degree=degree,
                    predictor_columns=['intercept',*AXES[axis][:degree],*NUISANCE],ordered_input_ids_sha256=array_digest(values['input_id'][members],'S64'),**dependence)
                assert all(d[k]==v for k,v in labels.items())
                matrix=np.asarray([[1.]+[float(values[k][i]) for k in labels['predictor_columns'][1:]] for i in members],dtype=float).reshape(len(members),len(labels['predictor_columns']))
                assert d['raw_design_sha256']==array_digest(matrix,'<f8');status=numeric_audit(matrix,degree,d)
                assert set(d)==set(labels)|{'raw_design_sha256','records','columns','active_columns','active_column_indices','exactly_zero_column_indices','column_maxabs','column_l2_after_maxabs','singular_values','rank_tolerance','rank_at_one_tenth_tolerance','rank','rank_at_ten_times_tolerance','normalized_condition_number','scaled_design_sha256','disposition'}
                design_count+=1;ds[status]+=1
                for outcome in OUTCOMES:
                    fit=json.loads(next(ff));y=np.asarray([float(values[outcome][i]) for i in members]);response=array_digest(y,'<f8');fid=fit_id(did,outcome,response)
                    disposition=status if status!='full_rank_design' else 'constant_response_requires_review' if len(set(y))==1 else 'ready_for_working_covariance_fit'
                    expected=dict(fit_input_id=fid,design_id=did,cohort_id=cid,outcome=outcome,response_sha256=response,records=len(y),
                        response_min=min(y.tolist()) if len(y) else None,response_max=max(y.tolist()) if len(y) else None,
                        disposition=disposition,trees=plan['trees'],source_contract=source['contract'])
                    assert fit==expected;fit_count+=1;fs[disposition]+=1;model_map[cid,order,axis,degree,outcome]=(did,fid,disposition)
            if number%25==0:print('independent_expanded_design_cohorts',number,'/',len(manifest),flush=True)
        assert next(df,None) is None and next(ff,None) is None
    seen=set();settings=0
    with gzip.open(root/'model_settings.tsv.gz','rt') as f:
        r=csv.DictReader(f,delimiter='\t');assert r.fieldnames==SETTING_FIELDS
        for row in r:
            key=tuple(row[k] for k in STRATUM);full=key+(row['outcome'],row['sequence_axis'],int(row['degree']))
            assert full not in seen;seen.add(full);cid=mapping[key[:-1]]
            did,fid,status=model_map[cid,key[-1],row['sequence_axis'],int(row['degree']),row['outcome']]
            assert row['cohort_id']==cid and row['design_id']==did and row['fit_input_id']==fid and row['disposition']==status
            assert int(row['records'])==len(cohorts[cid]['members']) and int(row['nominal_tree_fits'])==len(plan['trees'])
            settings+=1;ss[status]+=1
    expected_keys={tuple(r[k] for k in STRATUM)+(o,a,d) for r in source['counts'] for o,a,d in itertools.product(OUTCOMES,AXES,DEGREES)}
    assert seen==expected_keys and settings==plan['expected']['model_setting_rows']
    summary=dict(logical_cases=len(source['cases']),selected_records=source['config']['expected']['selected_records'],input_setting_rows=len(source['counts']),model_setting_rows=settings,
        unique_cohorts=len(cohorts),unique_designs=design_count,unique_fit_inputs=fit_count,nominal_tree_setting_fits=settings*len(plan['trees']),unique_tree_fit_inputs=fit_count*len(plan['trees']),
        design_status_counts=dict(ds),fit_input_status_counts=dict(fs),setting_status_counts=dict(ss),largest_cohort=max(len(c['members']) for c in cohorts.values()),
        cohort_member_occurrences=sum(len(c['members']) for c in cohorts.values()),trees=plan['trees'])
    assert all(receipt[k]==summary[k] for k in SUMMARY_FIELDS);verify(bindings)
    result=dict(status='passed_full_expanded_model_design_sql_qr_readback',plan_sha256=sha(path),producer_receipt_sha256=sha(rp),**summary,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.plan,a.output)
