#!/usr/bin/env python3
"""Verify all source recipes, latent error-contrast products and setting links."""
from collections import Counter
import argparse
import csv
import gzip
import itertools
import json
from pathlib import Path
import sqlite3
import tempfile
import numpy as np
from scipy import linalg
from covariance_basis_audit import _diagnostics
from covariance_basis_context import ComponentKernelProducts
from covariance_basis_independent import IndependentKernelProducts
from full_covariance_qualification_sources import load,jsonl,cohort_rows,design_matrix,folded_operators,audit_id,covariance_status,SUMMARY_FIELDS,LINK_EXTRA
from full_entity_operator_sources import MODES
from full_expanded_model_input_sources import ORDERS
from full_expanded_model_design_sources import AXES,DEGREES,SETTING_FIELDS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def numeric(saved,reference,names,records):
    assert saved['kernel_names']==names and saved['records']==records
    assert saved['covariance_model_selected'] is False
    condition=saved['normalized_design_condition_number'];expected=reference['normalized_design_condition_number']
    assert condition>=1 and np.isfinite(condition)
    assert abs(1/condition-1/expected)<=4*np.finfo(float).eps*max(records,saved['fixed_effect_columns'])
    for field,key in [('raw_gram','raw'),('projected_gram','projected')]:
        actual=np.asarray(saved[field]);target=reference[key]
        assert actual.shape==target.shape==(len(names),len(names)) and np.isfinite(actual).all()
        np.testing.assert_allclose(actual,target,rtol=3e-9,atol=2e-8)
    np.testing.assert_allclose(saved['raw_roundoff_envelope'],reference['raw_error'],rtol=2e-9,atol=1e-20)
    # The reciprocal-condition bound above verifies ill-conditioned cases
    # without an unjustified fixed relative tolerance on condition numbers.
    expected_error=reference['projected_error']*(condition/expected)
    np.testing.assert_allclose(saved['projected_roundoff_envelope'],expected_error,rtol=3e-8,atol=1e-18)
    independently_qualified=True;conservative_review=False
    for field,key,error in [('raw_diagnostics','raw','raw_roundoff_envelope'),
                            ('reml_diagnostics','projected','projected_roundoff_envelope')]:
        recorded=_diagnostics(np.asarray(saved[key+'_gram']),np.asarray(saved[error]),names)
        assert saved[field]==recorded
        independent=_diagnostics(reference[key],np.asarray(saved[error]),names)
        active=independent['resolved_kernel_indices']
        if active:
            gram=reference[key][np.ix_(active,active)];d=np.sqrt(np.diag(gram))
            normalized=gram/d[:,None]/d[None,:]
            singular=linalg.svd(normalized,compute_uv=False,lapack_driver='gesvd')
            ranks=[int(np.count_nonzero(singular>independent['rank_tolerance']*v)) for v in [.1,1.,10.]]
            assert ranks==[independent[k] for k in ['rank_at_one_tenth_tolerance','rank','rank_at_ten_times_tolerance']]
        independently_qualified &= independent['disposition']=='numerically_independent_covariance_bases'
        conservative_review |= independent['disposition']!=saved[field]['disposition']
    if covariance_status(saved)=='qualified_uniform_working_covariance_basis':
        assert independently_qualified,'Producer qualified a basis requiring independent review'
    return conservative_review,float(np.max(abs(np.asarray(saved['projected_gram'])-reference['projected']),initial=0))


def run(path,output):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);original=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_full_uniform_covariance_qualification_pending_independent_readback'
    assert receipt['plan_sha256']==sha(path) and receipt['source_contract']==source['contract']
    assert receipt['source_hashes']==original and receipt['scientific_eligibility'] is False
    bind(bindings,rp)
    assert set(receipt['artifacts'])=={'stage_plan.json','design_covariance_audits.jsonl.gz','setting_audit_links.tsv.gz'}
    for name,d in receipt['artifacts'].items():bind(bindings,root/name,d)
    verify(bindings)
    assert json.loads((root/'stage_plan.json').read_text())==dict(plan_sha256=sha(path),schema='expanded-uniform-covariance-qualification-v1')
    counts=Counter();linkcounts=Counter();design_count=audits=cohort_count=reviews=0;maximum_error=0.
    designs=jsonl(source['root']/'unique_designs.jsonl');exported=jsonl(root/'design_covariance_audits.jsonl.gz')
    with tempfile.TemporaryDirectory(prefix='uniform-covariance-readback-',dir=root) as directory:
        db=sqlite3.connect(str(Path(directory)/'audits.sqlite'))
        db.execute('CREATE TABLE audits(aid TEXT PRIMARY KEY,did TEXT,mode TEXT,tree TEXT,status TEXT,UNIQUE(did,mode,tree))')
        for cohort in source['cohorts']:
            rows=cohort_rows(source,cohort);contexts={};fallback={};seen=set()
            if len(rows):
                factors={tree:source['factors'][tree][rows] for tree in plan['trees']}
                for mode in MODES:
                    operators=folded_operators(source,rows,mode)
                    independent=IndependentKernelProducts(source['labels'][rows],operators)
                    contexts[mode]={tree:independent.tree(factors[tree]) for tree in plan['trees']}
                    fallback[mode]=(operators,factors)
            for _ in range(30):
                design=next(designs);matrix=design_matrix(source,cohort,rows,design);design_count+=1
                seen.add((design['order_contrast'],design['sequence_axis'],design['degree']))
                for mode,tree in itertools.product(MODES,plan['trees']):
                    record=next(exported);identifier=audit_id(source['contract'],design['design_id'],mode,tree)
                    expected=dict(audit_id=identifier,source_contract=source['contract'],cohort_id=cohort['cohort_id'],
                        design_id=design['design_id'],loading_mode=mode,tree=tree,residual_diagonal='uniform_one',
                        records=len(rows),source_design_disposition=design['disposition'],
                        active_column_indices=design['active_column_indices'],raw_design_sha256=design['raw_design_sha256'],
                        folded_terms={'target_node':'combined_with_uniform_residual','family':'zero_signed_endpoint_kernel' if mode=='signed' else
                            'combined_variance_family_intercept_plus_four_times_endpoint_family'})
                    assert all(record[k]==v for k,v in expected.items())
                    if design['disposition']!='full_rank_design':
                        assert record['disposition']==design['disposition'] and record['numerical_audit'] is None
                    elif record['disposition']=='numerical_covariance_qualification_requires_review':
                        # Reproduce the original algorithm's explicit failure; an independent
                        # solve may behave differently, but no case is promoted to qualified.
                        operators,factors=fallback[mode]
                        try:ComponentKernelProducts(source['labels'][rows],operators,np.ones(len(rows))).tree(factors[tree]).audit(matrix)
                        except (ValueError,ArithmeticError) as error:
                            assert record['error_type']==type(error).__name__ and record['error_message']==str(error)
                        else:raise AssertionError('False numerical-failure disposition')
                        assert record['numerical_audit'] is None
                    else:
                        reference=contexts[mode][tree].project(matrix);value=record['numerical_audit']
                        assert value['fixed_effect_columns']==matrix.shape[1] and value['residual_dimension']==len(rows)-matrix.shape[1]
                        assert value['family_components']==len(contexts[mode][tree].bank.parts)
                        assert value['exactly_zero_incidence_names']==[k for k,z in fallback[mode][0].items() if not z.nnz]
                        conservative,error=numeric(value,reference,contexts[mode][tree].bank.names,len(rows))
                        reviews+=int(conservative);maximum_error=max(maximum_error,error)
                        assert record['disposition']==covariance_status(value)
                    audits+=1;counts[record['disposition']]+=1
                    db.execute('INSERT INTO audits VALUES(?,?,?,?,?)',(identifier,design['design_id'],mode,tree,record['disposition']))
            assert seen==set(itertools.product(ORDERS,AXES,DEGREES));cohort_count+=1;db.commit()
            print('independent_full_covariance_cohort',cohort_count,'/',len(source['cohorts']),flush=True)
        assert next(designs,None) is None and next(exported,None) is None
        settings=links=0
        with gzip.open(source['root']/'model_settings.tsv.gz','rt') as f,gzip.open(root/'setting_audit_links.tsv.gz','rt') as g:
            original_settings=csv.DictReader(f,delimiter='\t');exported_links=csv.DictReader(g,delimiter='\t')
            assert original_settings.fieldnames==SETTING_FIELDS and exported_links.fieldnames==SETTING_FIELDS+LINK_EXTRA
            for row in original_settings:
                for mode,tree in itertools.product(MODES,plan['trees']):
                    record=next(exported_links);assert {k:record[k] for k in SETTING_FIELDS}==row
                    identifier=audit_id(source['contract'],row['design_id'],mode,tree)
                    selected=db.execute('SELECT did,mode,tree,status FROM audits WHERE aid=?',(identifier,)).fetchone()
                    assert selected is not None and selected[:3]==(row['design_id'],mode,tree)
                    status=selected[3];combined=row['disposition'] if row['disposition']!='ready_for_working_covariance_fit' else status
                    assert [record[k] for k in LINK_EXTRA]==[mode,tree,identifier,status,combined]
                    links+=1;linkcounts[combined]+=1
                settings+=1
            assert next(exported_links,None) is None
        assert db.execute('SELECT COUNT(*) FROM audits').fetchone()[0]==audits;db.close()
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=settings,unique_cohorts=cohort_count,unique_designs=design_count,
        audit_rows=audits,setting_audit_links=links,audit_status_counts=dict(counts),link_status_counts=dict(linkcounts),trees=plan['trees'],loading_modes=MODES)
    assert all(receipt[k]==summary[k] for k in SUMMARY_FIELDS)
    assert settings==source['design_completion']['model_setting_rows'] and design_count==source['design_completion']['unique_designs']
    assert audits==design_count*len(MODES)*len(plan['trees']) and links==settings*len(MODES)*len(plan['trees'])
    verify(bindings)
    result=dict(status='passed_full_uniform_covariance_latent_and_sql_readback',plan_sha256=sha(path),producer_receipt_sha256=sha(rp),
        source_contract=source['contract'],**summary,maximum_absolute_projected_gram_error=maximum_error,
        conservative_independent_classification_differences=reviews,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.output)
