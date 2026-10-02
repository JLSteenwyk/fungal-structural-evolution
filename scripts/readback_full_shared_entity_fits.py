#!/usr/bin/env python3
"""Full original-fit identity replay and independent spectral numerical audits."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import json
from pathlib import Path
import sqlite3
import numpy as np
from full_shared_entity_fit_sources import load,cohorts,cases,candidate_id,FIT_LINK_EXTRA,SUMMARY_FIELDS
from full_covariance_qualification_sources import folded_operators,LINK_EXTRA
from full_expanded_model_design_sources import SETTING_FIELDS,digest
from prepare_full_shared_entity_fits import candidate as original_candidate
from independent_shared_entity_likelihood import IndependentSharedEntityLikelihood,audit_candidate,audit_local_curvature
from independent_shared_entity_optimizer import compare_independent_optimizer
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def numeric(source,plan,rows,identity,matrix,response,exported,operators=None,source_audit=None):
    assert all(exported[k]==v for k,v in identity.items())
    expected_keys=set(identity)|{'disposition','fit','numerical_attempted'}
    if exported['disposition']=='shared_entity_fit_error_requires_review':expected_keys|={'error_type','error_message'}
    assert set(exported)==expected_keys
    fit=exported['fit'];source_ready=identity['source_combined_disposition']=='qualified_uniform_working_covariance_basis'
    if not source_ready:
        assert exported==dict(**identity,disposition=identity['source_combined_disposition'],fit=None,numerical_attempted=False)
        return dict(disposition='source_review_or_exclusion_retained',scientific_eligibility=False)
    assert exported['numerical_attempted'] is True
    if fit is None or 'variance_ratios' not in fit:
        # Preserve real failures without permitting an invented failure to hide
        # an original source row. Deterministic production replay verifies the
        # claimed error; it never independently promotes the failed candidate.
        replay=original_candidate(source,plan,rows,identity,matrix,response,operators)
        assert replay==exported
        return dict(disposition='original_numerical_failure_reproduced_requires_review',scientific_eligibility=False)
    assert exported['disposition']==fit['status'] and fit['method']==identity['method']
    if any(start['disposition']=='failed_search_requires_review' for start in fit['starts']):
        assert original_candidate(source,plan,rows,identity,matrix,response,operators)==exported
    if source_audit is not None:assert fit['source_covariance_qualification']==source_audit['numerical_audit']
    assert len(fit['starts'])==3 and fit['selected_start'] in range(3)
    if operators is None:operators=folded_operators(source,rows,identity['loading_mode'])
    independent=IndependentSharedEntityLikelihood(source['labels'][rows],operators,source['factors'][identity['tree']][rows],
        np.ones(len(rows)),matrix,response,column_batch=plan['independent_audit']['column_batch'])
    replay=audit_candidate(independent,fit,**plan['independent_audit']['replay'])
    assert replay['objective_absolute_error']<=plan['independent_audit']['replay']['objective_atol']
    assert replay['coefficient_and_conditional_covariance_comparison_passed'] and replay['profiled_scale_comparison_passed'] and replay['variance_components_comparison_passed']
    norms=np.asarray(fit['kernel_normalization']);upper=np.log1p(fit['maximum_scaled_variance'])
    assert fit['maximum_scaled_variance']==plan['optimizer']['maximum_scaled_variance']
    tolerance=64*np.finfo(float).eps*max(1.,upper);finite=[];start_checks=[]
    middle=min(1.,fit['maximum_scaled_variance']/2);rng=np.random.default_rng(20261002)
    initial_points=[np.zeros(len(norms)),np.full(len(norms),np.log1p(middle)),np.log1p(middle*rng.uniform(.1,1.8,len(norms)))]
    for number,start in enumerate(fit['starts']):
        assert start['start']==number and 0<start['search_evaluations']<=plan['optimizer']['max_evaluations']
        np.testing.assert_allclose(start['initial_coordinates'],initial_points[number],rtol=1e-13,atol=1e-14)
        if 'final_coordinates' not in start:
            assert start['disposition']=='failed_search_requires_review'
            assert start['error_type'] in ['ValueError','ArithmeticError','RuntimeError','LinAlgError']
            continue
        assert start['final_replay_evaluations']==1 and 0<=start['iterations']<=plan['optimizer']['max_iterations']
        point=np.asarray(start['final_coordinates'],dtype=float)
        assert point.shape==norms.shape and np.isfinite(point).all() and np.all((point>=0)&(point<=upper))
        value=independent.evaluate(np.expm1(point)/norms,fit['method'])
        error=abs(value['negative_profiled_likelihood']-start['objective'])
        assert error<=plan['independent_audit']['replay']['objective_atol']
        gradient=value['gradient']*np.exp(point)/norms;low=point<=tolerance;high=point>=upper-tolerance
        gradient[low]=np.minimum(gradient[low],0);gradient[high]=np.maximum(gradient[high],0)
        kkt=float(np.max(abs(gradient),initial=0))
        np.testing.assert_allclose(kkt,start['maximum_projected_gradient'],rtol=1e-3,atol=plan['independent_audit']['replay']['gradient_atol'])
        assert start['lower_boundary_indices']==np.flatnonzero(low).tolist() and start['upper_boundary_indices']==np.flatnonzero(high).tolist()
        expected='candidate_pending_independent_audit' if start['optimizer_success'] and start['maximum_projected_gradient']<=plan['optimizer']['gradient_tolerance'] and not np.any(high) else 'optimizer_or_variance_boundary_requires_review'
        assert start['disposition']==expected
        finite.append((start['objective'],number));start_checks.append(dict(start=number,objective_absolute_error=float(error),independent_projected_gradient=kkt))
    assert finite and min(finite)[1]==fit['selected_start']
    selected=fit['starts'][fit['selected_start']]
    np.testing.assert_allclose(selected['final_coordinates'],np.log1p(np.asarray(fit['variance_ratios'])*norms),rtol=1e-13,atol=1e-14)
    objectives=[f[0] for f in finite];spread=max(objectives)-min(objectives)
    np.testing.assert_allclose(fit['objective_spread'],spread,rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(fit['negative_profiled_likelihood'],min(objectives),rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(fit['multistart_objective_tolerance'],1e-7*max(1.,abs(min(objectives))),rtol=1e-13,atol=0)
    expected_status='optimized_shared_entity_candidate_pending_independent_audit' if all(s['disposition']=='candidate_pending_independent_audit' for s in fit['starts']) and spread<=fit['multistart_objective_tolerance'] else 'shared_entity_optimizer_candidate_requires_review'
    assert fit['status']==expected_status
    result=dict(disposition='producer_candidate_requires_review_independent_numeric_replay_passed',candidate_replay=replay,
        original_start_replays=start_checks,scientific_eligibility=False)
    if fit['status']=='optimized_shared_entity_candidate_pending_independent_audit':
        try:
            curvature=audit_local_curvature(independent,fit,**plan['independent_audit']['curvature'])
            search=compare_independent_optimizer(independent,fit,**plan['independent_audit']['optimizer'])
            result.update(local_curvature=curvature,independent_search=search)
            if curvature['status']=='independent_shared_entity_local_minimum_check_passed_pending_global_and_inferential_audits' and search['status']=='independent_multistart_searches_agree_pending_inferential_calibration':
                result['disposition']='independently_audited_working_candidate_pending_inferential_calibration'
            else:result['disposition']='independent_curvature_or_search_requires_review'
        except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
            result.update(disposition='independent_numerical_precision_requires_review',error_type=type(error).__name__,error_message=str(error))
    return result


def run(path,output):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);original=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    lock=(root/'readback.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    completed_marker=root/'independent_readback_completed.json'
    if completed_marker.exists():
        previous=json.loads(completed_marker.read_text());assert sha(previous['output'])==previous['output_sha256']
        raise AssertionError('Completed independent fitting readback cannot rerun under another output name')
    assert not Path(output).exists(),'Completed independent fitting readback cannot restart'
    assert receipt['status']=='complete_full_shared_entity_fits_pending_independent_readback'
    assert receipt['plan_sha256']==sha(path) and receipt['source_contract']==source['fit_contract']
    assert receipt['source_hashes']==original and receipt['scientific_eligibility'] is False
    bind(bindings,rp)
    for name,d in receipt['artifacts'].items():bind(bindings,root/name,d)
    verify(bindings)
    stage=dict(plan_sha256=sha(path),source_contract=source['fit_contract'],schema='full-shared-entity-uniform-fit-v1')
    assert json.loads((root/'stage_plan.json').read_text())==stage
    readback_state=dict(schema='full-shared-entity-spectral-readback-v1',plan_sha256=sha(path),producer_receipt_sha256=sha(rp))
    marker=root/'readback_stage.json'
    if marker.exists():assert json.loads(marker.read_text())==readback_state
    else:marker.write_text(json.dumps(readback_state,indent=2)+'\n')
    manifest=json.loads((root/'cohort_manifest.json').read_text())
    assert [r['cohort_id'] for r in manifest]==[r['cohort_id'] for r in source['cohorts']]
    assert len({r['cohort_id'] for r in manifest})==len(manifest)
    expected_artifacts={'stage_plan.json','cohort_manifest.json','setting_fit_links.tsv.gz'}|{r[k] for r in manifest for k in ['path','receipt_path']}
    assert set(receipt['artifacts'])==expected_artifacts
    database=root/'fit_readback.sqlite'
    if database.exists():database.unlink()  # Derived scratch; lock/state guard a full independent replay.
    db=sqlite3.connect(database);db.execute('PRAGMA cache_size=-65536')
    db.execute('CREATE TABLE candidates(id TEXT PRIMARY KEY,fit TEXT,mode TEXT,tree TEXT,method TEXT,audit TEXT,status TEXT,independent TEXT,UNIQUE(fit,mode,tree,method))')
    counts=Counter();audited=Counter();total=0;design_count=0
    audits_path=root/'independent_candidate_audits.jsonl.gz'
    with gzip.open(audits_path,'wt') as audit_output:
        for number,((cohort,rows,entries),part) in enumerate(zip(cohorts(source,plan),manifest),1):
            cid=cohort['cohort_id'];assert part['path']=='cohorts/'+cid+'.jsonl.gz' and part['receipt_path']=='cohorts/'+cid+'.receipt.json'
            saved=json.loads((root/part['receipt_path']).read_text());assert saved['stage']==stage and saved['cohort_id']==cid
            assert saved['candidate_sha256']==part['sha256']==receipt['artifacts'][part['path']]
            assert part['receipt_sha256']==receipt['artifacts'][part['receipt_path']]
            operators={mode:folded_operators(source,rows,mode) for mode in plan['loading_modes']} if len(rows) else {}
            source_audits={(d['design_id'],mode,tree):value for d,_,_,qualified in entries for (mode,tree),value in qualified.items()}
            local=Counter();seen=0
            with gzip.open(root/part['path'],'rt') as f:
                exported=(json.loads(line) for line in f)
                for expected,matrix,response in cases(source,plan,cohort,entries):
                    value=next(exported);assert value['scientific_eligibility'] is False
                    audit=numeric(source,plan,rows,expected,matrix,response,value,operators.get(expected['loading_mode']),
                        source_audits[expected['design_id'],expected['loading_mode'],expected['tree']])
                    record=dict(candidate_id=expected['candidate_id'],source_candidate_sha256=digest(value),**audit)
                    audit_output.write(json.dumps(record,sort_keys=True,allow_nan=False)+'\n')
                    db.execute('INSERT INTO candidates VALUES (?,?,?,?,?,?,?,?)',(value['candidate_id'],value['fit_input_id'],value['loading_mode'],value['tree'],value['method'],value['covariance_audit_id'],value['disposition'],audit['disposition']))
                    local[value['disposition']]+=1;audited[audit['disposition']]+=1;total+=1;seen+=1
                assert next(exported,None) is None
            assert part['candidate_rows']==saved['candidate_rows']==seen and saved['status_counts']==dict(local)
            counts.update(local);design_count+=len(entries);db.commit()
            print('independent_full_shared_entity_cohort',number,'/',len(manifest),'candidates',seen,flush=True)
    settings=Counter();links=0;link_audits=Counter()
    independent_links=root/'independent_setting_fit_links.tsv.gz'
    with gzip.open(source['qualification_root']/'setting_audit_links.tsv.gz','rt') as f,gzip.open(root/'setting_fit_links.tsv.gz','rt') as g,gzip.open(independent_links,'wt') as h:
        original_rows=csv.DictReader(f,delimiter='\t');reader=csv.DictReader(g,delimiter='\t')
        assert original_rows.fieldnames==SETTING_FIELDS+LINK_EXTRA and reader.fieldnames==SETTING_FIELDS+LINK_EXTRA+FIT_LINK_EXTRA
        writer=csv.DictWriter(h,reader.fieldnames+['independent_disposition'],delimiter='\t',lineterminator='\n');writer.writeheader()
        for row in original_rows:
            for method in plan['methods']:
                exported=next(reader);cid=candidate_id(source['fit_contract'],row['fit_input_id'],row['loading_mode'],row['tree'],method)
                saved=db.execute('SELECT audit,status,independent FROM candidates WHERE id=?',(cid,)).fetchone()
                assert saved is not None and saved[0]==row['audit_id']
                assert exported=={**row,'likelihood_method':method,'candidate_id':cid,'candidate_disposition':saved[1]}
                writer.writerow({**exported,'independent_disposition':saved[2]})
                settings[saved[1]]+=1;link_audits[saved[2]]+=1;links+=1
        assert next(reader,None) is None
    db.close();database.unlink()
    q=source['qualification_completion']
    assert total==design_count*2*len(plan['methods'])*len(plan['loading_modes'])*len(plan['trees'])
    assert design_count==q['unique_designs'] and links==q['setting_audit_links']*len(plan['methods'])
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=q['model_setting_rows'],unique_cohorts=len(manifest),
        unique_designs=design_count,unique_fit_inputs=design_count*2,candidate_rows=total,setting_fit_links=links,
        candidate_status_counts=dict(counts),setting_status_counts=dict(settings),methods=plan['methods'],trees=plan['trees'],loading_modes=plan['loading_modes'])
    assert all(receipt[k]==summary[k] for k in SUMMARY_FIELDS)
    for p in [audits_path,independent_links,marker]:bind(bindings,p)
    verify(bindings)
    result=dict(status='passed_full_shared_entity_original_grid_and_spectral_readback',plan_sha256=sha(path),producer_receipt_sha256=sha(rp),
        source_contract=source['fit_contract'],**summary,independent_candidate_status_counts=dict(audited),
        independent_setting_status_counts=dict(link_audits),source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    with completed_marker.open('x') as f:json.dump(dict(output=str(output),output_sha256=sha(output),
        plan_sha256=sha(path),producer_receipt_sha256=sha(rp)),f,indent=2);f.write('\n')
    print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.plan,a.output)
