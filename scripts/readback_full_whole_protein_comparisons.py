#!/usr/bin/env python3
"""Reconstruct every selected source/parameter independently and check all comparisons."""
import argparse
from collections import Counter
import gzip
from itertools import zip_longest
import json
from pathlib import Path
import numpy as np
from full_whole_protein_comparison_sources import load_sources,SUMMARY_FIELDS
from readback_whole_protein_ml_comparisons import verify_comparison
from run_whole_protein_ml import digest
from reference_measurement_union_sources import bind,verify
from screen_duplication_alignment_reuse import sha


def exact(a,b):
    # JSON distinguishes true from 1 and rejects NaN; Python dictionary equality does not.
    assert json.dumps(a,sort_keys=True,allow_nan=False)==json.dumps(b,sort_keys=True,allow_nan=False)


def reconstruct(entry,saved,follow,follow_root,support):
    identifier,tree=entry['fit_input_id'],entry['tree'];original=saved['payload'];spec=saved['specification'];state=original['status']
    assert saved['fit_input_id']==identifier and saved['tree']==tree and state==entry['status']
    assert saved['payload_sha256']==digest(original) and identifier==digest(spec)
    assert type(spec['records']) is int and spec['records']>len(spec['columns'])-1
    assert len(set(spec['columns']))==len(spec['columns']) and spec['columns'][1]=='intercept'
    source=entry['path'];source_hash=entry['sha256'];kind='original_production';effective=state
    candidate=original;raw=original;theta=original.get('log1p_ratios');case_path=case_hash=case_status=None
    if state=='ml_candidate_requires_optimization_review':
        r=follow[identifier,tree];case_path=str(Path(follow_root)/r['path']);case_hash=r['sha256']
        assert sha(case_path)==case_hash
        case=json.loads(Path(case_path).read_text());identity=case['identity'];result=case['result'];case_status=result['status']
        assert identity['fit_input_id']==identifier and identity['tree']==tree
        assert identity['original_fit']==entry['path'] and identity['original_fit_sha256']==entry['sha256'] and identity['input_sha256']==saved['input_sha256']
        assert case['result_sha256']==digest(result) and case_status==r['status'] and case['original_status']==state and case['scientific_eligibility'] is False
        if case_status=='numerical_followup_error_requires_review':
            kind='original_retained_after_followup_error'
        elif case_status in ['numerical_followup_passed_pending_independent_readback','numerical_followup_requires_review']:
            candidate=result['fitted'];raw=result;theta=result['selected_theta'];source=case_path;source_hash=case_hash
            kind='full_grid_refinement' if result['recovery'] is None else 'full_grid_recovery'
            effective={'numerical_followup_passed_pending_independent_readback':'ml_refinement_passed_numerical_checks',
                'numerical_followup_requires_review':'ml_full_grid_followup_requires_review'}[case_status]
        else:raise AssertionError('Unknown closed follow-up status')
    else:
        assert state in ['ml_candidate_passed_numerical_optimization_checks','fit_error_requires_review'] and (identifier,tree) not in follow
    missing=state=='fit_error_requires_review';width=len(spec['columns'])-1
    if not missing:
        beta=np.asarray(candidate['beta'],float);cov=np.asarray(candidate['conditional_beta_covariance'],float);scales=np.asarray(original['covariate_scales'],float)
        assert beta.shape==(width,) and cov.shape==(width,width) and scales.shape==(width-1,)
        assert np.isfinite(beta).all() and np.isfinite(cov).all() and np.isfinite(scales).all() and np.all(scales>0)
        conversion=np.concatenate(([1.],1/scales))
        np.testing.assert_allclose(raw['raw_unit_beta'],[beta[i]*conversion[i] for i in range(width)],rtol=1e-9,atol=1e-10)
        np.testing.assert_allclose(raw['raw_unit_conditional_beta_covariance'],[[cov[i,j]*conversion[i]*conversion[j] for j in range(width)] for i in range(width)],rtol=1e-9,atol=1e-10)
        assert np.isfinite(candidate['negative_profiled_ml']) and np.isfinite(candidate['profiled_scale']) and candidate['profiled_scale']>0
        assert type(candidate['variance_profile_denominator']) is int and candidate['variance_profile_denominator']==spec['records']
        assert len(theta)==3 and np.isfinite(theta).all() and np.all(np.asarray(theta)>=0)
    supported=support.resolve(identifier,saved['input_sha256'])
    s=dict(fit_input_id=identifier,tree=tree,input_sha256=saved['input_sha256'],original_status=state,effective_status=effective,
        original_negative_profiled_ml=original.get('negative_profiled_ml'),negative_profiled_ml=None if missing else candidate['negative_profiled_ml'],records=spec['records'],fixed_coefficients=width,
        original_production_fit=entry['path'],original_production_sha256=entry['sha256'],selected_source_kind=kind,selected_source_path=source,selected_source_sha256=source_hash,
        followup_status=case_status,followup_source_path=case_path,followup_source_sha256=case_hash,numerical_fit_verified=not missing,
        support_geometry_id=supported['geometry_id'],support_classification=supported['classification'],original_support_classification=supported['original_classification'],support_source_kind=supported['source_kind'])
    p=dict(fit_input_id=identifier,tree=tree,columns=spec['columns'],specification=spec,selected_source_path=source,selected_source_sha256=source_hash,effective_status=effective,numerical_fit_verified=not missing,
        covariate_scales=None if missing else original['covariate_scales'],scaled_beta=None if missing else candidate['beta'],
        scaled_conditional_beta_covariance=None if missing else candidate['conditional_beta_covariance'],raw_unit_beta=None if missing else raw['raw_unit_beta'],
        raw_unit_conditional_beta_covariance=None if missing else raw['raw_unit_conditional_beta_covariance'],selected_log1p_ratios=None if missing else theta,
        profiled_scale=None if missing else candidate['profiled_scale'],variance_profile_denominator=None if missing else candidate['variance_profile_denominator'],scientific_eligibility=False)
    return s,p


def run(path,output):
    plan=json.loads(Path(path).read_text());original,follow,follow_root,lr,support,bindings=load_sources(plan,path)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text());bind(bindings,rp)
    assert receipt['status']=='complete_full_whole_protein_comparisons_pending_independent_readback' and receipt['scientific_eligibility'] is False and receipt['plan_sha256']==sha(path)
    for p,h in receipt['source_hashes'].items():bind(bindings,p,h)
    for n,h in receipt['artifacts'].items():bind(bindings,root/n,h)
    verify(bindings)
    summaries={};old=Counter();states=Counter();cases=Counter();kinds=Counter()
    with (root/'fit_summaries.jsonl').open() as sf,gzip.open(root/'selected_fit_parameters.jsonl.gz','rt') as pf:
        for entry,sl,pl in zip_longest(original.values(),sf,pf):
            assert entry is not None and sl is not None and pl is not None
            saved=json.loads(Path(entry['path']).read_text());assert saved['plan_sha256']==sha(plan['fit_plan']) and sha(entry['path'])==entry['sha256']
            s,p=reconstruct(entry,saved,follow,follow_root,support);exact(json.loads(sl),s);exact(json.loads(pl),p)
            key=s['fit_input_id'],s['tree'];assert key not in summaries;summaries[key]=s
            old[s['original_status']]+=1;states[s['effective_status']]+=1;kinds[s['selected_source_kind']]+=1
            if s['followup_status'] is not None:cases[s['followup_status']]+=1
            if len(summaries)%25000==0:print('Read back selected full-grid sources and parameters',len(summaries),flush=True)
    trees=sorted({k[1] for k in summaries});ids={k[0] for k in summaries}
    assert len(summaries)==plan['expected']['fit_summaries'] and len(ids)==plan['expected']['unique_inputs'] and len(trees)==plan['expected']['trees']
    assert all((i,t) in summaries for i in ids for t in trees) and sum(cases.values())==len(follow)
    total=0;counts=Counter();relations=Counter();used=set();expected_artifacts={'fit_summaries.jsonl','selected_fit_parameters.jsonl.gz'}
    configuration=dict(plan_sha256=sha(path),summary_sha256=sha(root/'fit_summaries.jsonl'),parameters_sha256=sha(root/'selected_fit_parameters.jsonl.gz'))
    for tree in trees:
        n=0;local=Counter();rel=Counter();target=root/(tree+'.jsonl.gz');checkpoint=root/(tree+'.checkpoint.json');expected_artifacts.update([target.name,checkpoint.name])
        with gzip.open(target,'rt') as exported,(Path(plan['links'])/'comparison_input_map.jsonl').open() as links:
            for source_line,output_line in zip_longest(links,exported):
                assert source_line is not None and output_line is not None
                link=json.loads(source_line);row=json.loads(output_line);exact({k:row[k] for k in link},link);assert row['tree']==tree
                assert set(row)==set(link)|{'tree','left','right','log_likelihood_gain_right_vs_left','twice_log_likelihood_gain_right_vs_left','nested_model_log_likelihood_gain','nested_added_fixed_coefficients','numerical_likelihood_tolerance','nesting_status','comparison_status'}
                keys=[(link[k],tree) for k in ['left_input','right_input']];a,b=[summaries[k] for k in keys];used.update(keys)
                exact(row['left'],a);exact(row['right'],b);verify_comparison(row,a,b,link['relation'])
                for k in ['log_likelihood_gain_right_vs_left','twice_log_likelihood_gain_right_vs_left','nested_model_log_likelihood_gain','numerical_likelihood_tolerance']:
                    assert row[k] is None or type(row[k]) in [int,float] and np.isfinite(row[k])
                assert row['nested_added_fixed_coefficients'] is None or type(row['nested_added_fixed_coefficients']) is int
                n+=1;local[row['comparison_status']]+=1;rel[link['relation']]+=1
        c=json.loads(checkpoint.read_text());exact(c,dict(configuration=configuration,sha256=sha(target),comparisons=n,status_counts=dict(local),relation_counts=dict(rel)))
        assert n==plan['expected']['comparisons_per_tree'];exact(dict(rel),lr['relation_counts']);total+=n;counts.update(local);relations.update(rel)
        print('Read back complete optimized comparison tree',tree,total,flush=True)
    assert total==plan['expected']['comparisons'] and used==set(summaries) and set(receipt['artifacts'])==expected_artifacts
    observed=dict(fit_summaries=len(summaries),parameter_records=len(summaries),comparisons=total,trees=trees,original_fit_status_counts=dict(old),effective_fit_status_counts=dict(states),followup_status_counts=dict(cases),selected_source_kind_counts=dict(kinds),comparison_status_counts=dict(counts),relation_counts=dict(relations),followup_cases_integrated=len(follow))
    for k in SUMMARY_FIELDS:exact(receipt[k],observed[k])
    verify(bindings)
    result=dict(status='passed_full_whole_protein_comparisons_source_selection_arithmetic_readback',**observed,producer_receipt_sha256=sha(rp),plan_sha256=sha(path),source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    Path(output).parent.mkdir(exist_ok=True,parents=True)
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run(a.plan,a.output)
    print(json.dumps({k:v for k,v in r.items() if k!='source_hashes'}),flush=True)
