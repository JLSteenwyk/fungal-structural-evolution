#!/usr/bin/env python3
"""Check full timing census/restart contracts and actual numerical probes."""
import argparse
import copy
import gzip
import itertools
import json
from pathlib import Path
import tempfile
import numpy as np
from scipy import sparse
from full_retained_shared_entity_fit_fixtures_v2 import setup
from check_full_shared_entity_fits_v2 import settings
from check_retained_shared_entity_candidates import inputs,identity
from check_full_retained_shared_entity_timing_selection import main as selection_main
from unittest.mock import patch
import sys
from reference_measurement_union_sources import verify
from covariance_basis_context import ComponentKernelProducts
from full_retained_shared_entity_timing import probe,estimate
from prepare_full_retained_shared_entity_timing import run,validate_probes,atomic
from readback_full_retained_shared_entity_timing import run as readback
from run_ortholog_pair_guide_comparison import sha


def setup_timing(root):
    fitpath=setup(root);fit=json.loads(fitpath.read_text());fit['launch_state']='not_launched_or_queued';atomic(fitpath,fit)
    path=root/'timing-plan.json';atomic(path,dict(fit_plan=str(fitpath),fit_plan_sha256=sha(fitpath),
        output=str(root/'timing'),scaled_variance_points=[0.,1.],pins={},resources=dict(minimum_free_disk_gib=0),
        scope='Synthetic full timing workflow/source/provenance fixture only.'))
    return path


def numerical():
    fit=settings();timing=dict(scaled_variance_points=[0.,1.]);records=[];representatives={};all_rejected=[]
    for exception,mode in itertools.product([False,True],['signed','unsigned']):
        source,rows,x,y,operators,audit,old,cert=inputs(exception,mode)
        for method in ['ml','reml']:
            record_identity=identity(audit,method)
            group_id='group-'+str(exception)+'-'+mode+'-'+method
            representative=dict(identity=record_identity,matrix=x,response=y,source_audit=audit,
                original_audit=old,certificate=cert,
                rank=(x.shape[1],audit['numerical_audit']['normalized_design_condition_number'],record_identity['candidate_id']))
            record=probe(source,fit,timing,rows,operators,representative,7,group_id)
            assert record['status']=='timed_all_declared_points_agree_only',record
            validate_probes([record],{group_id:representative},{group_id:7},timing,fit)
            records.append(record);representatives[group_id]=representative
    planning=estimate(records,fit)
    assert planning['measured_candidate_coverage']==56 and planning['unmeasured_review_candidate_coverage']==0
    assert planning['production_finish_eta'] is None and planning['mathematical_runtime_bound'] is False
    assert {p['primary_evaluation_budget'] for p in planning['group_planning_costs']}=={1503}
    assert {p['independent_evaluation_budget'] for p in planning['group_planning_costs']}=={1522,1526}
    manual=0.
    for record in records:
        primary_setup=record['source_validation_seconds']+record['primary_constructor_seconds']+2*record['production_qualification_guard_seconds']
        reader_setup=record['source_validation_seconds']+record['reader_qualification_seconds']+record['independent_constructor_seconds']
        primary_time=max(p['primary_evaluation_seconds'] for p in record['points'])
        independent_time=max(p['independent_evaluation_seconds'] for p in record['points'])
        primary=primary_setup+1503*primary_time
        r=len(record['parameter_names'])
        reader=max(primary+reader_setup+4*independent_time,reader_setup+(5+2+4*r+1503)*independent_time)
        manual+=7*(primary+reader)
    np.testing.assert_allclose(planning['conditional_budget_weighted_seconds'],manual,rtol=1e-13,atol=1e-13)
    reviewed=copy.deepcopy(records);reviewed[0]['status']='timing_group_requires_review'
    partial=estimate(reviewed,fit)
    assert partial['measured_candidate_coverage']==49 and partial['unmeasured_review_candidate_coverage']==7
    original=records[0];group_id=original['group_id'];representative=representatives[group_id]
    for change in ['representative','eligible_count','missing_point','nonfinite_time','false_numeric_agreement','changed_condition',
                   'lost_retained_parameter','invalid_normalization','negative_source_time','nonfinite_reader_qualification_time']:
        altered=copy.deepcopy(original)
        if change=='representative':altered['representative']['candidate_id']='foreign'
        elif change=='eligible_count':altered['eligible_candidates']+=1
        elif change=='missing_point':altered['points'].pop()
        elif change=='nonfinite_time':altered['points'][0]['primary_evaluation_seconds']=float('inf')
        elif change=='false_numeric_agreement':altered['points'][0]['objective_absolute_error']=1.
        elif change=='changed_condition':altered['selection_condition']+=1.
        elif change=='lost_retained_parameter':altered['parameter_names'].pop()
        elif change=='invalid_normalization':altered['kernel_normalization'][0]=0.
        elif change=='negative_source_time':altered['source_validation_seconds']=-1.
        else:altered['reader_qualification_seconds']=float('inf')
        try:validate_probes([altered],{group_id:representative},{group_id:7},timing,fit)
        except AssertionError:all_rejected.append(change)
        else:raise AssertionError('Malformed retained numeric timing accepted: '+change)
    return records,all_rejected


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(exist_ok=False);assert not a.receipt.exists()
    selection=a.output/'selection.json'
    with patch.object(sys,'argv',['check_full_retained_shared_entity_timing_selection.py','--output',str(selection)]):selection_main()
    verify(json.loads(selection.read_text())['source_hashes'])
    records,numeric_rejections=numerical();grid_rejections=[]
    root=a.output/'fixture';root.mkdir()
    path=setup_timing(root);result=run(path);readback(path,root/'readback.json')
    assert result['candidate_rows']==6000 and result['unique_cohorts']==5 and result['model_setting_rows']==600
    try:run(path)
    except AssertionError:pass
    else:raise AssertionError('Completed timing restarted')
    output=root/'timing';receipt_path=output/'receipt.json';original_receipt=receipt_path.read_bytes()
    manifest_path=output/'cohort_manifest.json';original_manifest=manifest_path.read_bytes()
    first=json.loads(original_manifest)[0];fp=output/first['census_path'];original_census=fp.read_bytes()
    cp=output/first['receipt_path'];original_checkpoint=cp.read_bytes()
    for change in ['omit_case','duplicate_case','change_candidate','promote_source_status','invent_timing_group','false_scientific_acceptance']:
        with gzip.open(fp,'rt') as f:rows=[json.loads(line) for line in f]
        if change=='omit_case':rows.pop()
        elif change=='duplicate_case':rows[-1]=copy.deepcopy(rows[0])
        elif change=='change_candidate':rows[0]['candidate_id']='foreign'
        elif change=='promote_source_status':rows[0]['source_disposition']='qualified_exact_retained_uniform_covariance_basis'
        elif change=='invent_timing_group':rows[0]['group_id']='foreign'
        else:rows[0]['scientific_eligibility']=True
        with gzip.open(fp,'wt') as f:
            for row in rows:f.write(json.dumps(row)+'\n')
        checkpoint=json.loads(original_checkpoint);checkpoint['census_sha256']=sha(fp);checkpoint['census_rows']=len(rows);atomic(cp,checkpoint)
        manifest=json.loads(original_manifest);manifest[0].update(census_sha256=sha(fp),receipt_sha256=sha(cp));atomic(manifest_path,manifest)
        receipt=json.loads(original_receipt)
        for changed in [fp,cp,manifest_path]:receipt['artifacts'][str(changed.relative_to(output))]=sha(changed)
        atomic(receipt_path,receipt)
        try:readback(path,root/'malformed-readback.json')
        except AssertionError:grid_rejections.append(change)
        else:raise AssertionError('Rehashed timing census corruption accepted: '+change)
        fp.write_bytes(original_census);cp.write_bytes(original_checkpoint);manifest_path.write_bytes(original_manifest);receipt_path.write_bytes(original_receipt)
    interrupted=root/'interrupted';interrupted.mkdir();path2=setup_timing(interrupted)
    try:run(path2,stop_after_cohorts=2)
    except InterruptedError:pass
    else:raise AssertionError('Timing checkpoint interruption not exercised')
    checkpoints=list((interrupted/'timing/cohorts').glob('*.receipt.json'));assert len(checkpoints)==2
    hashes={str(p):sha(p) for p in (interrupted/'timing/cohorts').glob('*') if p.is_file()}
    run(path2);readback(path2,interrupted/'readback.json');assert all(sha(p)==d for p,d in hashes.items())
    paths=[Path(__file__),Path('scripts/full_retained_shared_entity_timing.py'),Path('scripts/prepare_full_retained_shared_entity_timing.py'),
        Path('scripts/readback_full_retained_shared_entity_timing.py'),Path('scripts/full_retained_shared_entity_fit_sources.py'),
        Path('scripts/shared_entity_likelihood.py'),Path('scripts/independent_shared_entity_likelihood.py'),
        Path('scripts/independent_shared_entity_likelihood_fast.py'),Path('scripts/covariance_basis_context.py')]
    result=dict(status='passed_full_retained_shared_entity_timing_grid_checkpoint_numerical_probe_contracts_v1',
        synthetic_cases=6000,synthetic_cohorts=5,synthetic_original_model_settings=600,
        ml_reml_two_point_q4_q5_probes_passed=True,actual_numeric_groups=8,actual_numeric_points=16,
        source_validation_and_both_producer_guards_and_reader_latent_guard_costs_included=True,
        synthetic_full_eligible_selection_census=1200,synthetic_selection_groups=40,full_fitting_budgets_preserved=True,
        unmeasured_review_groups_retained=True,completed_restart_refused=True,
        interrupted_checkpoints_reused_without_rewrite=True,rehashed_census_alterations_rejected=grid_rejections,
        malformed_numeric_probe_exports_rejected=numeric_rejections,source_and_journal_fixtures_synthetic=True,
        source_hashes={str(p):sha(p) for p in paths},scientific_eligibility=False,production_fitting_launched=False,
        scope='Complete declared synthetic6000candidate five-cohort source/provenance/checkpoint timing census and separate1200candidate40group eligible/mixed/excluded deterministic selection test; eight actual q4/q5 both-mode ML/REML two-point probes. Parent source/journal fixtures synthetic; unchanged budgets include separate inherited/fresh/latent qualification costs. No biological pilot, real full-scope production timing, runtime bound or inferential acceptance.')
    probe_file=a.output/'numeric_probes.json';atomic(probe_file,records)
    result['source_hashes'][str(selection)]=sha(selection);result['source_hashes'][str(probe_file)]=sha(probe_file)
    for fp in a.output.rglob('*'):
        if fp.is_file():result['source_hashes'][str(fp)]=sha(fp)
    verify(result['source_hashes']);atomic(a.receipt,result)
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
