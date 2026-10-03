#!/usr/bin/env python3
"""Check full timing census/restart contracts and actual numerical probes."""
import argparse
import copy
import gzip
import json
from pathlib import Path
import tempfile
import numpy as np
from scipy import sparse
from check_full_shared_entity_fits_v2 import setup,settings
from covariance_basis_context import ComponentKernelProducts
from full_shared_entity_timing import probe,estimate
from prepare_full_shared_entity_timing import run,validate_probes,atomic
from readback_full_shared_entity_timing import run as readback
from run_ortholog_pair_guide_comparison import sha


def setup_timing(root):
    fitpath=setup(root);fit=json.loads(fitpath.read_text());fit['launch_state']='not_launched_or_queued';atomic(fitpath,fit)
    path=root/'timing-plan.json';atomic(path,dict(fit_plan=str(fitpath),fit_plan_sha256=sha(fitpath),
        output=str(root/'timing'),scaled_variance_points=[0.,1.],pins={},resources=dict(minimum_free_disk_gib=0),
        scope='Synthetic full timing workflow/source/provenance fixture only.'))
    return path


def numerical():
    rng=np.random.default_rng(196);n=48;labels=np.repeat(np.arange(4),12);operators={}
    for name in ['background_node','model_pair','gene','model']:
        z=np.zeros((n,12))
        for column in range(12):z[labels==column%4,column]=rng.normal(size=12)
        operators[name]=sparse.csr_matrix(z)
    operators['family_intercept']=sparse.csr_matrix((np.ones(n),(np.arange(n),labels)),shape=(n,4))
    factor=rng.normal(size=(n,3));x=np.column_stack([np.ones(n),rng.normal(size=n)])
    y=operators['gene']@rng.normal(size=12)+rng.normal(size=n)
    guard=ComponentKernelProducts(labels,operators,np.ones(n)).tree(factor).audit(x)
    assert all(guard[k]['disposition']=='numerically_independent_covariance_bases' for k in ['raw_diagnostics','reml_diagnostics'])
    fit=settings();timing=dict(scaled_variance_points=[0.,1.]);source=dict(labels=labels,factors={'tree':factor});records=[]
    for method in ['ml','reml']:
        identity=dict(method=method,tree='tree',loading_mode='signed',candidate_id='synthetic-'+method,scientific_eligibility=False)
        representative=dict(identity=identity,matrix=x,response=y,source_audit=guard,rank=(2,guard['normalized_design_condition_number'],identity['candidate_id']))
        record=probe(source,fit,timing,np.arange(n),operators,representative,7,'group-'+method)
        assert record['status']=='timed_all_declared_points_agree_only',record
        validate_probes([record],{record['group_id']:representative},{record['group_id']:7},timing,fit);records.append(record)
    planning=estimate(records,fit);assert planning['measured_candidate_coverage']==14 and planning['unmeasured_review_candidate_coverage']==0
    assert planning['production_finish_eta'] is None and planning['mathematical_runtime_bound'] is False
    assert all(p['primary_evaluation_budget']==1503 and p['independent_evaluation_budget']==1534 for p in planning['group_planning_costs'])
    reviewed=copy.deepcopy(records);reviewed[0]['status']='timing_group_requires_review'
    partial=estimate(reviewed,fit);assert partial['measured_candidate_coverage']==partial['unmeasured_review_candidate_coverage']==7
    rejected=[];original=records[0]
    for change in ['representative','eligible_count','missing_point','nonfinite_time','false_numeric_agreement','changed_condition']:
        altered=copy.deepcopy(original)
        if change=='representative':altered['representative']['candidate_id']='foreign'
        elif change=='eligible_count':altered['eligible_candidates']+=1
        elif change=='missing_point':altered['points'].pop()
        elif change=='nonfinite_time':altered['points'][0]['primary_evaluation_seconds']=float('inf')
        elif change=='false_numeric_agreement':altered['points'][0]['objective_absolute_error']=1.
        else:altered['selection_condition']+=1.
        identity=original['representative'];rep=dict(identity=identity,matrix=x,rank=(2,guard['normalized_design_condition_number'],identity['candidate_id']))
        try:validate_probes([altered],{altered['group_id']:rep},{altered['group_id']:7},timing,fit)
        except AssertionError:rejected.append(change)
        else:raise AssertionError('Malformed numerical timing accepted: '+change)
    return records,rejected


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    records,numeric_rejections=numerical();grid_rejections=[]
    with tempfile.TemporaryDirectory(prefix='full-shared-entity-timing-contract-') as d:
        root=Path(d);path=setup_timing(root);result=run(path);readback(path,root/'readback.json')
        assert result['candidate_rows']==7200 and result['unique_cohorts']==6 and result['model_setting_rows']==720
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
            elif change=='promote_source_status':rows[0]['source_disposition']='qualified_uniform_working_covariance_basis'
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
    paths=[Path(__file__),Path('scripts/full_shared_entity_timing.py'),Path('scripts/prepare_full_shared_entity_timing.py'),
        Path('scripts/readback_full_shared_entity_timing.py'),Path('scripts/full_shared_entity_fit_sources.py'),
        Path('scripts/shared_entity_likelihood.py'),Path('scripts/independent_shared_entity_likelihood.py'),
        Path('scripts/independent_shared_entity_likelihood_fast.py'),Path('scripts/covariance_basis_context.py')]
    result=dict(status='passed_full_shared_entity_timing_grid_checkpoint_and_numerical_probe_contracts',
        synthetic_cases=7200,synthetic_cohorts=6,synthetic_original_model_settings=720,
        ml_reml_two_point_six_variance_probes_passed=True,full_fitting_budgets_preserved=True,
        unmeasured_review_groups_retained=True,completed_restart_refused=True,
        interrupted_checkpoints_reused_without_rewrite=True,rehashed_census_alterations_rejected=grid_rejections,
        malformed_numeric_probe_exports_rejected=numeric_rejections,source_and_journal_fixtures_synthetic=True,
        source_hashes={str(p):sha(p) for p in paths},scope='Software full-grid/source/provenance/checkpoint contracts and six-variance ML/REML probes only; no biological pilot, production timing result, runtime bound or inferential acceptance.')
    atomic(a.output,result);print(json.dumps(result))


if __name__=='__main__':main()
