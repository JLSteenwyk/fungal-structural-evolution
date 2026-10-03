#!/usr/bin/env python3
"""Full-grid checkpoint/link contracts and independent numerical-fit audit tests."""
import argparse
import copy
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
from scipy import sparse
from check_full_covariance_qualification import setup as qualification_setup,close_fixture,write
from prepare_full_covariance_qualification import run as qualification_run
from readback_full_covariance_qualification import run as qualification_readback
from full_covariance_qualification_sources import SUMMARY_FIELDS as QUALIFICATION_SUMMARY
from full_expanded_model_design_sources import array_digest,OUTCOMES,fit_id
from prepare_full_shared_entity_fits import run,candidate
from readback_full_shared_entity_fits_v2 import numeric
from shared_entity_likelihood import SharedEntityLikelihood
from fit_shared_entity_likelihood import fit_shared_entity
from run_ortholog_pair_guide_comparison import sha


def settings():
    return dict(independent_backend='component_spectral_v1',optimizer=dict(maximum_scaled_variance=1e6,max_iterations=200,max_evaluations=500,gradient_tolerance=3e-6),
        independent_audit=dict(column_batch=3,replay=dict(objective_atol=1e-7,gradient_atol=3e-6),
            curvature=dict(coordinate_step=1e-4,gradient_atol=3e-6),
            optimizer=dict(max_iterations=200,max_evaluations=500,gradient_tolerance=3e-6,objective_tolerance=1e-7)))


def setup(root):
    upstream=root/'qualification';upstream.mkdir();qpath=qualification_setup(upstream)
    qplan=json.loads(qpath.read_text());droot=Path(json.loads(Path(qplan['design_plan']).read_text())['output'])
    ip=Path(json.loads(Path(qplan['inputs_plan']).read_text())['output'])
    import pyarrow.parquet as pq
    arrays={}
    for part in json.loads((ip/'partition_manifest.json').read_text()):
        t=pq.read_table(ip/part['path']);arrays[part['mask'],part['order_contrast']]={k:np.asarray(t[k].to_pylist()) for k in OUTCOMES}
    manifest={r['cohort_id']:r for r in json.loads((droot/'cohort_manifest.json').read_text())}
    designs=[json.loads(line) for line in (droot/'unique_designs.jsonl').open()];fits=[]
    for design in designs:
        cohort=manifest[design['cohort_id']]
        with np.load(droot/cohort['path']) as a:rows=a['case_rows']
        for outcome in OUTCOMES:
            y=arrays[cohort['mask'],design['order_contrast']][outcome][rows];digest=array_digest(y,'<f8')
            status=design['disposition'] if design['disposition']!='full_rank_design' else 'constant_response_requires_review' if len(set(y.tolist()))==1 else 'ready_for_working_covariance_fit'
            fits.append(dict(fit_input_id=fit_id(design['design_id'],outcome,digest),design_id=design['design_id'],
                cohort_id=design['cohort_id'],outcome=outcome,response_sha256=digest,records=len(rows),disposition=status))
    (droot/'unique_fit_inputs.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in fits))
    receipt=json.loads((droot/'receipt.json').read_text());receipt['artifacts']['unique_fit_inputs.jsonl']=sha(droot/'unique_fit_inputs.jsonl')
    write(droot/'receipt.json',receipt)
    closed=json.loads(Path(qplan['design_completion']).read_text());summary={k:closed[k] for k in ['logical_cases','model_setting_rows','unique_cohorts','unique_designs','trees']}
    close_fixture(droot,'complete_verified_full_expanded_model_designs',droot/'receipt.json',summary,{})
    qplan['design_completion']=str(droot/'completion.json');write(qpath,qplan)
    produced=qualification_run(qpath);qualification_readback(qpath,upstream/'output/readback.json')
    summary={k:produced[k] for k in QUALIFICATION_SUMMARY}
    completed=close_fixture(upstream/'output','complete_verified_full_uniform_covariance_qualification',
        upstream/'output/receipt.json',summary,{str(upstream/'output/readback.json'):sha(upstream/'output/readback.json')})
    plan=root/'fit-plan.json';write(plan,dict(qualification_plan=str(qpath),qualification_completion=str(completed),
        expected=dict(model_setting_rows=produced['model_setting_rows']),methods=['ml','reml'],loading_modes=['signed','unsigned'],
        trees=produced['trees'],pins={},output=str(root/'output'),resources=dict(minimum_free_disk_gib=0),
        **settings(),scope='Complete synthetic full-grid/checkpoint/closure fixtures only; no production biological fit or pilot.'))
    return plan


def invoke(plan,output,good=True):
    result=subprocess.run([sys.executable,'scripts/readback_full_shared_entity_fits_v2.py','--plan',str(plan),'--output',str(output)],capture_output=True,text=True)
    if good and result.returncode:raise RuntimeError(result.stdout+result.stderr)
    if not good:assert result.returncode!=0,'Altered fit export accepted'


def mutate(root,name):
    output=root/'output';rp=output/'receipt.json';receipt=json.loads(rp.read_text())
    if name.startswith('link_'):
        fp=output/'setting_fit_links.tsv.gz'
        with gzip.open(fp,'rt') as f:r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;rows=list(r)
        if name=='link_omit_reml':rows=[r for r in rows if r['likelihood_method']!='reml']
        elif name=='link_wrong_candidate':rows[0]['candidate_id']='wrong'
        elif name=='link_wrong_status':rows[0]['candidate_disposition']='optimized_shared_entity_candidate_pending_independent_audit'
        else:raise AssertionError(name)
        with gzip.open(fp,'wt') as f:w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
        receipt['artifacts'][fp.name]=sha(fp)
    else:
        manifest=json.loads((output/'cohort_manifest.json').read_text());part=manifest[0];fp=output/part['path']
        with gzip.open(fp,'rt') as f:rows=[json.loads(line) for line in f]
        if name=='omit_last_tree':rows=[r for r in rows if r['tree']!='t5']
        elif name=='omit_unsigned':rows=[r for r in rows if r['loading_mode']!='unsigned']
        elif name=='duplicate_candidate':rows.append(copy.deepcopy(rows[0]))
        elif name=='wrong_fit_input':rows[0]['fit_input_id']='wrong'
        elif name=='wrong_response_hash':rows[0]['response_sha256']='wrong'
        elif name=='wrong_covariance_audit':rows[0]['covariance_audit_sha256']='wrong'
        elif name=='wrong_active_columns':rows[0]['active_column_indices']=[]
        elif name=='false_scientific_acceptance':rows[0]['scientific_eligibility']=True
        elif name=='promote_source_review':rows[0].update(disposition='optimized_shared_entity_candidate_pending_independent_audit',numerical_attempted=True)
        else:raise AssertionError(name)
        with gzip.open(fp,'wt') as f:f.write(''.join(json.dumps(r)+'\n' for r in rows))
        part['sha256']=sha(fp);part['candidate_rows']=len(rows)
        cp=output/part['receipt_path'];saved=json.loads(cp.read_text());saved['candidate_sha256']=sha(fp);saved['candidate_rows']=len(rows)
        from collections import Counter
        saved['status_counts']=dict(Counter(r['disposition'] for r in rows));write(cp,saved);part['receipt_sha256']=sha(cp)
        write(output/'cohort_manifest.json',manifest)
        for name in [part['path'],part['receipt_path'],'cohort_manifest.json']:receipt['artifacts'][name]=sha(output/name)
    write(rp,receipt)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    # Established two-basis optimizer fixture exercises full numerical-audit
    # acceptance and failures separately from the complete review-heavy grid.
    rng=np.random.default_rng(473);n=36;labels=np.repeat(np.arange(3),12);operators={}
    for name in ['gene','family']:
        z=np.zeros((n,9))
        for column in range(9):z[labels==column%3,column]=rng.normal(size=12)
        operators[name]=sparse.csr_matrix(z)
    factor=rng.normal(size=(n,3));x=np.column_stack([np.ones(n),rng.normal(size=n)])
    y=operators['gene']@rng.normal(size=9)+.3*(operators['family']@rng.normal(size=9))+rng.normal(size=n)
    source=dict(labels=labels,factors={'tree':factor});plan=settings();rows=np.arange(n)
    identity=dict(source_combined_disposition='qualified_uniform_working_covariance_basis',method='reml',tree='tree',
        loading_mode='signed',scientific_eligibility=False,candidate_id='explicit-synthetic-numerical-contract')
    value=candidate(source,plan,rows,identity,x,y,operators)
    audit=numeric(source,plan,rows,identity,x,y,value,operators)
    assert audit['disposition']=='independently_audited_working_candidate_pending_inferential_calibration',audit
    altered_numerical=[]
    for name in ['objective','start_objective','variance_component','false_failure','source_coefficient','injected_failed_start','changed_start_seed']:
        false=copy.deepcopy(value)
        if name=='objective':false['fit']['negative_profiled_likelihood']+=1.
        elif name=='start_objective':false['fit']['starts'][0]['objective']+=1.
        elif name=='variance_component':false['fit']['variance_components'][0]+=1.
        elif name=='source_coefficient':false['fit']['beta'][0]+=1.
        elif name=='injected_failed_start':false['fit']['starts'][0]['disposition']='failed_search_requires_review'
        elif name=='changed_start_seed':false['fit']['starts'][1]['initial_coordinates'][0]+=.1
        else:false.update(fit=None,disposition='shared_entity_fit_error_requires_review',error_type='ArithmeticError',error_message='fake')
        try:numeric(source,plan,rows,identity,x,y,false,operators)
        except (AssertionError,ValueError,ArithmeticError):altered_numerical.append(name)
        else:raise AssertionError('Altered numerical fit accepted: '+name)
    exhausted=copy.deepcopy(plan);exhausted['optimizer'].update(max_iterations=1,max_evaluations=1)
    failure=candidate(source,exhausted,rows,identity,x,y,operators)
    assert numeric(source,exhausted,rows,identity,x,y,failure,operators)['disposition']=='original_numerical_failure_reproduced_requires_review'
    # Full five-entity-plus-species variance recipe, including its naturally
    # unresolved strict optimization status. Never relax the tolerance to
    # manufacture a passing fit.
    rng=np.random.default_rng(196);n=48;labels=np.repeat(np.arange(4),12);five={}
    for name in ['background_node','model_pair','gene','model']:
        z=np.zeros((n,12))
        for column in range(12):z[labels==column%4,column]=rng.normal(size=12)
        five[name]=sparse.csr_matrix(z)
    five['family_intercept']=sparse.csr_matrix((np.ones(n),(np.arange(n),labels)),shape=(n,4))
    factor=rng.normal(size=(n,3));x=np.column_stack([np.ones(n),rng.normal(size=n)])
    y=five['gene']@rng.normal(size=12)+rng.normal(size=n);source=dict(labels=labels,factors={'tree':factor});rows=np.arange(n)
    full_basis_candidate=candidate(source,plan,rows,identity,x,y,five)
    from covariance_basis_context import ComponentKernelProducts
    source_audit=dict(numerical_audit=ComponentKernelProducts(labels,five,np.ones(n)).tree(factor).audit(x))
    full_basis_replay=numeric(source,plan,rows,identity,x,y,full_basis_candidate,five,source_audit)
    assert full_basis_replay['disposition']=='producer_candidate_requires_review_independent_numeric_replay_passed'
    changes=['omit_last_tree','omit_unsigned','duplicate_candidate','wrong_fit_input','wrong_response_hash',
        'wrong_covariance_audit','wrong_active_columns','false_scientific_acceptance','promote_source_review',
        'link_omit_reml','link_wrong_candidate','link_wrong_status']
    with tempfile.TemporaryDirectory(prefix='full-shared-entity-fitting-contract-') as directory:
        root=Path(directory);path=setup(root);result=run(path);invoke(path,root/'readback.json')
        invoke(path,root/'duplicate-readback.json',False)
        assert result['model_setting_rows']==720 and result['candidate_rows']==7200 and result['setting_fit_links']==14400
        try:run(path)
        except AssertionError:pass
        else:raise AssertionError('Completed fitting restarted')
        snapshot={str(p.relative_to(root/'output')):p.read_bytes() for p in (root/'output').rglob('*') if p.is_file()}
        for change in changes:
            mutate(root,change)
            (root/'output/readback_stage.json').unlink()  # Fresh software observer for the rehashed export.
            (root/'output/independent_readback_completed.json').unlink()
            invoke(path,root/'bad-readback.json',False)
            for p in (root/'output').rglob('*'):
                if p.is_file() and str(p.relative_to(root/'output')) not in snapshot:p.unlink()
            for name,data in snapshot.items():(root/'output'/name).write_bytes(data)
        interrupted=root/'interrupted';interrupted.mkdir();path2=setup(interrupted)
        try:run(path2,stop_after_cohorts=2)
        except InterruptedError:pass
        else:raise AssertionError('Checkpoint interruption not exercised')
        cpfiles=list((interrupted/'output/cohorts').glob('*.jsonl.gz'));assert len(cpfiles)==2
        checkpoint_hashes={str(p):sha(p) for p in cpfiles};run(path2);invoke(path2,interrupted/'readback.json')
        assert all(sha(p)==digest for p,digest in checkpoint_hashes.items())
    sources=['scripts/check_full_shared_entity_fits_v2.py','scripts/full_shared_entity_fit_sources.py',
        'scripts/prepare_full_shared_entity_fits.py','scripts/readback_full_shared_entity_fits_v2.py',
        'scripts/independent_shared_entity_optimizer.py','scripts/independent_shared_entity_likelihood.py',
        'scripts/independent_shared_entity_likelihood_fast.py',
        'scripts/shared_entity_likelihood.py','scripts/fit_shared_entity_likelihood.py',
        'scripts/full_covariance_qualification_sources.py','scripts/full_expanded_model_design_sources.py']
    receipt=dict(status='passed_full_shared_entity_fitting_grid_checkpoint_and_spectral_contracts',
        independent_backend='component_spectral_v1',synthetic_original_model_settings=720,synthetic_candidate_rows=7200,synthetic_setting_fit_links=14400,
        rehashed_altered_grid_exports_rejected=changes,altered_numerical_exports_rejected=altered_numerical,
        independent_candidate_curvature_and_multistart_path_passed=True,budget_failure_reproduced_without_promotion=True,
        completed_restart_refused=True,completed_reader_alternate_output_restart_refused=True,
        full_five_entity_plus_species_review_candidate_independently_replayed=True,
        interrupted_cohort_checkpoints_reused_without_rewrite=True,
        source_and_journal_fixtures_synthetic=True,source_hashes={p:sha(p) for p in sources},
        scope='Complete synthetic 24-case/six-cohort/720-setting source grid under both outcomes, '
            'ML/REML, signed/unsigned and all five trees, plus explicit positive numerical-fit and '
            'failure/export contracts. No biological pilot, production fitting launch, closed '
            'production numerical result or calibrated inference.')
    with a.output.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt))


if __name__=='__main__':main()
