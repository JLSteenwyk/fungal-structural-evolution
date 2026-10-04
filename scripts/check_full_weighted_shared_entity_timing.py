#!/usr/bin/env python3
"""Complete four-control timing grids and independent numeric probes, without mocks."""
import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha
from check_full_shared_entity_fits_v2 import settings
from check_full_weighted_shared_entity_fits_v2 import rejected
from check_retained_shared_entity_candidates import inputs
from check_weighted_shared_entity_candidates import source_grids
from full_weighted_fit_exports import atomic
from full_weighted_shared_entity_fit_sources import load, cohorts
from full_weighted_shared_entity_timing import groups, probe, estimate
from full_weighted_timing_contracts import validate_probes, semantic, same_numeric
from prepare_full_weighted_shared_entity_timing import run
from readback_full_weighted_shared_entity_timing import run as readback
from reference_measurement_union_sources import verify


def actual(root):
    path=Path('data/software_audits/weighted-shared-entity-candidates-20261004-v1/actual_weighted_candidate_cases.json')
    gate=Path('metadata/weighted_shared_entity_candidate_software_validation_20261004_v1.json')
    verify(json.loads(gate.read_text())['source_hashes'])
    fit=settings(); timing=dict(scaled_variance_points=[0.,1.]); records=[]; rejected_exports=[]
    for row in json.loads(path.read_text()):
        source,rows,x,y,ops,retained,old,cert=inputs(row['pair_exception'],row['loading_mode'])
        operators={} if row['policy']=='uniform' else {'target_node':sparse.eye(len(rows),format='csr')}
        operators.update(ops); y=y if row['outcome']=='rmsd_delta' else np.tanh(y)+.15*x[:,1]
        audit=row['audit']; identity={k:v for k,v in row['candidate'].items()
            if k not in ['disposition','fit','numerical_attempted','error_type','error_message']}
        route=dict(names=audit['retained_kernel_names'],certificate_sha256=audit['exact_certificate_sha256'],
            route=audit['basis_route'],exact_uniform_one=audit['residual_diagonal_is_exact_uniform_one'])
        r=dict(identity=identity,matrix=x,response=y,diagonal=np.asarray(row['diagonal']),source_audit=audit,
            route=route,rank=(x.shape[1],audit['numerical_audit']['normalized_design_condition_number'],identity['candidate_id']))
        gid=identity['candidate_id']; value=probe(source,fit,timing,rows,operators,r,7,gid)
        assert value['status']=='timed_all_declared_points_agree_only',value
        validate_probes([value],{gid:r},{gid:7},timing,fit)
        records.append(value)
        if len(records)==1:
            for change in ['foreign_identity','wrong_control','changed_D','count','condition','lost_point','negative_time',
                           'nan_time','false_agreement','extra_field','missing_parameter','negative_normalization','promote_science',
                           'invented_precision','invented_construction']:
                altered=deepcopy(value)
                if change=='foreign_identity':altered['representative']['candidate_id']='foreign'
                elif change=='wrong_control':altered['representative']['control_policy']='foreign'
                elif change=='changed_D':altered['representative']['diagonal_sha256']='foreign'
                elif change=='count':altered['eligible_candidates']+=1
                elif change=='condition':altered['selection_condition']+=1
                elif change=='lost_point':altered['points'].pop()
                elif change=='negative_time':altered['source_validation_seconds']=-1
                elif change=='nan_time':altered['reader_qualification_seconds']=float('nan')
                elif change=='false_agreement':altered['points'][0]['maximum_coordinate_gradient_error']=1
                elif change=='extra_field':altered['foreign']=False
                elif change=='missing_parameter':altered['parameter_names'].pop()
                elif change=='negative_normalization':altered['kernel_normalization'][0]=-1
                elif change=='promote_science':altered['scientific_eligibility']=True
                elif change=='invented_precision':
                    altered['points'][0]=dict(scaled_variance=0.,status='timing_probe_precision_requires_review',
                        error_type='ArithmeticError',error_message='invented',elapsed_seconds=.1)
                    altered['status']='timing_group_requires_review'
                else:
                    altered={k:v for k,v in altered.items() if k in ['group_id','representative','eligible_candidates',
                        'selection_active_columns','selection_condition','variance_points','scientific_eligibility']}
                    altered.update(status='timing_group_construction_requires_review',error_type='ValueError',error_message='invented')
                def action():
                    validate_probes([altered],{gid:r},{gid:7},timing,fit)
                    same_numeric(semantic(value),semantic(altered))
                rejected(action);rejected_exports.append(change)
    assert len(records)==64
    planning=estimate(records,fit)
    assert planning['measured_candidate_coverage']==448 and planning['unmeasured_review_candidate_coverage']==0
    assert {p['primary_evaluation_budget'] for p in planning['group_planning_costs']}=={1503}
    assert {p['independent_evaluation_budget'] for p in planning['group_planning_costs']}=={1522,1526,1530}
    manual=0.
    for r in records:
        primary=r['source_validation_seconds']+r['primary_constructor_seconds']+2*r['production_qualification_guard_seconds']
        reader=r['source_validation_seconds']+r['reader_qualification_seconds']+r['independent_constructor_seconds']
        primary+=1503*max(v['primary_evaluation_seconds'] for v in r['points'])
        it=max(v['independent_evaluation_seconds'] for v in r['points'])
        manual+=7*(primary+max(primary+reader+4*it,reader+(1503+7+4*len(r['parameter_names']))*it))
    np.testing.assert_allclose(planning['conditional_budget_weighted_seconds'],manual,rtol=1e-13)
    reviewed=deepcopy(records);reviewed[0]['status']='timing_group_requires_review'
    assert estimate(reviewed,fit)['unmeasured_review_candidate_coverage']==7
    atomic(root/'actual_64_probes.json',records);atomic(root/'actual_64_conditional_planning.json',planning)
    return {str(p):sha(p) for p in [path,gate,root/'actual_64_probes.json',root/'actual_64_conditional_planning.json']},rejected_exports


def corruptions(plan,root):
    rp=root/'receipt.json';mp=root/'cohort_manifest.json';manifest=json.loads(mp.read_text())
    part=next(m for m in manifest if json.loads((root/m['probes_path']).read_text()))
    fp=root/part['census_path'];pp=root/part['probes_path'];cp=root/part['receipt_path']
    paths=[rp,mp,fp,pp,cp,root/'conditional_planning.json']; originals={p:p.read_bytes() for p in paths}; failures=[]
    for name in ['omit_candidate','duplicate_candidate','control_identity','selected_candidate','eligible_count','missing_probe',
                 'false_success','invented_review','changed_diagonal','changed_planning','missing_artifact','extra_artifact',
                 'source_hash','stage_contract']:
        census=[json.loads(l) for l in gzip.decompress(originals[fp]).decode().splitlines()]
        probes=json.loads(originals[pp]);checkpoint=json.loads(originals[cp]); receipt=json.loads(originals[rp]);current=json.loads(originals[mp])
        if name=='omit_candidate':census.pop()
        elif name=='duplicate_candidate':census[-1]=deepcopy(census[0])
        elif name=='control_identity':census[0]['identity_sha256']='foreign'
        elif name=='selected_candidate':probes[0]['representative']['candidate_id']='foreign'
        elif name=='eligible_count':probes[0]['eligible_candidates']+=1
        elif name=='missing_probe':probes.pop()
        elif name=='false_success':
            probes[0]['points'][0]['maximum_coordinate_gradient_error']=1.
            # Make a structurally consistent but false numerical review.
            probes[0]['points'][0]['status']='timing_probe_numerical_agreement_requires_review'
            probes[0]['status']='timing_group_requires_review'
        elif name=='invented_review':
            probes[0]={k:v for k,v in probes[0].items() if k in ['group_id','representative','eligible_candidates',
                'selection_active_columns','selection_condition','variance_points','scientific_eligibility']}
            probes[0].update(status='timing_group_construction_requires_review',error_type='ValueError',error_message='invented')
        elif name=='changed_diagonal':probes[0]['representative']['diagonal_sha256']='foreign'
        elif name=='changed_planning':
            planning=json.loads(originals[root/'conditional_planning.json']);planning['conditional_budget_weighted_seconds']+=1
            (root/'conditional_planning.json').write_text(json.dumps(planning)+'\n')
        elif name=='missing_artifact':receipt['artifacts'].pop(part['census_path'])
        elif name=='extra_artifact':receipt['artifacts']['invented']='foreign'
        elif name=='source_hash':receipt['source_hashes'][next(iter(receipt['source_hashes']))]='foreign'
        else:receipt['fit_contract']='foreign'
        fp.write_bytes(gzip.compress((''.join(json.dumps(r)+'\n' for r in census)).encode(),mtime=0))
        pp.write_text(json.dumps(probes)+'\n');checkpoint.update(census_sha256=sha(fp),probes_sha256=sha(pp));cp.write_text(json.dumps(checkpoint)+'\n')
        entry=next(m for m in current if m['cohort_id']==part['cohort_id'])
        entry.update(census_sha256=sha(fp),probes_sha256=sha(pp),receipt_sha256=sha(cp));mp.write_text(json.dumps(current)+'\n')
        for p in [fp,pp,cp,mp,root/'conditional_planning.json']:
            key=p.relative_to(root).as_posix()
            if key in receipt['artifacts']:receipt['artifacts'][key]=sha(p)
        rp.write_text(json.dumps(receipt)+'\n')
        rejected(lambda:readback(plan,root/'rejected-readback.json'));failures.append(name)
        assert not (root/'rejected-readback.json').exists() and not (root/'reader_completed.json').exists()
        for p,b in originals.items():p.write_bytes(b)
    return failures


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    assert not a.receipt.exists();a.output.mkdir(parents=True,exist_ok=False)
    grids,bindings,source_rejections=source_grids(a.output);completed=[];output_rejections=[];selected_groups=0
    for number in range(3):
        fp=a.output/('fit-plan-'+str(number)+'.json');fit=json.loads(fp.read_text())
        fit.update(launch_state='not_launched_or_queued',output=str(a.output/('unlaunched-fits-'+str(number))))
        fp.write_text(json.dumps(fit)+'\n');bindings[str(fp)]=sha(fp)
        root=a.output/('timing-'+str(number));tp=a.output/('timing-plan-'+str(number)+'.json')
        atomic(tp,dict(fit_plan=str(fp),fit_plan_sha256=sha(fp),scaled_variance_points=[0.,1.],output=str(root),pins={},
            resources=dict(minimum_free_disk_gib=0),scope='Complete synthetic closed-source timing grid; real native numerical probes, no mocks or biological pilot.'))
        source,hashes=load(fit,fp);bindings.update(hashes)
        # Independently calculate each maximum tuple from all designs rather
        # than accepting the selected representative exported by the stage.
        from full_weighted_shared_entity_fit_sources import cases
        from weighted_shared_entity_candidate import READY
        for cohort,rows,entries in cohorts(source,fit):
            census,selected,counts=groups(source,fit,cohort,entries);reference={};counter=Counter()
            for identity,x,y,audit,route,d in cases(source,fit,cohort,entries):
                if identity['source_combined_disposition']!=READY:continue
                key=tuple(identity[k] for k in ['control_policy','loading_mode','tree','method','outcome'])
                rank=(x.shape[1],float(audit['numerical_audit']['normalized_design_condition_number']),identity['candidate_id'])
                reference[key]=max(rank,reference.get(key,rank));counter[key]+=1
            assert len(census)==4800 and len(selected)==len(reference)
            assert {r['key']:r['rank'] for r in selected.values()}==reference
            assert {r['key']:counts[g] for g,r in selected.items()}==counter
            selected_groups+=len(selected)
        if number==1:
            try:run(tp,stop_after_cohorts=2)
            except InterruptedError:pass
            else:raise AssertionError('Timing interruption not exercised')
            before={str(f):sha(f) for f in (root/'cohorts').glob('*')}
            assert len(list((root/'cohorts').glob('*.receipt.json')))==2
            # Rehashing a closed checkpoint cannot substitute a foreign D.
            pp=next(p for p in (root/'cohorts').glob('*.probes.json') if json.loads(p.read_text()))
            cp=pp.with_name(pp.name.replace('.probes.json','.receipt.json'));prior={p:p.read_bytes() for p in [pp,cp]}
            altered=json.loads(prior[pp]);altered[0]['representative']['diagonal_sha256']='foreign';pp.write_text(json.dumps(altered)+'\n')
            checkpoint=json.loads(prior[cp]);checkpoint['probes_sha256']=sha(pp);cp.write_text(json.dumps(checkpoint)+'\n')
            rejected(lambda:run(tp));output_rejections.append('rehashed_checkpoint_foreign_diagonal')
            for path,b in prior.items():path.write_bytes(b)
        producer=run(tp)
        if number==1:
            assert all(sha(f)==h for f,h in before.items());output_rejections+=corruptions(tp,root)
        reader=readback(tp,root/'readback.json')
        rejected(lambda:run(tp));rejected(lambda:readback(tp,root/'other-readback.json'))
        verify(reader['source_hashes']);bindings.update(reader['source_hashes'])
        bindings[str(root/'readback.json')]=sha(root/'readback.json');bindings[str(root/'reader_completed.json')]=sha(root/'reader_completed.json')
        completed.append(dict(grid=number,candidate_rows=producer['candidate_rows'],timing_groups=producer['timing_groups'],
            eligible_candidates=producer['eligible_candidates'],timing_status_counts=producer['timing_status_counts'],
            producer=str(root/'receipt.json'),reader=str(root/'readback.json')))
    actual_bindings,probe_rejections=actual(a.output);bindings.update(actual_bindings)
    modules=['check_full_weighted_shared_entity_timing','full_weighted_shared_entity_timing','full_weighted_timing_contracts',
        'prepare_full_weighted_shared_entity_timing','readback_full_weighted_shared_entity_timing','full_weighted_shared_entity_fit_sources',
        'weighted_shared_entity_candidate','independent_positive_diagonal_basis_context','shared_entity_likelihood',
        'independent_shared_entity_likelihood','independent_shared_entity_likelihood_fast','check_weighted_shared_entity_candidates',
        'check_full_weighted_shared_entity_fits_v2','full_weighted_fit_exports','reference_measurement_union_sources']
    bindings.update({'scripts/'+m+'.py':sha('scripts/'+m+'.py') for m in modules});verify(bindings)
    result=dict(status='passed_complete_four_control_timing_census_native_numeric_and_accounting_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),complete_synthetic_grids=completed,total_candidate_rows=72000,
        full_grid_selected_groups=selected_groups,actual_additional_numeric_groups=64,actual_additional_numeric_points=128,
        full_grid_native_probes_mocked=False,source_and_journal_fixtures_synthetic=True,
        rehashed_output_alterations_rejected=output_rejections,malformed_probe_exports_rejected=probe_rejections,
        synthetic_source_closure_alterations_rejected=source_rejections,closed_producer_chunks_replayed_without_rewrite=True,
        completed_producer_and_reader_restarts_refused=True,full_grid_selections_independently_recomputed=True,
        primary_evaluation_budget=1503,independent_evaluation_budgets=[1522,1526,1530],production_finish_eta=None,
        source_hashes=bindings,scientific_eligibility=False,fits_computed=0,production_timing_launched=False,
        scope='All three declared five-cohort source grids/72000candidate identities, four controls, both modes/outcomes/methods and all five trees. Every eligible group measured using actual D and independently derived latent/spectral numerics; no native mocks. Separate64common/pair-exception probes cover q4/q5/q6 and128points. Complete census/selection/cost readback, rehashed corruptions, checkpoint replay/refusal. Synthetic source/journal closures; no biological pilot, real production timing or effect accepted.')
    atomic(a.receipt,result);print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','scope']},indent=2),flush=True)


if __name__=='__main__':main()
