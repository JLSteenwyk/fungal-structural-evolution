#!/usr/bin/env python3
"""Full closed-grid fitting admission/capacity checks without native fitting."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import shutil

from ancestral_chain_attempt import sha
from check_full_reduced_covariance_qualification import closed
from check_full_weighted_shared_entity_fits_v2 import rejected
from full_weighted_fit_admission import closed_timing, capacity, production_contract, FIELDS, COMPLETED, FULL
from full_weighted_fit_exports import atomic
from reference_measurement_union_sources import verify


def copied_timing(root,number):
    """Change only fixture custody paths; immutable measured inputs stay bound."""
    original=Path('data/software_audits/full-weighted-shared-entity-timing-20261004-v3/timing-plan-'+str(number)+'.json')
    plan=json.loads(original.read_text());prior_root=Path(plan['output']);folder=root/('copied-timing-'+str(number))
    shutil.copytree(prior_root,folder)
    tp=root/('copied-timing-plan-'+str(number)+'.json');plan['output']=str(folder)
    plan['scope']='Synthetic admission custody fixture copying closed native timing outputs; hardware observations inherited, no new native run or original process closure claim.'
    atomic(tp,plan)
    stage=json.loads((folder/'stage_plan.json').read_text());stage['plan_sha256']=sha(tp);(folder/'stage_plan.json').write_text(json.dumps(stage)+'\n')
    manifest=json.loads((folder/'cohort_manifest.json').read_text())
    for part in manifest:
        cp=folder/part['receipt_path'];value=json.loads(cp.read_text());value['stage']=stage;cp.write_text(json.dumps(value)+'\n')
        part['receipt_sha256']=sha(cp)
    (folder/'cohort_manifest.json').write_text(json.dumps(manifest)+'\n')
    producer=json.loads((folder/'receipt.json').read_text());producer['plan_sha256']=sha(tp)
    producer['source_hashes'].pop(str(original));producer['source_hashes'][str(tp)]=sha(tp)
    producer['artifacts']={name:sha(folder/name) for name in producer['artifacts']}
    (folder/'receipt.json').write_text(json.dumps(producer)+'\n')
    reader=json.loads((folder/'readback.json').read_text());reader['plan_sha256']=sha(tp);reader['producer_receipt_sha256']=sha(folder/'receipt.json')
    reader['source_hashes'].pop(str(original));reader['source_hashes'][str(tp)]=sha(tp)
    for name in list(reader['source_hashes']):
        path=Path(name)
        if path.is_relative_to(prior_root):
            reader['source_hashes'].pop(name);target=folder/path.relative_to(prior_root);reader['source_hashes'][str(target)]=sha(target)
    (folder/'readback.json').write_text(json.dumps(reader)+'\n')
    from prepare_full_weighted_shared_entity_timing_v2 import READER
    (folder/'reader_completed.json').write_text(json.dumps(dict(output=str(folder/'readback.json'),sha256=sha(folder/'readback.json'),status=READER,plan_sha256=sha(tp)))+'\n')
    return tp


def completion(root,number,tp):
    plan=json.loads(tp.read_text());folder=Path(plan['output']);rp=folder/'receipt.json';rb=folder/'readback.json'
    producer=json.loads(rp.read_text());reader=json.loads(rb.read_text());summary={k:reader[k] for k in FIELDS}
    paths=[tp,rp,rb,*[folder/p for p in producer['artifacts']]]
    cp=closed(root,'timing-'+str(number),COMPLETED,paths,summary)
    value=json.loads(cp.read_text());value.update(producer_receipt=str(rp),producer_receipt_sha256=sha(rp),
        independent_readback=str(rb),independent_readback_sha256=sha(rb))
    cp.write_text(json.dumps(value)+'\n');return cp


def corrupt(root,number,tp,cp):
    plan=json.loads(tp.read_text());folder=Path(plan['output']);rp=folder/'receipt.json';rb=folder/'readback.json'
    mp=folder/'cohort_manifest.json';manifest=json.loads(mp.read_text());part=manifest[0]
    fp=folder/part['census_path'];pp=folder/part['probes_path'];ck=folder/part['receipt_path']
    archive=Path(json.loads(cp.read_text())['full_hash_archive'])
    paths=[cp,archive,rp,rb,mp,fp,pp,ck];originals={p:p.read_bytes() for p in paths};bad=[]
    for name in ['unfinished_closure','wrong_completion_status','foreign_fit_contract','changed_reader_source','omit_candidate',
        'duplicate_candidate','wrong_candidate_digest','foreign_representative','missing_group','missing_numeric_replay',
        'promoted_science','altered_hardware_claim','changed_budget','extra_artifact','wrong_stage_hash']:
        prod=json.loads(originals[rp]);read=json.loads(originals[rb]);complete=json.loads(originals[cp]);full=json.loads(originals[archive])
        current=json.loads(originals[mp]);checkpoint=json.loads(originals[ck])
        if name=='unfinished_closure':complete['exact_process_journals_checked']=1
        elif name=='wrong_completion_status':complete['status']='foreign'
        elif name=='foreign_fit_contract':prod['fit_contract']=read['fit_contract']='foreign'
        elif name=='changed_reader_source':read['source_hashes'][next(iter(read['source_hashes']))]='foreign'
        elif name in ['omit_candidate','duplicate_candidate','wrong_candidate_digest']:
            rows=[json.loads(l) for l in gzip.decompress(originals[fp]).decode().splitlines()]
            if name=='omit_candidate':rows.pop()
            elif name=='duplicate_candidate':rows[-1]=deepcopy(rows[0])
            else:rows[0]['identity_sha256']='foreign'
            fp.write_bytes(gzip.compress((''.join(json.dumps(r)+'\n' for r in rows)).encode(),mtime=0))
        elif name in ['foreign_representative','missing_group']:
            # For qualified grids target the third cohort; baseline has no
            # groups, and deleting its original empty export is still invalid.
            if number>0:
                part=manifest[2];pp=folder/part['probes_path'];ck=folder/part['receipt_path'];checkpoint=json.loads(ck.read_text())
                originals.setdefault(pp,pp.read_bytes());originals.setdefault(ck,ck.read_bytes())
                values=json.loads(pp.read_text())
                if name=='foreign_representative':values[0]['representative']['candidate_id']='foreign'
                else:values.pop()
                pp.write_text(json.dumps(values)+'\n')
            else:prod['timing_groups']=read['timing_groups']=1
        elif name=='missing_numeric_replay':read['fresh_numeric_probes_replayed']+=1
        elif name=='promoted_science':prod['scientific_eligibility']=read['scientific_eligibility']=True
        elif name=='altered_hardware_claim':read['independent_hardware_timing_reimplementation']=True
        elif name=='changed_budget':prod['conditional_budget_weighted_seconds']=read['conditional_budget_weighted_seconds']=complete['conditional_budget_weighted_seconds']=1.
        elif name=='extra_artifact':prod['artifacts']['invented']='foreign'
        else:prod['plan_sha256']=read['plan_sha256']='foreign'
        # Rehash every changed local export/checkpoint/manifest and both
        # receipt lineages in the synthetic closure, retaining positive source
        # bindings. The gate must inspect semantics beyond local hashes.
        part=current[2] if number>0 and name in ['foreign_representative','missing_group'] else current[0]
        checkpoint.update(census_sha256=sha(folder/part['census_path']),probes_sha256=sha(folder/part['probes_path']))
        ck.write_text(json.dumps(checkpoint)+'\n')
        for key in ['census','probes','receipt']:part[key+'_sha256']=sha(folder/part[key+'_path'])
        mp.write_text(json.dumps(current)+'\n')
        for path in [fp,pp,ck,mp]:
            key=path.relative_to(folder).as_posix()
            if key in prod['artifacts']:prod['artifacts'][key]=sha(path)
            if str(path) in read['source_hashes']:read['source_hashes'][str(path)]=sha(path)
        rp.write_text(json.dumps(prod)+'\n');read['producer_receipt_sha256']=sha(rp);read['source_hashes'][str(rp)]=sha(rp)
        rb.write_text(json.dumps(read)+'\n');complete.update(producer_receipt_sha256=sha(rp),independent_readback_sha256=sha(rb))
        for path in paths+list(originals):
            if str(path) in full['source_hashes']:full['source_hashes'][str(path)]=sha(path)
        if name=='changed_budget':full['summary']['conditional_budget_weighted_seconds']=1.
        archive.write_text(json.dumps(full)+'\n');complete['full_hash_archive_sha256']=sha(archive);cp.write_text(json.dumps(complete)+'\n')
        rejected(lambda:closed_timing(plan['fit_plan'],tp,cp));bad.append([number,name])
        for path,value in originals.items():path.write_bytes(value)
        part=manifest[0];fp=folder/part['census_path'];pp=folder/part['probes_path'];ck=folder/part['receipt_path']
    return bad


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);gate=Path('metadata/full_weighted_timing_software_validation_20261004_v3.json')
    parent=json.loads(gate.read_text());verify(parent['source_hashes']);bindings={str(gate):sha(gate)};cases=[];bad=[]
    try:
        for number in range(3):
            tp=copied_timing(a.output,number)
            plan=json.loads(tp.read_text());cp=completion(a.output,number,tp)
            source,hashes,admitted=closed_timing(plan['fit_plan'],tp,cp)
            fit=json.loads(Path(plan['fit_plan']).read_text());fit['resources']=dict(minimum_free_disk_gib=3172,output_scratch_reserve_gib=3072)
            resources=capacity(admitted,fit)
            assert resources['producer_native_cpu_allocation_seconds']>=604800 and resources['reader_native_cpu_allocation_seconds']>=604800
            assert resources['production_finish_eta'] is None and resources['resources_installed'] is False
            assert admitted['fresh_admission_candidate_rows']==24000 and admitted['admission_repeats_numeric_probes'] is False
            rejected(lambda:production_contract(fit,admitted))
            bad+=corrupt(a.output,number,tp,cp)
            _,restored,restored_result=closed_timing(plan['fit_plan'],tp,cp)
            assert restored_result==admitted and restored==hashes;verify(restored);bindings.update(restored)
            cases.append(dict(grid=number,candidate_rows=24000,timing_groups=admitted['timing_groups'],
                measured_candidate_coverage=admitted['measured_candidate_coverage'],unmeasured_review_candidate_coverage=admitted['unmeasured_review_candidate_coverage'],
                immutable_measured_fit_plan=str(plan['fit_plan']),source_fit_contract=source['fit_contract'],resources=resources))
        # Production scope/tolerances are separate from synthetic source
        # acceptance. Exercise full headers and reject every scope relaxation.
        pp=Path('metadata/full_weighted_shared_entity_fit_draft_plan_20261004_v1.json');production=json.loads(pp.read_text());verify(production['pins'])
        header=dict(logical_cases=75188,unique_cohorts=4340,candidate_rows=20832000,model_setting_rows=622080)
        production_contract(production,header);scope_bad=[]
        for name in ['partial_expected','partial_admission','drop_policy','drop_mode','drop_tree','drop_method','changed_primary_budget',
            'relaxed_primary_gradient','changed_variance_cap','changed_column_batch','relaxed_independent_gradient','relaxed_replay','changed_curvature']:
            fit=deepcopy(production);h=deepcopy(header)
            if name=='partial_expected':fit['expected']['candidate_rows']-=1
            elif name=='partial_admission':h['candidate_rows']-=1
            elif name=='drop_policy':fit['policies'].pop()
            elif name=='drop_mode':fit['loading_modes'].pop()
            elif name=='drop_tree':fit['trees'].pop()
            elif name=='drop_method':fit['methods'].pop()
            elif name=='changed_primary_budget':fit['optimizer']['max_evaluations']+=1
            elif name=='relaxed_primary_gradient':fit['optimizer']['gradient_tolerance']=3e-6
            elif name=='changed_variance_cap':fit['optimizer']['maximum_scaled_variance']=1e5
            elif name=='changed_column_batch':fit['independent_audit']['column_batch']=3
            elif name=='relaxed_independent_gradient':fit['independent_audit']['optimizer']['gradient_tolerance']=3e-6
            elif name=='relaxed_replay':fit['independent_audit']['replay']['gradient_atol']=3e-6
            else:fit['independent_audit']['curvature']['coordinate_step']=.001
            rejected(lambda:production_contract(fit,h));scope_bad.append(name)
        actual_tp=Path('metadata/full_weighted_shared_entity_timing_plan_20261004_v2.json');actual=json.loads(actual_tp.read_text())
        # Both actual full prerequisite receipts are still absent. Refuse
        # admission before loading source arrays or creating any fit output.
        assert not Path(actual['completion']).exists()
        rejected(lambda:closed_timing(pp,actual_tp,actual['completion']))
        assert not Path(production['output']).exists()
        bindings[str(pp)]=sha(pp);bindings[str(actual_tp)]=sha(actual_tp)
        modules=['full_weighted_fit_admission','check_full_weighted_fit_admission','full_weighted_shared_entity_fit_sources',
            'full_weighted_shared_entity_timing','full_weighted_timing_contracts_v2','prepare_full_weighted_shared_entity_timing_v2',
            'reference_measurement_union_sources','prepare_full_weighted_fit_execution']
        bindings.update({'scripts/'+m+'.py':sha('scripts/'+m+'.py') for m in modules});verify(bindings)
        result=dict(status='passed_complete_closed_weighted_timing_fit_admission_and_capacity_contracts_v1',
            checked_utc=datetime.now(timezone.utc).isoformat(),complete_synthetic_grids=cases,total_admission_candidate_rows=72000,
            rehashed_closure_and_accounting_alterations_rejected=bad,production_scope_and_tolerance_relaxations_rejected=scope_bad,
            actual_missing_prerequisites_refused=True,source_and_journal_closure_fixtures_synthetic=True,
            parent_numeric_probes_closed_not_repeated=True,original_parent_bytes_never_modified=True,own_corrupted_fixture_bytes_restored=True,synthetic_scope_cannot_admit_production=True,
            immutable_measured_configuration_unchanged=True,source_hashes=bindings,production_resources_installed=False,
            production_fitting_launched_or_queued=False,scientific_eligibility=False,new_cost_usd=0,gpu=False,
            scope='All three original five-cohort timing grids/72000candidate identities reconstruct before admission. Actual parent timing numerics already qualified; no fresh probe, native fit or biological pilot. Synthetic compact closures/journals only;45rehashed accounting/closure alterations and13production scope/tolerance relaxations rejected, positive original bytes restored. Actual missing timing prerequisite refuses admission without fit output or resource installation. Capacity budgets remain allocations, not ETAs or numerical acceptance.')
        atomic(a.receipt,result);print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','scope','complete_synthetic_grids']},indent=2),flush=True)
    finally:
        # Negative tests touch only own copied fixture outputs. Published
        # parent bytes remain immutable even if this checker fails.
        verify(parent['source_hashes'])


if __name__=='__main__':main()
