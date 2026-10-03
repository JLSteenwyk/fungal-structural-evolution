"""Complete original design/response grid bound to closed exact retained audits."""
import itertools
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from background_measurement_union_sources import closed_source
from full_covariance_qualification_sources import jsonl, folded_operators
from full_expanded_model_design_sources import digest, OUTCOMES, SETTING_FIELDS
from full_reduced_covariance_sources_v3 import load as load_retained, NEW_LINK_EXTRA
from full_shared_entity_fit_sources import load as load_original, cohorts as original_cohorts
from full_shared_entity_fit_sources import identity as original_identity, FIT_LINK_EXTRA, SUMMARY_FIELDS
from reduced_covariance_basis import readback_record
from reference_measurement_union_sources import bind, verify

SCHEMA='full-exact-retained-shared-entity-uniform-fit-v1'
PRODUCER_STATUS='complete_full_retained_shared_entity_fits_pending_independent_readback_v1'
READER_STATUS='passed_full_retained_shared_entity_grid_and_spectral_readback_v1'


def candidate_id(contract,fid,mode,tree,method):
    return digest([SCHEMA,contract,fid,mode,tree,method])


def load(plan,path):
    source,bindings=load_original(plan,path)
    reduced_path=Path(plan['retained_plan']);reduced_plan=json.loads(reduced_path.read_text())
    assert reduced_plan['qualification_plan']==plan['qualification_plan']
    assert reduced_plan['qualification_completion']==plan['qualification_completion']
    assert reduced_plan['completion']==plan['retained_completion']
    retained,extra=load_retained(reduced_plan,reduced_path)
    for p,h in extra.items():bind(bindings,p,h)
    status='complete_verified_full_exact_retained_uniform_covariance_qualification_v3'
    completion=closed_source(plan['retained_completion'],status,status+'_archive',2,bindings)
    root=Path(reduced_plan['output']);rp=root/'receipt.json';rb=root/'readback.json'
    producer=json.loads(rp.read_text());reader=json.loads(rb.read_text())
    assert completion['producer_receipt']==str(rp) and completion['producer_receipt_sha256']==sha(rp)
    assert completion['independent_readback']==str(rb) and completion['independent_readback_sha256']==sha(rb)
    assert producer['status']=='complete_full_exact_retained_uniform_covariance_qualification_pending_readback_v3'
    assert reader['status']=='passed_full_exact_retained_uniform_covariance_qualification_independent_readback_v3'
    assert producer['plan_sha256']==reader['plan_sha256']==sha(reduced_path)
    assert producer['source_contract']==reader['source_contract']==retained['contract']
    assert producer['scientific_eligibility'] is reader['scientific_eligibility'] is False
    assert reader['producer_receipt_sha256']==sha(rp)
    for key in ['logical_cases','model_setting_rows','unique_cohorts','unique_designs','audit_rows','setting_audit_links','trees','loading_modes']:
        assert completion[key]==producer[key]==reader[key]==source['qualification_completion'][key]
    assert producer['retained_basis_audit_counts']==reader['retained_basis_audit_counts']==completion['retained_basis_audit_counts']
    assert set(producer['artifacts'])=={'stage_plan.json','design_covariance_audits.jsonl.gz','setting_audit_links.tsv.gz'}
    for name,h in producer['artifacts'].items():bind(bindings,root/name,h)
    bind(bindings,rp);bind(bindings,rb);bind(bindings,reduced_path)
    source.update(retained_root=root,retained_contract=retained['contract'],certificates=retained['certificates'],
                  original_qualification_root=source['qualification_root'],qualification_root=root,
                  original_qualification_completion=source['qualification_completion'],qualification_completion=completion)
    assert sorted(c['cohort_id'] for c in source['cohorts'])==retained['cohort_ids']
    source['fit_contract']=digest(dict(schema=SCHEMA,original_qualification_completion=sha(plan['qualification_completion']),
        retained_qualification_completion=sha(plan['retained_completion']),retained_plan=sha(reduced_path),
        retained_source_contract=retained['contract'],exact_completion=sha(reduced_plan['exact_completion']),
        methods=plan['methods'],optimizer=plan['optimizer'],independent_audit=plan['independent_audit'],
        independent_backend=plan['independent_backend'],loading_modes=plan['loading_modes'],trees=plan['trees'],
        fit_backend_pins=plan['pins']))
    verify(bindings)
    return source,bindings


def cohorts(source,plan):
    # The frozen original generator still verifies every design, actual response
    # array/hash, fit identity and original audit, including all non-ready rows.
    view=dict(source,qualification_root=source['original_qualification_root'])
    reduced=jsonl(source['retained_root']/'design_covariance_audits.jsonl.gz')
    for cohort,rows,entries in original_cohorts(view,plan):
        updated=[]
        for design,matrix,responses,original in entries:
            qualified={}
            for mode,tree in itertools.product(plan['loading_modes'],plan['trees']):
                old=original[mode,tree];audit=next(reduced)
                cert=source['certificates'][cohort['cohort_id'],mode]
                readback_record(audit,old,cert,source['retained_contract'])
                assert (audit['cohort_id'],audit['design_id'],audit['loading_mode'],audit['tree'])==(cohort['cohort_id'],design['design_id'],mode,tree)
                qualified[mode,tree]=(audit,old,cert)
            updated.append((design,matrix,responses,qualified))
        yield cohort,rows,updated
    assert next(reduced,None) is None


def identity(source,cohort,design,fit,audit,original,certificate,mode,tree,method):
    record=original_identity(source,cohort,design,fit,audit,mode,tree,method)
    record.update(candidate_id=candidate_id(source['fit_contract'],fit['fit_input_id'],mode,tree,method),
        original_covariance_audit_id=original['audit_id'],original_covariance_audit_sha256=digest(original),
        exact_covariance_certificate_sha256=digest(certificate),retained_kernel_names=audit['retained_kernel_names'],
        retained_source_contract=source['retained_contract'],component_variance_attribution_accepted=False,
        nonuniform_weighting_accepted=False)
    return record


def cases(source,plan,cohort,entries):
    for design,matrix,responses,qualified in entries:
        for mode,tree in itertools.product(plan['loading_modes'],plan['trees']):
            audit,original,cert=qualified[mode,tree]
            for outcome,method in itertools.product(OUTCOMES,plan['methods']):
                fit,response=responses[outcome]
                yield identity(source,cohort,design,fit,audit,original,cert,mode,tree,method),matrix,response,audit,original,cert


def operators_for(source,rows,audit):
    names=audit['retained_kernel_names']
    assert names in [['residual','background_node','family_intercept','species'],
                     ['residual','background_node','model_pair','family_intercept','species']]
    original=folded_operators(source,rows,audit['loading_mode'])
    return {name:original[name] for name in names[1:-1]}
