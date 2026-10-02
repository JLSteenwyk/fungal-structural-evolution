"""Closed original full-grid inputs and identities for uniform ML/REML fits."""
import itertools
import json
from pathlib import Path
import numpy as np
import pyarrow.parquet as pq
from background_measurement_union_sources import closed_source
from full_covariance_qualification_sources import load as load_covariance,jsonl,cohort_rows,design_matrix,audit_id
from full_expanded_model_design_sources import array_digest,digest,OUTCOMES,AXES,DEGREES,fit_id,SETTING_FIELDS
from full_expanded_model_input_sources import ORDERS
from full_entity_operator_sources import MODES
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

SUMMARY_FIELDS=['logical_cases','model_setting_rows','unique_cohorts','unique_designs','unique_fit_inputs',
    'candidate_rows','setting_fit_links','candidate_status_counts','setting_status_counts','methods','trees','loading_modes']
FIT_LINK_EXTRA=['likelihood_method','candidate_id','candidate_disposition']


def candidate_id(contract,fid,mode,tree,method):
    return digest(['full-shared-entity-uniform-fit-v1',contract,fid,mode,tree,method])


def load(plan,path):
    assert plan['methods']==['ml','reml'] and plan['loading_modes']==MODES
    bindings=dict(plan['pins']);bind(bindings,path)
    q=closed_source(plan['qualification_completion'],'complete_verified_full_uniform_covariance_qualification',
        'complete_verified_full_uniform_covariance_qualification_archive',2,bindings)
    qpath=Path(plan['qualification_plan']);qplan=json.loads(qpath.read_text());bind(bindings,qpath)
    source,qbindings=load_covariance(qplan,qpath)
    for p,d in qbindings.items():bind(bindings,p,d)
    root=Path(qplan['output'])
    assert q['producer_receipt']==str(root/'receipt.json') and q['trees']==plan['trees']==qplan['trees']
    assert q['scientific_eligibility'] is False and q['loading_modes']==MODES
    assert q['model_setting_rows']==plan['expected']['model_setting_rows']
    rp=json.loads((root/'receipt.json').read_text())
    assert rp['plan_sha256']==sha(qpath) and rp['source_contract']==source['contract']
    assert rp['audit_rows']==q['audit_rows']==q['unique_designs']*len(MODES)*len(plan['trees'])
    input_config=json.loads(Path(qplan['inputs_plan']).read_text());input_root=Path(input_config['output'])
    for part in json.loads((input_root/'partition_manifest.json').read_text()):
        fp=input_root/part['path'];assert sha(fp)==part['sha256']==bindings[str(fp)]
        values=pq.read_table(fp,columns=['case_id',*OUTCOMES])
        assert values['case_id'].to_pylist()==source['ids']
        source['arrays'][part['mask'],part['order_contrast']].update({name:np.asarray(values[name].to_pylist(),dtype=float) for name in OUTCOMES})
    source['fit_contract']=digest(dict(schema='full-shared-entity-uniform-fit-v1',
        qualification_completion=sha(plan['qualification_completion']),qualification_plan=sha(qpath),
        qualification_contract=source['contract'],methods=plan['methods'],optimizer=plan['optimizer'],
        independent_audit=plan['independent_audit'],loading_modes=MODES,trees=plan['trees'],fit_backend_pins=plan['pins']))
    source.update(qualification_root=root,qualification_receipt=rp,qualification_completion=q)
    verify(bindings);return source,bindings


def cohorts(source,plan):
    designs=jsonl(source['root']/'unique_designs.jsonl');fits=jsonl(source['root']/'unique_fit_inputs.jsonl')
    audits=jsonl(source['qualification_root']/'design_covariance_audits.jsonl.gz')
    for cohort in source['cohorts']:
        rows=cohort_rows(source,cohort);entries=[];seen=set()
        for _ in range(len(ORDERS)*len(AXES)*len(DEGREES)):
            design=next(designs);assert design['cohort_id']==cohort['cohort_id']
            seen.add((design['order_contrast'],design['sequence_axis'],design['degree']))
            matrix=design_matrix(source,cohort,rows,design);responses={}
            for outcome in OUTCOMES:
                fit=next(fits);response=source['arrays'][cohort['mask'],design['order_contrast']][outcome][rows]
                response_hash=array_digest(response,'<f8')
                assert fit['design_id']==design['design_id'] and fit['cohort_id']==cohort['cohort_id']
                assert fit['outcome']==outcome and fit['response_sha256']==response_hash
                assert fit['fit_input_id']==fit_id(design['design_id'],outcome,response_hash)
                assert fit['records']==len(rows) and np.isfinite(response).all()
                expected=design['disposition'] if design['disposition']!='full_rank_design' else 'constant_response_requires_review' if len(set(response.tolist()))==1 else 'ready_for_working_covariance_fit'
                assert fit['disposition']==expected
                responses[outcome]=(fit,response)
            qualified={}
            for mode,tree in itertools.product(MODES,plan['trees']):
                record=next(audits)
                assert (record['cohort_id'],record['design_id'],record['loading_mode'],record['tree'])==(cohort['cohort_id'],design['design_id'],mode,tree)
                assert record['audit_id']==audit_id(source['contract'],design['design_id'],mode,tree)
                assert record['raw_design_sha256']==design['raw_design_sha256'] and record['records']==len(rows)
                qualified[mode,tree]=record
            entries.append((design,matrix,responses,qualified))
        assert seen==set(itertools.product(ORDERS,AXES,DEGREES))
        yield cohort,rows,entries
    assert next(designs,None) is None and next(fits,None) is None and next(audits,None) is None


def identity(source,cohort,design,fit,audit,mode,tree,method):
    disposition=fit['disposition'] if fit['disposition']!='ready_for_working_covariance_fit' else audit['disposition']
    return dict(candidate_id=candidate_id(source['fit_contract'],fit['fit_input_id'],mode,tree,method),
        source_contract=source['fit_contract'],cohort_id=cohort['cohort_id'],design_id=design['design_id'],
        fit_input_id=fit['fit_input_id'],outcome=fit['outcome'],loading_mode=mode,tree=tree,method=method,
        covariance_audit_id=audit['audit_id'],covariance_audit_sha256=digest(audit),
        source_fit_disposition=fit['disposition'],source_covariance_disposition=audit['disposition'],
        source_combined_disposition=disposition,records=fit['records'],response_sha256=fit['response_sha256'],
        raw_design_sha256=design['raw_design_sha256'],ordered_case_ids_sha256=cohort['ordered_case_ids_sha256'],
        active_column_indices=design['active_column_indices'],exactly_zero_column_indices=design['exactly_zero_column_indices'],
        coefficient_columns=[design['predictor_columns'][i] for i in design['active_column_indices']],
        scientific_eligibility=False)


def cases(source,plan,cohort,entries):
    for design,matrix,responses,qualified in entries:
        for mode,tree in itertools.product(MODES,plan['trees']):
            audit=qualified[mode,tree]
            for outcome,method in itertools.product(OUTCOMES,plan['methods']):
                fit,y=responses[outcome]
                yield identity(source,cohort,design,fit,audit,mode,tree,method),matrix,y
