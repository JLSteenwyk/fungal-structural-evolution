"""Retained-basis sources with an explicit separately closed parallel readback."""
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import closed_subset
from full_expanded_model_design_sources import digest, SETTING_FIELDS
from full_covariance_qualification_sources import jsonl, LINK_EXTRA
from reduced_covariance_basis import SCHEMA, validate_certificate
from reference_measurement_union_sources import verify
from run_parallel_covariance_readback_v2 import STATUS as PARALLEL_STATUS

NEW_LINK_EXTRA=LINK_EXTRA+['source_audit_id','source_covariance_disposition']
SUMMARY_FIELDS=['logical_cases','model_setting_rows','unique_cohorts','unique_designs','audit_rows','setting_audit_links',
    'audit_status_counts','link_status_counts','trees','loading_modes','retained_basis_audit_counts']


def load(plan,path):
    bindings=dict(plan['pins']);bindings[str(path)]=sha(path)
    old_plan_path=Path(plan['qualification_plan']);old_plan=json.loads(old_plan_path.read_text());old_root=Path(old_plan['output'])
    parallel_path=Path(plan['parallel_readback_plan']);parallel=json.loads(parallel_path.read_text())
    assert parallel['source_plan']==str(old_plan_path)
    independent_path=Path(parallel['output'])/'readback.json'
    assert independent_path!=old_root/'readback.json','Original reader output may never be replaced'
    paths=[old_plan_path,parallel_path,independent_path]+[old_root/n for n in
        ['receipt.json','stage_plan.json','design_covariance_audits.jsonl.gz','setting_audit_links.tsv.gz']]
    old=closed_subset(plan['qualification_completion'],'complete_verified_full_uniform_covariance_qualification',paths,bindings)
    assert old['independent_readback']==str(independent_path) and old['independent_readback_sha256']==sha(independent_path)
    assert old['producer_receipt']==str(old_root/'receipt.json') and old['producer_receipt_sha256']==sha(old_root/'receipt.json')
    reader=json.loads(independent_path.read_text());producer=json.loads((old_root/'receipt.json').read_text())
    assert reader['status']==PARALLEL_STATUS and reader['scientific_eligibility'] is False
    assert reader['execution_backend']==parallel['resources']['execution_backend']=='fork_process'
    assert reader['failure_capture_enabled'] is reader['original_numerical_arithmetic_unchanged'] is True
    assert isinstance(reader['worker_pids'],list) and 1<=len(reader['worker_pids'])<=8
    assert len(set(reader['worker_pids']))==len(reader['worker_pids'])
    assert all(type(pid) is int and pid>0 for pid in reader['worker_pids'])
    failures=Path(parallel['output'])/'failures'
    assert failures.is_dir() and not list(failures.iterdir()),'Rejected arithmetic forbids downstream admission'
    assert reader['plan_sha256']==producer['plan_sha256']==sha(old_plan_path)
    assert reader['parallel_plan_sha256']==sha(parallel_path) and reader['producer_receipt_sha256']==sha(old_root/'receipt.json')
    assert reader['source_contract']==producer['source_contract']
    assert reader['parallel_workers']==parallel['resources']['workers']==8
    assert reader['maximum_observed_pending_cohorts']<=parallel['resources']['maximum_pending_cohorts']==16
    assert reader['actual_cgroup_limits']=={'cpu.max':'800000 100000','memory.max':str(64*2**30),'memory.swap.max':'0'}
    expected=plan['expected']
    for key in ['logical_cases','model_setting_rows','unique_cohorts','unique_designs','audit_rows','setting_audit_links']:
        assert old[key]==reader[key]==producer[key]==expected[key]
    assert old['trees']==reader['trees']==plan['trees'] and old['loading_modes']==reader['loading_modes']==['signed','unsigned']
    assert old['audit_rows']==old['unique_designs']*len(old['trees'])*len(old['loading_modes'])
    assert old['setting_audit_links']==old['model_setting_rows']*len(old['trees'])*len(old['loading_modes'])
    exact_plan_path=Path(plan['exact_plan']);exact_plan=json.loads(exact_plan_path.read_text());exact_root=Path(exact_plan['output'])
    exact_paths=[exact_plan_path]+[exact_root/n for n in ['receipt.json','readback.json','cohort_certificates.jsonl']]
    exact=closed_subset(plan['exact_completion'],'complete_verified_full_exact_uniform_covariance_folds_v2',exact_paths,bindings)
    assert exact['logical_cases']==old['logical_cases'] and exact['cohorts']==old['unique_cohorts']
    certificates={}
    for row in jsonl(exact_root/'cohort_certificates.jsonl'):
        validate_certificate(row);key=(row['cohort_id'],row['certificate']['loading_mode'])
        assert key not in certificates;certificates[key]=row
    assert len(certificates)==exact['certificates']==2*old['unique_cohorts']
    cohort_ids=sorted({k[0] for k in certificates})
    assert set(certificates)=={(c,m) for c in cohort_ids for m in old['loading_modes']}
    contract=digest(dict(schema=SCHEMA,original_qualification_completion=sha(plan['qualification_completion']),
        exact_completion=sha(plan['exact_completion']),original_audits_sha256=bindings[str(old_root/'design_covariance_audits.jsonl.gz')],
        certificates_sha256=bindings[str(exact_root/'cohort_certificates.jsonl')],
        numerical_envelopes='copy_original_principal_submatrices_unchanged',residual='uniform_one'))
    verify(bindings)
    return dict(old=old,old_root=old_root,certificates=certificates,cohort_ids=cohort_ids,contract=contract,
                settings_fields=SETTING_FIELDS),bindings
