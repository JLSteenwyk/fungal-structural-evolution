"""Complete synthetic source closures for retained full-grid fitting contracts."""
from datetime import datetime, timezone
import csv
import gzip
import itertools
import json
from pathlib import Path
from unittest.mock import patch

from ancestral_chain_attempt import sha
from check_full_covariance_qualification import setup as qualification_setup, close_fixture, write
from check_full_shared_entity_fits_v2 import setup as original_fit_setup
from covariance_exact_folds_v2 import certificate, variance_map
from full_covariance_qualification_sources import load as load_covariance, cohort_rows, jsonl, SUMMARY_FIELDS
from full_reduced_covariance_sources_v3 import SUMMARY_FIELDS as RETAINED_SUMMARY
from reference_measurement_union_sources import verify
from run_full_reduced_covariance_qualification_v3 import run as qualify_retained
from run_parallel_covariance_readback_v2 import STATUS as PARALLEL_STATUS


def nonempty_qualification_setup(root):
    # Certificates require positive membership. The declared complete software
    # grid contains all five nonempty cohorts, including too-small/rank-invalid
    # designs and constants. No production cohort is removed or redefined.
    pp=qualification_setup(root);plan=json.loads(pp.read_text())
    droot=Path(json.loads(Path(plan['design_plan']).read_text())['output'])
    cohorts=json.loads((droot/'cohort_manifest.json').read_text())
    empty={c['cohort_id'] for c in cohorts if c['records']==0};assert len(empty)==1
    cohorts=[c for c in cohorts if c['cohort_id'] not in empty]
    write(droot/'cohort_manifest.json',cohorts)
    for name in ['unique_designs.jsonl','unique_fit_inputs.jsonl']:
        rows=[r for r in jsonl(droot/name) if r['cohort_id'] not in empty]
        (droot/name).write_text(''.join(json.dumps(r)+'\n' for r in rows))
    with gzip.open(droot/'model_settings.tsv.gz','rt') as f:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames
        settings=[r for r in reader if r['cohort_id'] not in empty]
    with gzip.open(droot/'model_settings.tsv.gz','wt') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(settings)
    receipt=json.loads((droot/'receipt.json').read_text())
    receipt.update(model_setting_rows=len(settings),unique_cohorts=len(cohorts),unique_designs=30*len(cohorts))
    for name in ['cohort_manifest.json','unique_designs.jsonl','unique_fit_inputs.jsonl','model_settings.tsv.gz']:
        receipt['artifacts'][name]=sha(droot/name)
    write(droot/'receipt.json',receipt)
    summary={k:receipt[k] for k in ['logical_cases','model_setting_rows','unique_cohorts','unique_designs','trees']}
    close_fixture(droot,'complete_verified_full_expanded_model_designs',droot/'receipt.json',summary,{})
    plan['expected']['model_setting_rows']=len(settings);write(pp,plan)
    return pp


def setup(root):
    with patch('check_full_shared_entity_fits_v2.qualification_setup',side_effect=nonempty_qualification_setup):
        fit_path=original_fit_setup(root)
    fit_plan=json.loads(fit_path.read_text());qp=Path(fit_plan['qualification_plan']);qplan=json.loads(qp.read_text())
    source,bindings=load_covariance(qplan,qp)
    qroot=Path(qplan['output']);producer=json.loads((qroot/'receipt.json').read_text())
    summary={k:producer[k] for k in SUMMARY_FIELDS}
    parallel_root=root/'parallel';parallel_root.mkdir();(parallel_root/'failures').mkdir()
    parallel_plan=root/'parallel-plan.json'
    write(parallel_plan,dict(source_plan=str(qp),output=str(parallel_root),
        resources=dict(workers=8,maximum_pending_cohorts=16,execution_backend='fork_process')))
    reader=json.loads((qroot/'readback.json').read_text())
    reader.update(status=PARALLEL_STATUS,parallel_plan_sha256=sha(parallel_plan),
        parallel_workers=8,maximum_observed_pending_cohorts=4,execution_backend='fork_process',
        worker_pids=[1001,1002],failure_capture_enabled=True,original_numerical_arithmetic_unchanged=True,
        actual_cgroup_limits={'cpu.max':'800000 100000','memory.max':str(64*2**30),'memory.swap.max':'0'})
    independent=parallel_root/'readback.json';write(independent,reader)
    completed=close_fixture(qroot,'complete_verified_full_uniform_covariance_qualification',qroot/'receipt.json',summary,
        {str(qp):sha(qp),str(parallel_plan):sha(parallel_plan),str(independent):sha(independent)})
    compact=json.loads(completed.read_text());compact.update(independent_readback=str(independent),independent_readback_sha256=sha(independent));write(completed,compact)
    exact=root/'exact';exact.mkdir();exact_plan=root/'exact-plan.json';write(exact_plan,dict(output=str(exact)))
    certificates=[]
    for cohort in source['cohorts']:
        rows=cohort_rows(source,cohort)
        for mode in ['signed','unsigned']:
            operators={kind:source['operators'][mode,kind][rows].tocsr() for kind in
                ['target_node','background_node','model_pair','gene','model','family']}
            family=source['operators']['contrast','family_intercept'][rows].tocsr()
            value=certificate(operators,family,mode,rows)
            assert value['relations']['gene_equals_half_target_plus_background']['exact']
            assert value['relations']['model_equals_half_pair']['exact']
            certificates.append(dict(cohort_id=cohort['cohort_id'],certificate=value,
                variance_map=variance_map(value['relations'],mode),cohort_rows_sha256=cohort['case_rows_sha256'],
                ordered_case_ids_sha256=cohort['ordered_case_ids_sha256'],
                raw_reml_basis_qualification_complete=False,scientific_eligibility=False))
    cert_path=exact/'cohort_certificates.jsonl';cert_path.write_text(''.join(json.dumps(c)+'\n' for c in certificates))
    for name in ['receipt.json','readback.json']:write(exact/name,dict(explicit_synthetic_exact_closure_fixture=True))
    from check_full_reduced_covariance_qualification import closed
    exact_completed=closed(root,'exact','complete_verified_full_exact_uniform_covariance_folds_v2',
        [exact_plan,cert_path,exact/'receipt.json',exact/'readback.json'],
        dict(logical_cases=len(source['ids']),cohorts=len(source['cohorts']),certificates=len(certificates)))
    reduced_path=root/'retained-plan.json';reduced_root=root/'retained'
    expected={k:producer[k] for k in ['logical_cases','model_setting_rows','unique_cohorts','unique_designs','audit_rows','setting_audit_links']}
    reduced_completed=root/'retained'/'completion.json'
    reduced_plan=dict(qualification_plan=str(qp),qualification_completion=str(completed),parallel_readback_plan=str(parallel_plan),
        exact_plan=str(exact_plan),exact_completion=str(exact_completed),expected=expected,trees=fit_plan['trees'],
        output=str(reduced_root),pins={},resources=dict(minimum_free_disk_gib=0),completion=str(reduced_completed),
        scope='Declared complete five-cohort software grid. Closed parent/journal fixtures are synthetic; no biological acceptance.')
    write(reduced_path,reduced_plan)
    produced=qualify_retained(reduced_path);read=qualify_retained(reduced_path,reader=True)
    reduced_summary={k:produced[k] for k in RETAINED_SUMMARY}
    close_fixture(reduced_root,'complete_verified_full_exact_retained_uniform_covariance_qualification_v3',reduced_root/'receipt.json',
        reduced_summary,{str(reduced_root/'readback.json'):sha(reduced_root/'readback.json'),str(reduced_path):sha(reduced_path)})
    compact=json.loads(reduced_completed.read_text());compact.update(independent_readback=str(reduced_root/'readback.json'),
        independent_readback_sha256=sha(reduced_root/'readback.json'));write(reduced_completed,compact)
    fit_plan.update(qualification_completion=str(completed),retained_plan=str(reduced_path),retained_completion=str(reduced_completed),
        output=str(root/'output'),scope='Complete declared synthetic full source/grid checkpoint and link contracts; no production fitting or biological pilot.')
    write(fit_path,fit_plan)
    return fit_path
