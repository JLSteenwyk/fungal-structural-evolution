#!/usr/bin/env python3
"""Fresh whole-stage hashes, serialized census accounting and original terminals."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from full_covariance_qualification_sources import jsonl
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from reference_measurement_union_sources import verify


def main():
    p = argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    cp=Path('metadata/full_weighted_covariance_source_census_completed_20261004_v1.json');c=json.loads(cp.read_text())
    status='complete_verified_full_four_control_covariance_source_census_v1'
    assert c['status']==status and c['scientific_eligibility'] is False and c['exact_process_journals_checked']==2
    ap=Path(c['full_hash_archive']);assert sha(ap)==c['full_hash_archive_sha256'];archive=json.loads(ap.read_text())
    assert archive['status']==status+'_archive' and len(archive['services'])==2
    assert len(archive['source_hashes'])==c['bound_source_hashes'];verify(archive['source_hashes'])
    assert all(c[k]==v for k,v in archive['summary'].items())
    producer=json.loads(Path(c['producer_receipt']).read_text());reader=json.loads(Path(c['independent_readback']).read_text())
    assert sha(c['producer_receipt'])==c['producer_receipt_sha256']==reader['producer_receipt_sha256']
    assert sha(c['independent_readback'])==c['independent_readback_sha256']
    assert producer['source_contract']==reader['source_contract']
    for r in [producer,reader]:
        assert r['scientific_eligibility'] is r['raw_reml_basis_qualification_complete'] is False
        assert r['numerical_audits_computed']==r['working_model_fits_computed']==0
        assert all(r[k]==v for k,v in archive['summary'].items())
        assert all(archive['source_hashes'][k]==h for k,h in r['source_hashes'].items())
    root=Path(c['producer_receipt']).parent
    assert set(producer['artifacts'])=={'stage_plan.json','cohort_source_census.jsonl.gz'}
    assert producer['artifacts']==reader['artifacts']
    for n,h in producer['artifacts'].items():assert archive['source_hashes'][str(root/n)]==h
    pp=Path('metadata/full_weighted_covariance_source_census_plan_20261004_v1.json');plan=json.loads(pp.read_text())
    assert producer['plan_sha256']==reader['plan_sha256']==sha(pp);verify(plan['pins'])
    designroot=Path(json.loads(Path(plan['design_plan']).read_text())['output'])
    cohorts=json.loads((designroot/'cohort_manifest.json').read_text());rows_n=count=0;design_ids=set();fit_ids=set();policies=Counter()
    exports=iter(jsonl(root/'cohort_source_census.jsonl.gz'))
    for cohort in cohorts:
        row=next(exports);assert row['cohort_id']==cohort['cohort_id'] and row['records']==cohort['records']
        assert row['cohort_rows_sha256']==cohort['case_rows_sha256'] and row['ordered_case_ids_sha256']==cohort['ordered_case_ids_sha256']
        assert row['guide']==cohort['guide'] and row['mask']==cohort['mask']
        assert len(row['design_ids'])==len(set(row['design_ids']))==30
        assert len(row['fit_input_ids'])==len(set(row['fit_input_ids']))==60
        assert not design_ids.intersection(row['design_ids']) and not fit_ids.intersection(row['fit_input_ids'])
        design_ids.update(row['design_ids']);fit_ids.update(row['fit_input_ids'])
        choices=row['policy_basis_choices'];assert len(choices)==8
        assert {(r['policy'],r['loading_mode']) for r in choices}=={(p,m) for p in c['policies'] for m in c['loading_modes']}
        for r in choices:
            assert r['route']==('closed_exact_uniform_named_fold' if r['exact_uniform_one'] else 'closed_positive_diagonal_cone')
            assert ('target_node' not in r['names'])==r['exact_uniform_one']
            policies[r['policy']+':'+r['loading_mode']+':q'+str(len(r['names']))]+=1
        assert row['prospective_numerical_audits']==1200
        assert row['raw_reml_basis_qualification_complete'] is row['scientific_eligibility'] is False
        count+=1;rows_n+=row['records']
    assert next(exports,None) is None
    assert (count,len(design_ids),len(fit_ids),rows_n)==(4340,130200,260400,34110120)
    assert dict(policies)==c['policy_basis_counts']
    inventory=json.loads(Path('metadata/full_weighted_covariance_source_census_launches_20261004_v1.json').read_text())
    terminals=[]
    for lp in inventory['launches']:
        r=json.loads(Path(lp).read_text());r['launch']=lp
        assert fingerprint(r) is None;terminals.append(journal_terminal(r))
    result=dict(status='verified_completed_full_four_control_source_census_whole_stage_hashes_and_original_terminals_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),completion=str(cp),completion_sha256=sha(cp),
        full_hash_archive_sha256=sha(ap),fresh_stage_bindings_checked=len(archive['source_hashes']),
        original_terminal_handles_checked=len(terminals),original_terminals=terminals,
        logical_cases=75188,cohorts=count,case_row_occurrences=rows_n,designs=len(design_ids),fit_inputs=len(fit_ids),
        original_settings=c['settings'],policy_basis_counts=dict(policies),numerical_audits_computed=0,
        working_model_fits_computed=0,scientific_eligibility=False,
        scope='Every actual closed stage binding freshly rehashed; complete serialized cohort/design/response/basis census and three exact original terminal handles verified. Complete X/y/setting reconstruction rests on the closed source producer/reader; not rerun here. No broader unconsumed parent replay, covariance numeric audit, fit or biological acceptance.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='original_terminals'}),flush=True)


if __name__=='__main__':main()
