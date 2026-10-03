#!/usr/bin/env python3
"""Serial/parallel agreement and rehashed malformed full synthetic exports."""
import argparse
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import resource
import time
from unittest.mock import patch

from ancestral_chain_attempt import sha
from check_full_covariance_qualification import setup, mutate
from full_covariance_qualification_sources import SUMMARY_FIELDS
from prepare_full_covariance_qualification import run as produce
from readback_full_covariance_qualification import run as serial
from reference_measurement_union_sources import verify
import run_parallel_covariance_readback as workflow


def rejected(action):
    try:action()
    except (AssertionError,ValueError,KeyError,RuntimeError,FileNotFoundError,FileExistsError,StopIteration):return
    raise AssertionError('Altered parallel readback accepted')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args();root=a.output.resolve();root.mkdir(exist_ok=False)
    assert not a.receipt.exists();started=time.monotonic()
    source=root/'fixture';source.mkdir();sp=setup(source);producer=produce(sp)
    original=serial(sp,root/'serial-reference.json')
    plans=[];results=[]
    def plan(name,workers):
        pp=root/(name+'-plan.json');value=dict(source_plan=str(sp),output=str(root/name),pins={},
            resources=dict(cpus=2,memory_gib=8,workers=workers,maximum_pending_cohorts=2*workers,minimum_free_disk_gib=1),
            scope='Synthetic complete source closures/journals and numerical recipes only; no biological samples or production acceptance.')
        pp.write_text(json.dumps(value,indent=2)+'\n');plans.append(pp);return pp
    for workers in [1,8]:
        pp=plan('valid-'+str(workers),workers);r=workflow.run(pp);results.append(r)
        assert all(r[k]==original[k] for k in SUMMARY_FIELDS)
        assert r['maximum_absolute_projected_gram_error']==original['maximum_absolute_projected_gram_error']
        assert r['conservative_independent_classification_differences']==original['conservative_independent_classification_differences']
        assert (r['unique_cohorts'],r['audit_rows'],r['setting_audit_links'])==(6,1800,7200)
        assert r['maximum_observed_pending_cohorts']<=2*workers
        rejected(lambda:workflow.run(pp))
    snapshot={p.name:p.read_bytes() for p in (source/'output').iterdir() if p.is_file()}
    changes=['omit_final_tree','omit_unsigned','wrong_raw_gram','wrong_projected_gram','wrong_envelope','wrong_rank',
        'promote_dependent_basis','undo_target_fold','wrong_family_fold','wrong_design_hash','wrong_source_contract','false_failure',
        'link_omit_final_tree','link_omit_empty','link_wrong_fit','link_wrong_audit','link_wrong_status']
    failed=[]
    for name in changes:
        mutate(source,name);pp=plan('rejected-'+name,8)
        rejected(lambda:workflow.run(pp));failed.append(name)
        for name,raw in snapshot.items():(source/'output'/name).write_bytes(raw)
    serialization=[];real=workflow.verify_cohort
    for name in ['missing_index','changed_index_status','changed_design_batch_hash','promoted_scientific_flag']:
        def altered(*args):
            value=real(*args)
            if name=='missing_index':value['audit_index'].pop()
            elif name=='changed_index_status':value['audit_index'][0][-1]='invented'
            elif name=='changed_design_batch_hash':value['source_design_batch_sha256']='0'*64
            else:value['scientific_eligibility']=True
            return value
        pp=plan('rejected-serialized-'+name,8)
        with patch.object(workflow,'verify_cohort',side_effect=altered):rejected(lambda:workflow.run(pp))
        serialization.append(name)
    bindings={str(p):sha(p) for p in [Path(__file__),Path('scripts/parallel_covariance_readback.py'),
        Path('scripts/run_parallel_covariance_readback.py'),Path('scripts/readback_full_covariance_qualification.py'),sp,*plans]}
    for p in root.rglob('*'):
        if p.is_file():bindings[str(p)]=sha(p)
    verify(bindings)
    result=dict(status='passed_complete_parallel_covariance_readback_software_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        fixture_cohorts=6,fixture_audits=1800,fixture_setting_links=7200,trees=5,loading_modes=2,
        worker_counts=[1,8],serial_numerical_and_summary_agreement=True,all_source_review_empty_constant_states_retained=True,
        completed_restart_refused=True,partial_roots_preserved=True,rehashed_source_cases_rejected=failed,
        serialized_result_cases_rejected=serialization,source_and_journal_fixtures_synthetic=True,
        numerical_arithmetic_and_original_tolerances_unchanged=True,scientific_eligibility=False,
        elapsed_seconds=time.monotonic()-started,self_cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+
            resource.getrusage(resource.RUSAGE_SELF).ru_stime,self_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        source_hashes=bindings,scope='Complete synthetic original grid through actual unchanged latent/QR/gesvd calculations and one/eight-worker scheduling, all exports/SQL identities. Synthetic parent journals are not production proof.21 altered source/result cases and completed restarts rejected; no fungal pilot, source modification or biological inference.')
    a.receipt.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2),flush=True)


if __name__=='__main__':main()
