#!/usr/bin/env python3
"""Prepare the complete failed-role matrix only after actual first-input/native/readback closure."""
from datetime import datetime,timezone
import json
from pathlib import Path

import psutil

from ancestral_chain_attempt import sha
from baliphy_joint_fasta_v7_failure_jobs_v1 import failed_sources,validate_jobs
from reference_measurement_union_sources import bind,verify
from run_baliphy_joint_fasta_v7_failure_grid_v1 import gates


def save(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def main():
    pp=Path('metadata/baliphy_joint_fasta_v7_failure_grid_plan_20261004_v1.json')
    rp=Path('metadata/baliphy_joint_fasta_v7_failure_grid_resources_20261004_v1.json')
    assert not pp.exists() and not rp.exists()
    descriptors=[]
    for prefix,version,status in [
        ('baliphy_joint_fasta_v7_failure_jobs_software','v1','passed_exact_all24_failed_role_candidate_job_contracts_v1'),
        ('baliphy_joint_fasta_v7_full_comparison_readback','v1','passed_independent_full_v7_candidate_comparison_readback')]:
        receipt='metadata/'+prefix+('_validation' if 'software' in prefix else '')+'_20261004_'+version+'.json'
        transport='metadata/'+prefix+'_transport_20261004_'+version+'.json'
        descriptors.append(dict(receipt=receipt,transport=transport,status=status))
    pins=gates(dict(gates=descriptors));verify(pins)
    actual=json.loads(Path(descriptors[1]['receipt']).read_text())
    assert actual['scalar_rows']==21 and actual['mapped_scalar_values']==903 and actual['joint_frames']==3
    assert actual['native_tips']==622 and actual['native_ancestors']==621
    qualified=json.loads(Path(descriptors[0]['receipt']).read_text())
    jobs_path=Path(qualified['jobs']);jobs=json.loads(jobs_path.read_text())
    originals,source_pins,mapping=failed_sources();scope=validate_jobs(jobs,originals)
    for path,h in source_pins.items():bind(pins,path,h)
    root='results/ancestral/baliphy-joint-fasta-v7-all-failed-roles-20261004-v1'
    assert not Path(root).exists()
    native_cpu_caps=[int(next(v.split('=')[1] for v in j['config']['command'] if v.startswith('--cpu='))) for j in jobs]
    native_wall_caps=[j['config']['timeout_seconds'] for j in jobs]
    completed=json.loads(Path('metadata/baliphy_joint_fasta_v7_full_comparison_review_20261004_v2.json').read_text())
    observed=completed['native_elapsed_seconds']
    resources=dict(prepared_utc=datetime.now(timezone.utc).isoformat(),cpus=4,workers=4,memory_gib=200,
        reservation_capacity_gib=192,controller_soft_address_space_gib=8,address_space_gib=192,
        swap_gib=0,blas_threads=1,cpu_seconds_per_stage=43200,wall_seconds_per_stage=86400,
        per_file_limit_mib=2048,minimum_available_ram_gib=200,minimum_free_disk_gib=256,
        available_ram_gib=psutil.virtual_memory().available/2**30,free_disk_gib=psutil.disk_usage('.').free/2**30,
        original_native_as_gib=48,original_stack_soft_mib=8,original_stack_hard_unlimited=True,
        sum_original_native_cpu_cap_seconds=sum(native_cpu_caps),sum_original_role_wall_timeout_seconds=sum(native_wall_caps),
        maximum_original_native_cpu_seconds=max(native_cpu_caps),maximum_original_role_wall_seconds=max(native_wall_caps),
        actual_first_candidate_native_wall_seconds=observed,
        one_case_linear_worker_hours_for24=24*observed/3600,
        one_case_linear_four_worker_wall_hours=6*observed/3600,
        prospective_wall_hours=[2,12],runtime_uncalibrated_beyond_one_case=True,
        output_allowance_gib=128,new_cost_usd=0,gpu=False,
        scope='All24failed-role comparison,4workers,original48GiB/native CPU/wall/file/8MiB stack limits. '
              'Four48GiB AS leases reserve192GiB, with8GiB controller headroom under200GiB/noSwap. '
              'Parent8GiB soft/192GiB hard allows native prlimit to lower hard AS to48GiB. '
              'One-case24minute linear estimate is planning only, not measured grid throughput or posterior adequacy.')
    assert resources['available_ram_gib']>=200 and resources['free_disk_gib']>=256
    save(rp,resources)
    for path in [Path(__file__),Path('scripts/run_baliphy_joint_fasta_v7_failure_grid_v1.py'),
                 Path('scripts/baliphy_joint_fasta_v7_failure_jobs_v1.py'),
                 Path('scripts/run_baliphy_scalar_v6_sampler_v1.py'),Path('scripts/reference_sampler_memory_budget.py'),
                 jobs_path,Path(mapping),rp,Path('metadata/baliphy_joint_fasta_v7_full_comparison_review_20261004_v2.json')]:bind(pins,path)
    verify(pins)
    plan=dict(status='prepared_full_all24_original_failure_comparison_after_actual_native_and_reader_closure',
        prepared_utc=datetime.now(timezone.utc).isoformat(),jobs=str(jobs_path.resolve()),mapping=mapping,
        output=root,job_scope=scope,gates=descriptors,resources=resources,pins=pins,
        scope='All24original nativeSIGSEGV roles,2full622-tip effective inputs,all3priors and4chains '
              'per quartet,under the qualified pure row-wise logger candidate. Each old native '
              'seed,CPU/wall/AS/file cap and original8MiB stack is retained in a fresh output root. '
              'Native failure/nonfinite/invalid outcomes stay explicit; no automatic retry or dropped '
              'role. All successes require the existing independent source/tree/scalar/joint readers '
              'and serialized replay. This compares a corrected logger on the full failure set; '
              'it does not rerun every other1596role or establish adequate ancestral posterior, '
              'fullV7grid acceptance, biological inference, predictor calibration or any completed aim.')
    save(pp,plan)
    print(json.dumps({k:v for k,v in plan.items() if k not in ['pins']},indent=2))


if __name__=='__main__':main()
