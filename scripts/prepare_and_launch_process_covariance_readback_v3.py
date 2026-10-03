#!/usr/bin/env python3
"""Launch complete process-isolated source readback in a new immutable root."""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys

import psutil

from ancestral_chain_attempt import sha
from full_covariance_qualification_sources import SUMMARY_FIELDS
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from record_project_runtime_checkpoint_v4 import journal_terminal
from reference_measurement_union_sources import verify
from run_parallel_covariance_readback_v2 import STATUS


def main():
    gp=Path('metadata/process_covariance_readback_software_validation_20261003_v2.json');gate=json.loads(gp.read_text())
    assert gate['status']=='passed_complete_process_covariance_readback_software_contracts_v2'
    assert (gate['fixture_cohorts'],gate['fixture_audits'],gate['fixture_setting_links'])==(6,1800,7200)
    assert gate['worker_counts']==[1,8] and gate['serial_numerical_and_summary_agreement'] is True
    assert len(gate['rehashed_source_cases_rejected'])==17 and len(gate['serialized_result_cases_rejected'])==4
    assert gate['execution_backend']=='fork_process' and gate['worker_process_identity_checked'] is True
    assert gate['exact_rejected_audit_matrix_reference_capture_checked'] is True
    verify(gate['source_hashes'])
    transport_path=Path('metadata/process_covariance_readback_software_transport_20261003_v2.json')
    transport=json.loads(transport_path.read_text());verify(transport['source_hashes'])
    assert transport['actual_tool_terminal_exit_code']==0 and transport['tool_session_id']==98973
    unit='fungal-process-covariance-readback-software-20261003-v2.service';inv=transport['invocation_id']
    raw=subprocess.check_output(['journalctl','--user','-u',unit,'-o','json','--no-pager'],text=True)
    jr=[json.loads(l) for l in raw.splitlines()]
    execution=json.loads(Path('metadata/process_covariance_readback_software_execution_20261003_v2.json').read_text())
    assert execution['status']=='exited_zero_with_receipt' and execution['exit_code']==0
    assert execution['invocation_id']==inv and execution['receipt_sha256']==sha(gp)
    verify(execution['source_hashes']);verify(execution['artifacts'])
    command=execution['wrapper']['cmdline']
    # Pipe-mode software stdout is captured by the actual tool. Its original
    # wrapper identity/configuration and invocation-specific manager command,
    # time and completion are checked; absent wrapper journal rows are explicit.
    starts=[r for r in jr if r.get('USER_INVOCATION_ID')==inv and r.get('MESSAGE','').startswith('Started ')
        and ' '.join(command) in r['MESSAGE']]
    assert len(starts)==1
    assert abs(float(starts[0]['__REALTIME_TIMESTAMP'])/1e6-execution['wrapper']['created'])<2
    resource=[r for r in jr if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')]
    assert resource and not any('Failed with result' in r.get('MESSAGE','') or 'Main process exited' in r.get('MESSAGE','') for r in jr if r.get('USER_INVOCATION_ID')==inv)
    assert transport['wrapper']==execution['wrapper'] and transport['original_command_start_records']==1
    ep='metadata/process_covariance_readback_software_execution_20261003_v2.json'
    sp='metadata/full_uniform_covariance_qualification_plan_20261002.json';source=json.loads(Path(sp).read_text())
    rp=Path(source['output'])/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_full_uniform_covariance_qualification_pending_independent_readback'
    assert (receipt['logical_cases'],receipt['model_setting_rows'],receipt['unique_cohorts'],receipt['unique_designs'],
        receipt['audit_rows'],receipt['setting_audit_links'])==(75188,622080,4340,130200,1302000,6220800)
    assert receipt['scientific_eligibility'] is False and receipt['plan_sha256']==sha(sp)
    launch_inventory=json.loads(Path('metadata/full_uniform_covariance_qualification_launches_20261002.json').read_text())
    producer_launch=launch_inventory['launches'][0];original=json.loads(Path(producer_launch).read_text());original['launch']=producer_launch
    terminal=journal_terminal(original)
    assert terminal['status']=='verified_original_terminal_success_with_bound_completed_artifacts'
    output='results/phylogeny/full-uniform-covariance-process-readback-20261003-v2';assert not Path(output).exists()
    assert psutil.virtual_memory().available>=64*2**30 and shutil.disk_usage('.').free>=128*2**30
    resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(),cpus=8,memory_gib=64,swap_gib=0,workers=8,
        maximum_pending_cohorts=16,blas_threads=1,address_space_gib_per_process=48,cpu_seconds_per_process=65536,maximum_numeric_processes=9,
        maximum_sum_numeric_cpu_limits_seconds=589824,per_file_limit_gib=2,execution_backend='fork_process',
        planning_output_gib=16,minimum_free_disk_gib=128,source_audit_compressed_bytes=(rp.parent/'design_covariance_audits.jsonl.gz').stat().st_size,
        source_links_compressed_bytes=(rp.parent/'setting_audit_links.tsv.gz').stat().st_size,
        existing_serial_reader_sampled_rss_gib=3.1841659545898438,existing_serial_reader_recent_completed_cohorts=714,
        finish_eta=None,runtime_calibrated=False,new_cost_usd=0,gpu=False,
        available_memory_gib=psutil.virtual_memory().available/2**30,available_disk_gib=shutil.disk_usage('.').free/2**30,
        caveat='Eight forked workers inherit the verified source, with copy-on-write overhead, rather than receiving full source copies through IPC. At most16cohort payloads pending. The64GiB0swap cgroup bounds aggregate memory;48GiB AS and65536CPU seconds are per numeric process, not aggregate. At most9 numeric processes gives589824sum CPU allowance. Planning estimate: parent source~3.2GiB; at most~26GiB inherited-private growth plus~32GiB eight-worker arrays/temporary products, under64GiB but not measured peaks. Existing source/RSS and targeted2700replay6GiB peak are observations, not full-grid bounds. No calibrated ETA.16GiB output planning includes cohort reports and possible rejected-case matrices, not a disk quota. Source hash verification finishes before worker creation; failures propagate and cannot close. Scheduling is not an established repair of the original mismatch.')
    own=['parallel_covariance_readback_v2','run_parallel_covariance_readback_v2','check_parallel_covariance_readback_v2',
        'prepare_and_launch_process_covariance_readback_v3','full_covariance_qualification_sources','readback_full_covariance_qualification',
        'covariance_basis_independent','covariance_basis_context','covariance_basis_audit','full_entity_operator_sources',
        'full_expanded_model_design_sources','full_expanded_model_input_sources','background_measurement_union_sources',
        'reference_measurement_union_sources','run_after_verified_dependencies_v2','close_full_triad_sequence_stage',
        'record_completed_process_handoffs_v2']
    pins={str(p):sha(p) for p in [*[Path('scripts/'+n+'.py') for n in own],Path(sp),rp,gp,Path(ep),
        Path('metadata/process_covariance_readback_software_resources_20261003_v2.json'),transport_path,Path(producer_launch),Path('/usr/bin/prlimit'),Path(sys.executable)]}
    scope=('Complete independent readback of the unchanged original75188cases/622080settings/4340cohorts/130200designs/'
        '1302000audits/6220800setting links, both loading modes and all five trees. Eight bounded forked cohort processes change only '
        'scheduling; original latent independent products, normalizedQR, reciprocal-condition checks, entire Gram/error '
        'envelopes,gesvd rank validation and failure reproduction remain frozen. All source hashes, recipes, review/empty/'
        'constant/failed states and exact SQL identities retained. At most16pending cohort payloads; new immutable root, '
        'capture rejected exact audits/matrices/reference products without allowing completion; no original edits/restarts/output substitution. Original producer plus this complete independent reader require '
        'both exact original process/invocation journals and complete hashes for separate alternate closure. Original serial '
        'reader and queued descendants remain unchanged. New downstream consumers need a separately versioned explicit '
        'alternate-readback source adapter; old paths are never written. No numerical tolerance relaxation, rank deletion, '
        'accepted covariance/variance model, production optimizer or biological effect. Existing CPU resources only; '
        'no GPU, new charges or completed biological aims.')
    pp='metadata/process_covariance_readback_plan_20261003_v2.json'
    plan=dict(source_plan=sp,output=output,pins=pins,resources=resources,source_producer_launch=producer_launch,
        software_validation=str(gp),completion='metadata/full_uniform_covariance_process_completed_20261003_v2.json',
        launch_inventory='metadata/process_covariance_readback_launches_20261003_v2.json',scope=scope)
    create(pp,plan)
    prefix=['/usr/bin/prlimit','--as='+str(48*2**30),'--cpu=65536','--fsize='+str(2*2**30),'--']
    reader=launch('process-covariance-readback-v2',prefix+[sys.executable,'scripts/run_parallel_covariance_readback_v2.py','--plan',pp],
        [producer_launch],pp,cpus=8,memory=64)
    cp='metadata/process_covariance_readback_completion_plan_20261003_v2.json'
    create(cp,dict(source_plan=sp,producer_receipt=str(rp),independent_readback=output+'/readback.json',
        producer_status='complete_full_uniform_covariance_qualification_pending_independent_readback',reader_status=STATUS,
        completed_status='complete_verified_full_uniform_covariance_qualification',summary_fields=SUMMARY_FIELDS,
        launches=[producer_launch,reader],pins={p:sha(p) for p in [sp,pp,producer_launch,reader,
            'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'],scope=scope))
    closer=launch('process-covariance-readback-v2-closure',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],
        [producer_launch,reader],cp,cpus=2,memory=32)
    create(plan['launch_inventory'],dict(status='launched_full_process_isolated_covariance_readback_v2',source_plan=pp,
        source_plan_sha256=sha(pp),launches=[reader,closer],original_producer_launch=producer_launch,
        original_producer_terminal_observation=terminal,full_cohorts=4340,full_audits=1302000,full_setting_links=6220800,
        original_serial_jobs_changed=False,original_readback_path_replaced=False,production_fitting_launched=False,gpu=False,new_cost_usd=0))
    print(json.dumps(resources),flush=True)


if __name__=='__main__':main()
