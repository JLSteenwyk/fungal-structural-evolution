#!/usr/bin/env python3
"""Launch a complete parallel source readback without editing original jobs."""
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
from run_parallel_covariance_readback import STATUS


def main():
    gp=Path('metadata/parallel_covariance_readback_software_validation_20261003_v1.json');gate=json.loads(gp.read_text())
    assert gate['status']=='passed_complete_parallel_covariance_readback_software_contracts'
    assert (gate['fixture_cohorts'],gate['fixture_audits'],gate['fixture_setting_links'])==(6,1800,7200)
    assert gate['worker_counts']==[1,8] and gate['serial_numerical_and_summary_agreement'] is True
    assert len(gate['rehashed_source_cases_rejected'])==17 and len(gate['serialized_result_cases_rejected'])==4
    verify(gate['source_hashes'])
    unit='fungal-parallel-covariance-readback-software-20261003-v1.service';inv='c866f80316444ff29d3b977576c01885'
    raw=subprocess.check_output(['journalctl','--user','-u',unit,'-o','json','--no-pager'],text=True)
    jr=[json.loads(l) for l in raw.splitlines()]
    command=[sys.executable,'scripts/check_parallel_covariance_readback.py','--output',
        'data/software_audits/parallel-covariance-readback-20261003-v1','--receipt',str(gp)]
    exact=[r for r in jr if r.get('_CMDLINE')==' '.join(command) and r.get('_SYSTEMD_INVOCATION_ID')==inv]
    assert exact and len({r['_PID'] for r in exact})==1
    resource=[r for r in jr if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')]
    assert resource and not any('Failed with result' in r.get('MESSAGE','') or 'Main process exited' in r.get('MESSAGE','') for r in jr if r.get('USER_INVOCATION_ID')==inv)
    journal=Path('data/software_audits/parallel-covariance-readback-20261003-v1/original-invocation-journal.jsonl')
    assert not journal.exists();journal.write_text(raw)
    ep='metadata/parallel_covariance_readback_software_execution_20261003_v1.json'
    create(ep,dict(status='verified_original_parallel_covariance_software_wait_exited_zero',actual_tool_terminal_exit_code=0,
        tool_session_id=58029,invocation_id=inv,original_journal_pid=int(exact[0]['_PID']),original_command=command,
        original_process_messages=len(exact),original_completion_resource_records=len(resource),
        software_self_cpu_seconds=gate['self_cpu_seconds'],software_self_peak_rss_bytes=gate['self_peak_rss_bytes'],
        source_hashes={str(gp):sha(gp),str(journal):sha(journal)},
        scope='Observed original tool wait, exact invocation/PID/command journal and original completion resources. Self RSS is not charged group/native peak. Source/journal fixtures inside numerical software runs remain synthetic.'))
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
    output='results/phylogeny/full-uniform-covariance-parallel-readback-20261003-v1';assert not Path(output).exists()
    assert psutil.virtual_memory().available>=64*2**30 and shutil.disk_usage('.').free>=128*2**30
    resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(),cpus=8,memory_gib=64,swap_gib=0,workers=8,
        maximum_pending_cohorts=16,blas_threads=1,address_space_gib=48,cpu_seconds=604800,per_file_limit_gib=2,
        planning_output_gib=16,minimum_free_disk_gib=128,source_audit_compressed_bytes=(rp.parent/'design_covariance_audits.jsonl.gz').stat().st_size,
        source_links_compressed_bytes=(rp.parent/'setting_audit_links.tsv.gz').stat().st_size,
        existing_serial_reader_sampled_rss_gib=3.1841659545898438,existing_serial_reader_recent_completed_cohorts=714,
        finish_eta=None,runtime_calibrated=False,new_cost_usd=0,gpu=False,
        available_memory_gib=psutil.virtual_memory().available/2**30,available_disk_gib=shutil.disk_usage('.').free/2**30,
        caveat='Eight cohort workers share read-only source arrays; at most16cohort payloads are pending.64GiB/48GiB limits bound whole reader, not individual thread peaks. Existing serial RSS is a checkpoint observation, not a complete peak or reliable scaling model.604800CPU seconds is an upper allowance, not an ETA.16GiB covers new reports/source receipts/closure planning, not a disk quota. All original source hashes and numerical checks retained.')
    own=['parallel_covariance_readback','run_parallel_covariance_readback','check_parallel_covariance_readback',
        'prepare_and_launch_parallel_covariance_readback','full_covariance_qualification_sources','readback_full_covariance_qualification',
        'covariance_basis_independent','covariance_basis_context','covariance_basis_audit','full_entity_operator_sources',
        'full_expanded_model_design_sources','full_expanded_model_input_sources','background_measurement_union_sources',
        'reference_measurement_union_sources','run_after_verified_dependencies_v2','close_full_triad_sequence_stage',
        'record_completed_process_handoffs_v2']
    pins={str(p):sha(p) for p in [*[Path('scripts/'+n+'.py') for n in own],Path(sp),rp,gp,Path(ep),
        Path('metadata/parallel_covariance_readback_software_resources_20261003_v1.json'),Path(producer_launch),Path('/usr/bin/prlimit'),Path(sys.executable)]}
    scope=('Complete independent readback of the unchanged original75188cases/622080settings/4340cohorts/130200designs/'
        '1302000audits/6220800setting links, both loading modes and all five trees. Eight bounded cohort threads change only '
        'scheduling; original latent independent products, normalizedQR, reciprocal-condition checks, entire Gram/error '
        'envelopes,gesvd rank validation and failure reproduction remain frozen. All source hashes, recipes, review/empty/'
        'constant/failed states and exact SQL identities retained. At most16pending cohort payloads; new immutable root, '
        'no original edits/restarts/output substitution. Original producer plus this complete independent reader require '
        'both exact original process/invocation journals and complete hashes for separate alternate closure. Original serial '
        'reader and queued descendants remain unchanged. New downstream consumers need a separately versioned explicit '
        'alternate-readback source adapter; old paths are never written. No numerical tolerance relaxation, rank deletion, '
        'accepted covariance/variance model, production optimizer or biological effect. Existing CPU resources only; '
        'no GPU, new charges or completed biological aims.')
    pp='metadata/parallel_covariance_readback_plan_20261003_v1.json'
    plan=dict(source_plan=sp,output=output,pins=pins,resources=resources,source_producer_launch=producer_launch,
        software_validation=str(gp),completion='metadata/full_uniform_covariance_parallel_completed_20261003_v1.json',
        launch_inventory='metadata/parallel_covariance_readback_launches_20261003_v1.json',scope=scope)
    create(pp,plan)
    prefix=['/usr/bin/prlimit','--as='+str(48*2**30),'--cpu=604800','--fsize='+str(2*2**30),'--']
    reader=launch('parallel-covariance-readback-v1',prefix+[sys.executable,'scripts/run_parallel_covariance_readback.py','--plan',pp],
        [producer_launch],pp,cpus=8,memory=64)
    cp='metadata/parallel_covariance_readback_completion_plan_20261003_v1.json'
    create(cp,dict(source_plan=sp,producer_receipt=str(rp),independent_readback=output+'/readback.json',
        producer_status='complete_full_uniform_covariance_qualification_pending_independent_readback',reader_status=STATUS,
        completed_status='complete_verified_full_uniform_covariance_qualification',summary_fields=SUMMARY_FIELDS,
        launches=[producer_launch,reader],pins={p:sha(p) for p in [sp,pp,producer_launch,reader,
            'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'],scope=scope))
    closer=launch('parallel-covariance-readback-v1-closure',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],
        [producer_launch,reader],cp,cpus=2,memory=32)
    create(plan['launch_inventory'],dict(status='launched_full_alternate_parallel_covariance_readback',source_plan=pp,
        source_plan_sha256=sha(pp),launches=[reader,closer],original_producer_launch=producer_launch,
        original_producer_terminal_observation=terminal,full_cohorts=4340,full_audits=1302000,full_setting_links=6220800,
        original_serial_jobs_changed=False,original_readback_path_replaced=False,production_fitting_launched=False,gpu=False,new_cost_usd=0))
    print(json.dumps(resources),flush=True)


if __name__=='__main__':main()
