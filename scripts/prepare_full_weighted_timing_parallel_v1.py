#!/usr/bin/env python3
"""Freeze full parallel timing and fitting-source drafts without bypassing closure."""
from datetime import datetime,timezone
import json
from pathlib import Path
import sys

import psutil

from ancestral_chain_attempt import sha
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind,verify


def write(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2);f.write('\n')


def main():
    pins={}
    for stem,session,status in [
        ('full_weighted_timing_parallel',5573,'passed_complete_parallel_four_control_timing_and_numeric_readback_v1'),
        ('weighted_timing_parallel_memory_guards',56888,'passed_private_timing_source_and_submitted_task_memory_guards_v1')]:
        paths=[Path('metadata/'+stem+'_software_'+kind+'_20261004_v1.json') for kind in ['validation','execution','transport']]
        v,e,t=[json.loads(path.read_text()) for path in paths]
        assert v['status']==status and e['status']=='exited_zero_with_receipt' and e['exit_code']==0 and not e['timed_out']
        assert sha(paths[0])==e['receipt_sha256'] and t['original_tool_session_id']==session
        assert t['original_tool_terminal_exit_code']==0 and t['wrapper']==e['wrapper'] and t['invocation_id']==e['invocation_id']
        assert t['whole_wrapper_initial_and_terminal_payloads_matched'] and t['manager_start_records']==t['manager_completion_records']==1
        for record in [v,e,t]:verify(record['source_hashes']);verify(record.get('artifacts',{}))
        for path in paths:bind(pins,path)
        if stem=='full_weighted_timing_parallel':
            assert v['total_candidate_rows']==72000 and v['full_grid_selected_groups']==320
            assert v['serial_parallel_numeric_probe_groups_compared']==v['independent_numeric_groups_replayed']==320
            assert v['actual_variance_points_per_group']==2 and len(v['private_rehashed_corruptions_rejected'])==17
            assert not v['full_grid_native_probes_mocked'] and v['fits_computed']==0
        else:
            assert [(case['case'],case['rejected']) for case in v['cases']]==[('control',False),('source_before',True),('source_after',True),('task_after',True)]
            assert v['parent_source_and_submitted_task_arrays_preserved']
    oldfitpath=Path('metadata/full_weighted_shared_entity_fit_draft_plan_20261004_v1.json')
    oldfit=json.loads(oldfitpath.read_text());verify(oldfit['pins']);bind(pins,oldfitpath)
    qp=Path('metadata/full_weighted_parallel_qualification_plan_20261004_v1.json')
    q=json.loads(qp.read_text());verify(q['pins']);bind(pins,qp)
    assert q['expected']['cohorts']==4340 and q['expected']['designs']==130200
    assert oldfit['expected']['candidate_rows']==20832000 and oldfit['expected']['setting_fit_links']==49766400
    assert oldfit['optimizer']['gradient_tolerance']==1e-6 and oldfit['independent_audit']['column_batch']==32
    assert oldfit['launch_state']=='not_launched_or_queued' and not Path(oldfit['output']).exists()
    modules=project_sources(pins,[Path(__file__),Path('scripts/full_weighted_timing_parallel_v1.py'),
        Path('scripts/close_full_weighted_timing_parallel_v1.py'),Path('scripts/record_full_weighted_parallel_launch_v1.py')])
    fp=Path('metadata/full_weighted_parallel_source_fit_draft_plan_20261004_v1.json')
    fit=dict(oldfit,prepared_utc=datetime.now(timezone.utc).isoformat(),qualification_plan=str(qp),
        qualification_completion=q['completion'],output='results/phylogeny/full-parallel-source-four-control-shared-entity-fits-20261004-v1',
        pins={**oldfit['pins'],**pins},source_adapter_module='full_weighted_shared_entity_fit_sources_parallel_v1',
        full_timing_implemented=True,full_timing_complete=False,
        fitting_export_runner_integration_pending=True,launch_state='not_launched_or_queued',
        status='parallel_checkpoint_source_qualified_full_fit_draft_pending_real_numeric_closure_timing_and_runner_integration',
        scope=oldfit['scope']+' A separately qualified checkpoint-aware source adapter consumes the full parallel numerical closure. '
            'Complete parallel timing is software qualified; its actual full run remains gated. Fitting export/reader runner '
            'integration and installed calibrated resources remain required. No fitting worker is queued.')
    assert not Path(fit['output']).exists();write(fp,fit);bind(pins,fp)
    oldr=Path('metadata/full_weighted_shared_entity_timing_resources_20261004_v2.json')
    resource=json.loads(oldr.read_text());bind(pins,oldr)
    control=Path('metadata/retained_factor_watchpoint_debugger_execution_20261004_v1.json')
    observed=json.loads(control.read_text());bind(pins,control)
    assert observed['exit_code']==1 and observed['wall_seconds']>0
    # V8 has forty completed real probes but a retained post-exit debugger failure.
    # Its timing is a contextual observation, not a repaired full weighted benchmark.
    resource.update(status='prospective_complete_parallel_weighted_timing_resources',
        prepared_utc=datetime.now(timezone.utc).isoformat(),cpus=16,memory_gib=200,workers=16,
        worker_address_space_gib=12,reservation_capacity_gib=192,parent_headroom_gib=8,
        parent_soft_address_space_gib=8,parent_hard_address_space_gib=24,
        native_cpu_seconds_per_process=4838400,native_cpu_cap_is_eta=False,
        runtime_planning_core_hours_per_stage=[1000,12000],
        planning_wall_hours_per_stage_at_16_cores=[62.5,750],runtime_calibrated=False,
        original_v8_40_probe_wall_seconds=observed['wall_seconds'],
        full_parallel_timing_scaling_measured=False,
        available_ram_gib=psutil.virtual_memory().available/2**30,free_disk_gib=psutil.disk_usage('.').free/2**30,
        scope='Full20.832millioncandidate census, maximum694400eligible groups, two actual variance '
            'points and full native replay, with16bounded fork workers/12GiB AS each,8GiB parent '
            'soft AS and200GiB/no swap/BLAS1. Review-array custody is serialized across workers. '
            'CPU allocation of eight weeks per process is a cap, not a finish ETA or measured '
            'runtime. Planning range1000–12000core hours per stage uses original V8 context '
            'as a rough order-of-magnitude reference, includes no optimizer fits and is not '
            'calibrated for four-control full data. V8 debugger exit1 remains retained. '
            'Original768GiB scratch/868GiB minimum-free-disk arithmetic retained; no compressed '
            'or cross-cohort sharing assumed. Full numerical closure gates native admission. '
            'No new charge, GPU, old native restart or biological acceptance.')
    assert resource['maximum_original_candidates']==20832000 and resource['maximum_eligible_groups']==694400
    assert resource['output_scratch_allowance_gib']==768 and resource['minimum_free_disk_gib']==868
    assert psutil.virtual_memory().available>=200*2**30 and psutil.disk_usage('.').free>=868*2**30
    rp=Path('metadata/full_weighted_timing_parallel_resources_20261004_v1.json');write(rp,resource);bind(pins,rp)
    dependency='metadata/full_weighted_parallel_closure_launch_20261004_v1.json';bind(pins,Path(dependency))
    pp=Path('metadata/full_weighted_timing_parallel_plan_20261004_v1.json')
    plan=dict(status='qualified_complete_parallel_timing_pending_original_full_parallel_numerical_closure',
        fit_plan=str(fp),fit_plan_sha256=sha(fp),scaled_variance_points=[0.,1.],
        numerical_closure_launch=dependency,resources=resource,pins=pins,
        output='results/phylogeny/full-four-control-shared-entity-parallel-timing-20261004-v1',
        completion='metadata/full_weighted_timing_parallel_completed_20261004_v1.json',
        prepared_utc=datetime.now(timezone.utc).isoformat(),project_modules=sorted(map(str,modules)),
        source_adapter='full_weighted_shared_entity_fit_sources_parallel_v1',production_fitting_launched=False,
        scope='Every original4340cohorts/130200designs/260400response inputs and20832000candidate '
            'identities, allfourpolicies/twoloadingmodes/fivetrees/twooutcomes/MLandREML. Complete '
            'deterministic representative census, actual-D component/latent qualification, '
            'component/spectral likelihoods atvariance0and1, every numeric group independently '
            'replayed. All original exclusions and construction/precision reviews retained with '
            'content-addressed original arrays. Cached source and submitted array bytes checked '
            'around every cohort; bounded task queue stops new admissions on failure. Both '
            'original journals and all producer/reader checkpoints gate full closure. Hardware '
            'timings are contextual observations and conditional fit-budget planning, not an '
            'optimizer runtime bound. Exact original full parallel numerical closer gates '
            'native execution. No original serial/timing plan changed, no native retry, pilot, '
            'GPU or paid infrastructure. Fitting/calibration, adequate ancestry and all8aims '
            'remain required.')
    assert not Path(plan['output']).exists();verify(pins);write(pp,plan)
    command=['/usr/bin/prlimit','--as='+str(8*2**30)+':'+str(24*2**30),
        '--cpu=4838400','--fsize='+str(2*2**30),'--',sys.executable,
        'scripts/full_weighted_timing_parallel_v1.py','--plan',str(pp)]
    wp=Path('metadata/full_weighted_timing_parallel_producer_wait_plan_20261004_v1.json')
    write(wp,dict(dependencies=[dependency],command=command,
        pins={str(path):sha(path) for path in [pp,Path(dependency),Path('scripts/run_after_verified_dependencies_v2.py')]}))
    print(json.dumps(dict(status=plan['status'],plan=str(pp),plan_sha256=sha(pp),
        wait_plan=str(wp),candidate_rows=20832000,maximum_groups=694400,
        resources=resource),indent=2))


if __name__=='__main__':main()
