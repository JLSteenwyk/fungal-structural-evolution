#!/usr/bin/env python3
"""Qualify unchanged fitting scope and actual CPU/file/address-space boundaries.

No model fit or production closure is synthesized. Header mutations are
synthetic contracts; four real children are explicitly non-fit software probes.
Previously closed full-grid numerical/export/admission proofs stay immutable.
"""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import shutil
from unittest.mock import patch

from ancestral_chain_attempt import sha
from full_weighted_parallel_fit_admission_v1 import capacity, FULL
from full_weighted_fit_exports import atomic, PRODUCER, READER, SUMMARY
from reference_measurement_union_sources import bind, verify
from weighted_parallel_fit_execution_runtime_v2 import (validate_models, role_command, native_limits,
    execute_native, validate_native_receipt)
from close_full_weighted_parallel_fit_execution_v2 import native_custody


def refused(call):
    try:call()
    except (AssertionError,ValueError,KeyError,TypeError,FileExistsError,FileNotFoundError):return
    raise AssertionError('Invalid contract unexpectedly accepted')


def models(fit,fp):
    admission=dict(logical_cases=75188,unique_cohorts=4340,candidate_rows=20832000,model_setting_rows=622080,
        fresh_admission_candidate_rows=20832000,fit_contract='synthetic-full-header-only',source_contract='synthetic-full-header-only',
        original_fit_plan_sha256=sha(fp),scientific_eligibility=False,admission_repeats_numeric_probes=False,
        unmeasured_review_candidate_coverage=10,conditional_budget_weighted_seconds=5.,
        planning=dict(group_planning_costs=[dict(conditional_budget_weighted_producer_seconds=2.,conditional_budget_weighted_reader_seconds=3.)]))
    operation=dict(expected=deepcopy(FULL),fit_contract=admission['fit_contract'],fit_plan_sha256=sha(fp),fit_plan=str(fp),
        producer_status=PRODUCER,reader_status=READER,summary_fields=deepcopy(SUMMARY),
        cached_and_native_input_memory_guards_required=True,
        producer_command=role_command(fp,fit,'producer'),reader_command=role_command(fp,fit,'reader'))
    resources=capacity(admission,fit)
    validate_models(operation,fit,admission,resources)
    return operation,admission,resources


def scope_mutations(fit,operation,admission,resources):
    rejected=[]
    for name in ['partial_expected','partial_admission','partial_fresh_admission','foreign_contract','foreign_source_contract',
        'guard_required_false','legacy_producer_command','legacy_reader_command',
        'foreign_original_plan_hash','drop_policy','drop_mode','drop_tree','drop_method','primary_gradient','primary_evaluations',
        'primary_iterations','variance_cap','column_batch','independent_gradient','replay_objective','curvature_step',
        'promoted_science','repeat_numeric','producer_command','reader_command','producer_status','reader_status','missing_summary',
        *['capacity_'+k for k in ['cpus','memory_gib','swap_gib','blas_threads','address_space_gib','per_file_limit_gib',
            'minimum_available_ram_gib','minimum_free_disk_gib','output_scratch_reserve_gib','producer_native_cpu_allocation_seconds',
            'reader_native_cpu_allocation_seconds','resources_installed','gpu','new_cost_usd','complete_runtime_calibration']]]:
        f,o,a,r=map(deepcopy,[fit,operation,admission,resources])
        if name=='partial_expected':o['expected']['candidate_rows']-=1
        elif name=='partial_admission':a['candidate_rows']-=1
        elif name=='partial_fresh_admission':a['fresh_admission_candidate_rows']-=1
        elif name=='foreign_contract':o['fit_contract']='foreign'
        elif name=='foreign_source_contract':a['source_contract']='foreign'
        elif name=='foreign_original_plan_hash':a['original_fit_plan_sha256']='foreign'
        elif name=='guard_required_false':o['cached_and_native_input_memory_guards_required']=False
        elif name=='legacy_producer_command':o['producer_command'][1]='scripts/prepare_full_weighted_shared_entity_fits.py'
        elif name=='legacy_reader_command':o['reader_command'][1]='scripts/readback_full_weighted_shared_entity_fits.py'
        elif name.startswith('drop_'):f[{'drop_policy':'policies','drop_mode':'loading_modes','drop_tree':'trees','drop_method':'methods'}[name]].pop()
        elif name in ['primary_gradient','primary_evaluations','primary_iterations','variance_cap']:
            k={'primary_gradient':'gradient_tolerance','primary_evaluations':'max_evaluations','primary_iterations':'max_iterations','variance_cap':'maximum_scaled_variance'}[name];f['optimizer'][k]*=2
        elif name=='column_batch':f['independent_audit']['column_batch']=3
        elif name=='independent_gradient':f['independent_audit']['optimizer']['gradient_tolerance']=3e-6
        elif name=='replay_objective':f['independent_audit']['replay']['objective_atol']=1e-5
        elif name=='curvature_step':f['independent_audit']['curvature']['coordinate_step']=.01
        elif name=='promoted_science':a['scientific_eligibility']=True
        elif name=='repeat_numeric':a['admission_repeats_numeric_probes']=True
        elif name in ['producer_command','reader_command']:o[name][-1]='foreign'
        elif name in ['producer_status','reader_status']:o[name]='foreign'
        elif name=='missing_summary':o['summary_fields'].pop()
        else:
            k=name.removeprefix('capacity_');r[k]=not r[k] if type(r[k]) is bool else r[k]+1
        refused(lambda:validate_models(o,f,a,r));rejected.append(name)
    for role in ['producer','reader']:
        for bad in [0,-1,2**63,1.5,True]:
            r=deepcopy(resources);r[role+'_native_cpu_allocation_seconds']=bad
            refused(lambda:native_limits(r,role));rejected.append([role,'invalid_cpu',bad])
    return rejected


def boundary_probes(root,bindings):
    # Private boundary overrides apply only to these explicitly non-fit
    # children. Production CLI admission still requires exact measured limits.
    resources=dict(cpus=2,memory_gib=32,swap_gib=0,blas_threads=1,
        minimum_available_ram_gib=16,minimum_free_disk_gib=128)
    limits=dict(cpu_seconds=3,address_space_bytes=128*2**20,per_file_bytes=64*2**10)
    rows=[];restart_checks=[]
    for name,body in [('observe','pass'),('cpu','while True: pass'),
        ('file',"with open('"+str(root/'file-limit-probe.bin')+"','wb') as f: f.write(b'x'*2**20)"),
        ('address_space','bytearray(256*2**20)')]:
        program=root/('software-boundary-'+name+'.py')
        program.write_text('import resource, time, json, os\ntime.sleep(.2)\n'+
            "print(json.dumps(dict(cpu=resource.getrlimit(resource.RLIMIT_CPU),address_space=resource.getrlimit(resource.RLIMIT_AS),file=resource.getrlimit(resource.RLIMIT_FSIZE),blas={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']})),flush=True)\n"+body+'\n')
        local=dict(bindings);bind(local,program)
        directory=root/name;command=[sys.executable,str(program)]
        record=execute_native(command,limits,resources,directory,local,'software_boundary_'+name)
        assert record['child']['created'] is not None
        native=json.loads((directory/'stdout.log').read_text().splitlines()[0])
        assert native['cpu']==[3,3] and native['address_space']==[128*2**20]*2 and native['file']==[64*2**10]*2
        assert set(native['blas'].values())=={'1'}
        if name=='observe':assert record['exit_code']==0 and record['status']=='original_native_weighted_fit_role_exited_zero'
        else:
            assert record['exit_code']!=0 and record['status']=='original_native_weighted_fit_role_failed_preserved'
            if name=='cpu':assert record['exit_code'] in [-9,-24]
            if name=='file':assert 'File too large' in (directory/'stderr.log').read_text() or record['exit_code']==-25
            if name=='address_space':assert 'MemoryError' in (directory/'stderr.log').read_text()
        before={str(p):sha(p) for p in directory.iterdir()}
        refused(lambda:execute_native(command,limits,resources,directory,local,'software_boundary_'+name))
        assert before=={str(p):sha(p) for p in directory.iterdir()}
        restart_checks.append(name);rows.append(dict(case=name,receipt=str(directory/'receipt.json'),exit_code=record['exit_code'],native_limits=limits,nonfit_software_probe=True))
    good=json.loads((root/'observe'/'receipt.json').read_text());config=json.loads((root/'observe'/'configuration.json').read_text());process=json.loads((root/'observe'/'process.json').read_text())
    validate_native_receipt(good,config,process,good['child']['command'],limits,resources)
    native_custody(root/'observe',good['role'],good['child']['command'],limits,resources,good['source_hashes'],{})
    bad=[]
    for name in ['nonzero_exit','failure_status','promoted_science','restarted','config_science','config_role','cpu_limit','command',
        'child_identity','wrapper_identity','missing_invocation','different_invocation','blas','config_source','config_cgroup',
        'cgroup_path','cpu_quota','memory_quota','swap_quota','swap_used','oom_event','foreign_wrapper_binary']:
        r,c,p=map(deepcopy,[good,config,process])
        if name=='nonzero_exit':r['exit_code']=1
        elif name=='failure_status':r['status']='failure'
        elif name=='promoted_science':r['scientific_eligibility']=True
        elif name=='restarted':r['restarted']=True
        elif name=='config_science':c['scientific_eligibility']=True
        elif name=='config_role':c['role']='foreign'
        elif name=='cpu_limit':r['native_limits']['cpu_seconds']+=1
        elif name=='command':r['child']['command'][-1]='foreign'
        elif name=='child_identity':p['created']+=1
        elif name=='wrapper_identity':c['wrapper']['pid']+=1
        elif name=='missing_invocation':r['invocation_id']=None
        elif name=='different_invocation':c['invocation_id']='0'*32
        elif name=='blas':r['blas_environment']['OMP_NUM_THREADS']='2'
        elif name=='config_source':c['source_hashes']['foreign']='foreign'
        elif name=='config_cgroup':c['cgroup_before']['path']='foreign'
        elif name=='cgroup_path':r['cgroup_after']['path']='foreign'
        elif name in ['cpu_quota','memory_quota','swap_quota']:r['cgroup_after']['limits'][{'cpu_quota':'cpu.max','memory_quota':'memory.max','swap_quota':'memory.swap.max'}[name]]='foreign'
        elif name=='swap_used':r['cgroup_after']['memory_bytes']['memory.swap.current']=1
        elif name=='oom_event':r['cgroup_after']['memory_events']['oom']=1
        else:r['wrapper']['cmdline'][0]='foreign';c['wrapper']=deepcopy(r['wrapper'])
        refused(lambda:validate_native_receipt(r,c,p,good['child']['command'],limits,resources));bad.append(name)
    # Corrupt only independent copies of this turn's software attempt.
    # The actual original receipt/logs and all published parents stay intact.
    for name in ['omit_artifact','extra_artifact','changed_stdout','changed_config','foreign_source_hash','wrong_role']:
        directory=root/('copied-custody-'+name);shutil.copytree(root/'observe',directory)
        record=deepcopy(good);record['artifacts']={str(directory/Path(q).name):sha(directory/Path(q).name) for q in good['artifacts']}
        if name=='omit_artifact':record['artifacts'].pop(str(directory/'stderr.log'))
        elif name=='extra_artifact':record['artifacts']['foreign']='foreign'
        elif name=='changed_stdout':(directory/'stdout.log').write_text('changed\n')
        elif name=='changed_config':
            changed=json.loads((directory/'configuration.json').read_text());changed['native_limits']['cpu_seconds']+=1
            (directory/'configuration.json').write_text(json.dumps(changed)+'\n');record['artifacts'][str(directory/'configuration.json')]=sha(directory/'configuration.json')
        elif name=='foreign_source_hash':record['source_hashes']['foreign']='foreign'
        else:record['role']='foreign'
        (directory/'receipt.json').write_text(json.dumps(record)+'\n')
        refused(lambda:native_custody(directory,good['role'],good['child']['command'],limits,resources,good['source_hashes'],{}));bad.append(name)
    return rows,restart_checks,bad


def successful_preparation_control(root, request, fit, admission):
    """Exercise operational preparation with explicitly mocked full closures.

    Actual qualified admission/export proofs are consumed unchanged. No
    source acceptance or completed biological stage is fabricated: the
    resulting private operation must fail real runtime configuration because
    its actual numerical/timing closures are still absent.
    """
    from prepare_full_weighted_parallel_fit_execution_v2 import run as prepare_operations, headers
    from weighted_parallel_fit_execution_runtime_v2 import configuration
    private = deepcopy(request)
    for key in ['operational_output','resources_output','admission_output']:
        private[key] = str(root/('synthetic-'+key+'.json'))
    pp = root/'synthetic-successful-preparation-request.json'
    atomic(pp, private)
    actual_fit, actual_timing = headers(private)
    synthetic = deepcopy(admission)
    synthetic.update(original_fit_plan=private['fit_plan'], timing_plan=private['timing_plan'],
        timing_completion=private['timing_completion'], timing_completion_sha256='synthetic-absent-closure')
    with patch('prepare_full_weighted_parallel_fit_execution_v2.requirements',
               return_value=(actual_fit,actual_timing,[])), \
         patch('prepare_full_weighted_parallel_fit_execution_v2.closed_timing',
               return_value=(dict(fit_contract=synthetic['fit_contract']),{},synthetic)), \
         patch('prepare_full_weighted_parallel_fit_execution_v2.fingerprint',return_value=None), \
         patch('prepare_full_weighted_parallel_fit_execution_v2.journal_terminal',
               return_value=dict(synthetic_terminal_control=True)):
        operation = prepare_operations(pp)
    assert operation['cached_and_native_input_memory_guards_required'] is True
    assert operation['launch_state']=='not_launched_or_queued' and operation['fits_computed']==0
    assert operation['producer_command']==role_command(private['fit_plan'],fit,'producer')
    assert operation['reader_command']==role_command(private['fit_plan'],fit,'reader')
    assert not Path(fit['output']).exists()
    controller = root/'synthetic-prepared-controller.json'
    atomic(controller, dict(operation=private['operational_output'], expected=deepcopy(FULL),
        pins={private['operational_output']:sha(private['operational_output'])}))
    try:configuration(controller)
    except FileNotFoundError as error:
        assert str(Path(private['timing_completion'])) in str(error)
    else:raise AssertionError('Synthetic successful branch admitted an actual fit')
    atomic(root/'successful-preparation-control-scope.json',dict(
        successful_branch_exercised=True, actual_qualified_software_proofs_consumed=True,
        source_closures_and_terminal_checks_explicitly_mocked=True,
        actual_runtime_refused_missing_timing_closure=True, production_root_created=False,
        scientific_eligibility=False))
    return True


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    assert not a.receipt.exists();a.output.mkdir(parents=True,exist_ok=False)
    parent_paths=['metadata/weighted_parallel_fit_admission_software_validation_20261004_v1.json',
        'metadata/weighted_parallel_fit_admission_software_transport_20261004_v1.json',
        'metadata/weighted_parallel_source_fitting_software_validation_20261004_v3.json',
        'metadata/weighted_parallel_source_fitting_software_transport_20261004_v3.json']
    bindings={}
    for path in parent_paths:
        parent=json.loads(Path(path).read_text());verify(parent['source_hashes']);bind(bindings,path)
        for path,d in parent['source_hashes'].items():bind(bindings,path,d)
    for path in ['scripts/'+n+'.py' for n in ['weighted_parallel_fit_execution_runtime_v2','close_full_weighted_parallel_fit_execution_v2',
        'launch_full_weighted_parallel_fit_execution_v2','check_full_weighted_parallel_fit_controller_v2','prepare_full_weighted_parallel_fit_execution_v2','run_weighted_fit_controller_software_stage']]:bind(bindings,path)
    fp=Path('metadata/full_weighted_parallel_source_fit_draft_plan_20261004_v1.json');fit=json.loads(fp.read_text());verify(fit['pins']);bind(bindings,fp)
    operation,admission,resources=models(fit,fp)
    bad=scope_mutations(fit,operation,admission,resources)
    atomic(a.output/'synthetic_full_header_contract.json',dict(operation=operation,admission=admission,resources=resources,synthetic_only=True))
    probes,restarts,custody=boundary_probes(a.output,bindings)
    from launch_full_weighted_parallel_fit_execution_v2 import prepare, REQUEST, CONTROLLER
    request=json.loads(Path(REQUEST).read_text());verify(request['pins']);bind(bindings,REQUEST)
    absent=[CONTROLLER,fit['output'],request['operational_output'],request['resources_output'],request['admission_output']]
    assert all(not Path(q).exists() for q in absent)
    # Actual pending production path: any submit attempt is a test failure.
    with patch('launch_full_weighted_parallel_fit_execution_v2.launch',side_effect=AssertionError('No fit launch authorized before closure')) as submit:
        pending=prepare();assert not submit.called
    assert pending['status']=='pending_original_numerical_and_timing_closure_no_fit_prepared'
    assert pending['missing_prerequisites'] and not pending['fit_resources_installed'] and not pending['production_fitting_launched_or_queued']
    assert all(not Path(q).exists() for q in absent)
    atomic(a.output/'actual_pending_production_refusal.json',pending)
    successful_preparation = successful_preparation_control(a.output,request,fit,admission)
    for q in a.output.rglob('*'):
        if q.is_file():bind(bindings,q)
    verify(bindings)
    result=dict(status='passed_full_weighted_parallel_fit_controller_scope_and_actual_native_boundary_contracts_v2',
        checked_utc=datetime.now(timezone.utc).isoformat(),expected_production=FULL,original_fit_plan_sha256=sha(fp),
        actual_native_boundary_probes=4,actual_successful_boundary_probes=1,actual_preserved_boundary_failures=3,
        actual_boundary_cases=probes,rejected_scope_and_capacity_changes=len(bad),scope_and_capacity_rejections=bad,
        rejected_native_custody_changes=len(custody),native_custody_rejections=custody,
        completed_and_failed_attempt_restarts_refused=True,restart_cases=restarts,actual_missing_prerequisite_refusal=True,
        successful_operational_preparation_branch_qualified=successful_preparation,
        successful_branch_closures_explicitly_mocked=True,
        production_fitting_launched_or_queued=False,fits_computed=0,scientific_eligibility=False,gpu=False,new_cost_usd=0,
        source_hashes=bindings,
        scope='Unchanged full20.832million candidate/49.7664million link model header and actual pending production refusal; synthetic header/capacity/custody changes, four actual non-fit CPU/file/address-space/resource-observation children, original exit/identity/log preservation and restart refusal. Prior complete-grid fitting math/export/admission proofs inherited and rehashed, not repeated; no real fit or production closure invented. Actual production launch and fit/readback closure remain unexecuted behind both original gates; all eight aims and inferential calibration open.')
    atomic(a.receipt,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','scope_and_capacity_rejections','native_custody_rejections']},indent=2))


if __name__=='__main__':main()
