"""Exact-role fitting commands, operational admission and native resource custody."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import psutil

from ancestral_chain_attempt import sha
from full_weighted_fit_admission import capacity, production_contract, FULL
from full_weighted_fit_exports import atomic, PRODUCER, READER, SUMMARY
from reference_measurement_union_sources import bind, verify


def validate_models(operation,fit,admission,resources):
    production_contract(fit,admission)
    assert operation['expected']==FULL and operation['fit_contract']==admission['fit_contract']==admission['source_contract']
    assert operation['fit_plan_sha256']==admission['original_fit_plan_sha256']
    assert admission['fresh_admission_candidate_rows']==FULL['candidate_rows']
    assert admission['scientific_eligibility'] is False and admission['admission_repeats_numeric_probes'] is False
    assert resources==capacity(admission,fit)
    assert resources['resources_installed'] is False and resources['gpu'] is False and resources['new_cost_usd']==0
    assert (resources['cpus'],resources['memory_gib'],resources['swap_gib'],resources['blas_threads'],
        resources['address_space_gib'],resources['per_file_limit_gib'])==(2,32,0,1,24,2)
    assert resources['minimum_free_disk_gib']==3172 and resources['output_scratch_reserve_gib']==3072
    assert operation['producer_status']==PRODUCER and operation['reader_status']==READER and operation['summary_fields']==SUMMARY
    for role in ['producer','reader']:assert operation[role+'_command']==role_command(operation['fit_plan'],fit,role)


def role_command(fit_path,fit,role):
    if role=='producer':return [sys.executable,'scripts/prepare_full_weighted_shared_entity_fits.py','--plan',str(fit_path)]
    assert role=='reader'
    return [sys.executable,'scripts/readback_full_weighted_shared_entity_fits.py','--plan',str(fit_path),
        '--output',str(Path(fit['output'])/'readback.json')]


def cgroup_state():
    group=next(l[3:] for l in Path('/proc/self/cgroup').read_text().splitlines() if l.startswith('0::'))
    root=Path('/sys/fs/cgroup')/group.lstrip('/')
    return dict(path=str(root),limits={k:(root/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']},
        memory_bytes={k:int((root/k).read_text()) for k in ['memory.current','memory.peak','memory.swap.current']},
        memory_events={k:int(v) for k,v in (l.split() for l in (root/'memory.events').read_text().splitlines())})


def validate_limits(state,resources,environment):
    assert state['limits']=={'cpu.max':str(resources['cpus']*100000)+' 100000',
        'memory.max':str(resources['memory_gib']*2**30),'memory.swap.max':'0'}
    assert resources['swap_gib']==0 and resources['blas_threads']==1
    assert all(environment.get(k)=='1' for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'])
    assert state['memory_bytes']['memory.swap.current']==0
    assert all(state['memory_events'][k]==0 for k in ['oom','oom_kill','oom_group_kill'])


def configuration(plan_path):
    plan_path=Path(plan_path);plan=json.loads(plan_path.read_text());verify(plan['pins'])
    op_path=Path(plan['operation']);operation=json.loads(op_path.read_text());verify(operation['pins'])
    assert sha(op_path)==plan['pins'][str(op_path)]
    assert operation['status']=='prepared_full_weighted_fit_operations_not_launched'
    assert operation['launch_state']=='not_launched_or_queued' and operation['fits_computed']==0
    assert operation['expected']==plan['expected']==FULL
    assert all(operation[k] is False for k in ['scientific_eligibility','nonuniform_weighting_accepted','component_variance_attribution_accepted','gpu'])
    fp=Path(operation['fit_plan']);fit=json.loads(fp.read_text());verify(fit['pins'])
    assert sha(fp)==operation['fit_plan_sha256'] and fit['launch_state']=='not_launched_or_queued'
    ap=Path(operation['admission']);admission=json.loads(ap.read_text());verify(admission['source_hashes'])
    assert admission['original_fit_plan']==str(fp) and admission['original_fit_plan_sha256']==sha(fp)
    resources=json.loads(Path(operation['resources']).read_text());validate_models(operation,fit,admission,resources)
    timing_path=Path(admission['timing_plan']);timing=json.loads(timing_path.read_text())
    completion_path=Path(admission['timing_completion']);completion=json.loads(completion_path.read_text())
    assert timing['fit_plan']==str(fp) and timing['fit_plan_sha256']==sha(fp)
    assert timing['completion']==str(completion_path) and sha(completion_path)==admission['timing_completion_sha256']
    assert completion['status']=='complete_verified_full_four_control_input_timing_accounting_v2'
    assert completion['logical_cases']==FULL['logical_cases'] and completion['candidate_rows']==FULL['candidate_rows']
    assert completion['unique_cohorts']==FULL['cohorts'] and completion['exact_process_journals_checked']==2
    assert sha(completion['full_hash_archive'])==completion['full_hash_archive_sha256']
    numeric=json.loads(Path(fit['qualification_completion']).read_text())
    assert numeric['status']=='complete_verified_full_four_control_covariance_numerical_qualification_v1'
    assert numeric['cohorts']==FULL['cohorts'] and numeric['designs']==FULL['designs'] and numeric['exact_process_journals_checked']==2
    assert sha(numeric['full_hash_archive'])==numeric['full_hash_archive_sha256']
    bindings=dict(plan['pins'])
    for mapping in [operation['pins'],fit['pins'],admission['source_hashes']]:
        for p,d in mapping.items():bind(bindings,p,d)
    for p in [plan_path,op_path,fp,ap,operation['resources']]:bind(bindings,p)
    verify(bindings)
    return plan,operation,fit,resources,bindings


def native_limits(resources,role):
    assert role in ['producer','reader']
    cpu=resources[role+'_native_cpu_allocation_seconds']
    assert type(cpu) is int and 0<cpu<2**63
    return dict(cpu_seconds=cpu,address_space_bytes=resources['address_space_gib']*2**30,
        per_file_bytes=resources['per_file_limit_gib']*2**30)


def validate_native_receipt(record,config,process,command,limits,resources):
    """Check original custody semantics separately from file hash checks."""
    assert record['status']=='original_native_weighted_fit_role_exited_zero' and record['exit_code']==0
    assert record['scientific_eligibility'] is record['restarted'] is False
    assert config['scientific_eligibility'] is False and config['role']==record['role']
    assert record['native_limits']==config['native_limits']==limits
    assert record['child']['command']==config['command']==command
    assert process==record['child'] and config['wrapper']==record['wrapper']
    wrapper=record['wrapper'];child=record['child']
    assert type(wrapper['pid']) is int and type(child['pid']) is int and 0<child['pid']!=wrapper['pid']>0
    assert wrapper['created']>0 and child['created']>=wrapper['created']
    assert isinstance(wrapper['cmdline'],list) and wrapper['cmdline'][0]==sys.executable
    import re
    assert re.fullmatch('[0-9a-f]{32}',record['invocation_id'])
    assert config['invocation_id']==record['invocation_id']
    assert config['blas_environment']==record['blas_environment']
    assert config['source_hashes']==record['source_hashes']
    assert config['cgroup_before']==record['cgroup_before']
    assert record['cgroup_before']['path']==record['cgroup_after']['path']
    for state in [record['cgroup_before'],record['cgroup_after']]:validate_limits(state,resources,record['blas_environment'])


def execute_native(command,limits,resources,directory,bindings,role):
    """One immutable attempt; caller supplies a validated exact-role command.

    Tests may exercise this boundary with explicitly labeled non-fit children.
    The production CLI validates configuration/commands before reaching it.
    """
    assert all(type(limits[k]) is int and 0<limits[k]<2**63 for k in ['cpu_seconds','address_space_bytes','per_file_bytes'])
    assert isinstance(command,list) and command[0]==sys.executable and all(isinstance(p,str) for p in command)
    verify(bindings);before=cgroup_state();validate_limits(before,resources,os.environ)
    assert psutil.virtual_memory().available>=resources['minimum_available_ram_gib']*2**30
    assert shutil.disk_usage('.').free>=resources['minimum_free_disk_gib']*2**30
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=False)
    wrapper=psutil.Process();identity=dict(pid=wrapper.pid,created=wrapper.create_time(),cmdline=wrapper.cmdline())
    bounded=['/usr/bin/prlimit','--as='+str(limits['address_space_bytes']),'--cpu='+str(limits['cpu_seconds']),
        '--fsize='+str(limits['per_file_bytes']),'--',*command]
    input_pins=dict(bindings);bind(input_pins,sys.executable);bind(input_pins,'/usr/bin/prlimit')
    started=datetime.now(timezone.utc).isoformat();tick=time.monotonic();usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    blas={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']}
    atomic(directory/'configuration.json',dict(role=role,command=bounded,wrapper=identity,
        invocation_id=os.environ.get('INVOCATION_ID'),started_utc=started,native_limits=limits,
        cgroup_before=before,blas_environment=blas,source_hashes=input_pins,scientific_eligibility=False))
    print(json.dumps(dict(original_fit_execution_wrapper=identity,role=role,invocation_id=os.environ.get('INVOCATION_ID'))),flush=True)
    with (directory/'stdout.log').open('x') as out,(directory/'stderr.log').open('x') as err:
        child=subprocess.Popen(bounded,stdout=out,stderr=err,start_new_session=True)
        try:created=psutil.Process(child.pid).create_time()
        except psutil.NoSuchProcess:created=None
        atomic(directory/'process.json',dict(pid=child.pid,created=created,command=bounded))
        code=child.wait()
    after=cgroup_state();validate_limits(after,resources,os.environ);verify(input_pins)
    last=resource.getrusage(resource.RUSAGE_CHILDREN)
    result=dict(status=('original_native_identity_missing_requires_review' if created is None else
        'original_native_weighted_fit_role_exited_zero' if code==0 else 'original_native_weighted_fit_role_failed_preserved'),
        role=role,wrapper=identity,child=dict(pid=child.pid,created=created,command=bounded),
        invocation_id=os.environ.get('INVOCATION_ID'),exit_code=code,native_limits=limits,
        checked_utc=datetime.now(timezone.utc).isoformat(),wall_seconds=time.monotonic()-tick,
        native_cpu_seconds=last.ru_utime-usage.ru_utime+last.ru_stime-usage.ru_stime,
        children_cumulative_peak_rss_bytes=last.ru_maxrss*1024,cgroup_before=before,cgroup_after=after,blas_environment=blas,
        source_hashes=input_pins,artifacts={str(p):sha(p) for p in directory.iterdir() if p.is_file()},
        scientific_eligibility=False,restarted=False,
        scope='Exact wrapper/native identity, supplied role command and actual cgroup/BLAS/native limits; original exit and log/input bytes preserved. RSS is the cumulative child high-water statistic for this wrapper, not an isolated per-child peak. CPU allocation is capacity, not ETA or convergence guarantee. No automatic retry, successful-output rewrite or biological acceptance. Software boundary probes must not be represented as native model fits.')
    atomic(directory/'receipt.json',result)
    print(json.dumps(dict(role=role,exit_code=code,wall_seconds=result['wall_seconds'],native_cpu_seconds=result['native_cpu_seconds'])),flush=True)
    return result


def run(plan_path,role):
    plan,operation,fit,resources,bindings=configuration(plan_path)
    root=Path(fit['output'])
    if role=='producer':assert not (root/'receipt.json').exists(), 'Completed fit producer cannot restart'
    else:
        assert role=='reader' and not (root/'independent_readback_completed.json').exists()
        producer=json.loads((root/'receipt.json').read_text())
        assert producer['status']==PRODUCER and producer['source_contract']==operation['fit_contract']
        bind(bindings,root/'receipt.json')
    result=execute_native(role_command(operation['fit_plan'],fit,role),native_limits(resources,role),resources,
        Path(plan['runtime_root'])/role,bindings,role)
    if result['exit_code']!=0:raise SystemExit(result['exit_code'] if result['exit_code']>0 else 1)
    assert result['status']=='original_native_weighted_fit_role_exited_zero'
    expected=root/('receipt.json' if role=='producer' else 'readback.json')
    record=json.loads(expected.read_text())
    assert record['status']==(PRODUCER if role=='producer' else READER)
    assert record['source_contract']==operation['fit_contract'] and record['plan_sha256']==operation['fit_plan_sha256']
    assert record['candidate_rows']==FULL['candidate_rows'] and record['setting_fit_links']==FULL['setting_fit_links']
    assert record['scientific_eligibility'] is record['nonuniform_weighting_accepted'] is record['component_variance_attribution_accepted'] is False
    atomic(Path(plan['runtime_root'])/role/'role_output.json',dict(receipt=str(expected),sha256=sha(expected),status=record['status'],
        fit_contract=operation['fit_contract'],scientific_eligibility=False))


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--role',choices=['producer','reader'],required=True)
    a=p.parse_args();run(a.plan,a.role)
