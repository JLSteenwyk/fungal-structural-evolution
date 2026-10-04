#!/usr/bin/env python3
"""Exercise actual timing workers and source/task byte guards in private fork memory."""
import argparse
from datetime import datetime,timezone
import json
import multiprocessing
from pathlib import Path

import numpy as np
import psutil

import full_weighted_timing_parallel_v1 as target
from full_weighted_shared_entity_fit_sources_parallel_v1 import cohorts
from full_weighted_shared_entity_timing import groups
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind,verify


def write(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2);f.write('\n')


def worker(source,fit,plan,root,task,case,custody):
    root=Path(root);root.mkdir();(root/'cohorts').mkdir();(root/'checkpoints').mkdir()
    target.worker_init(source,fit,plan,str(root),False,custody)
    factor=source['factors'][next(iter(source['factors']))]
    matrix=task[3][0][1]
    before_factor=target.array_bytes(factor);before_matrix=target.array_bytes(matrix)
    def flip(value):
        scalar=np.asarray([value.flat[0]],dtype=value.dtype)
        byte=scalar.view(np.uint8);byte.flat[0]=int(byte.flat[0])^1
        value.flat[0]=scalar[0]
    phase={'source_before':'before','source_after':'after','task_after':'submitted-task-after'}.get(case)
    original=target.source_guard
    def injected_guard(root,cohort,refs,before,observed_phase):
        if observed_phase==phase:flip(matrix if case=='task_after' else factor)
        return original(root,cohort,refs,before,observed_phase)
    target.source_guard=injected_guard
    try:value=target.cohort_job(task)
    except AssertionError as error:
        assert case!='control' and str(error)=='Cached timing source changed '+phase,(case,str(error))
        outcome=dict(rejected=True,error=str(error))
    else:
        assert case=='control' and value['timing_groups']>0
        assert value['cached_numeric_inputs_preserved'] and value['submitted_numeric_inputs_preserved']
        outcome=dict(rejected=False,timing_groups=value['timing_groups'])
    factor_after=target.array_bytes(factor);matrix_after=target.array_bytes(matrix)
    assert (factor_after==before_factor)==(case in ['control','task_after'])
    assert (matrix_after==before_matrix)==(case!='task_after')
    if case!='control':
        failures=list((root/'failures').glob('*/failure.json'));assert len(failures)==1
        assert list(failures[0].parent.glob('*.npy'))
    process=psutil.Process()
    write(root/'guard.json',dict(case=case,**outcome,worker=dict(pid=process.pid,
        created=process.create_time(),cmdline=process.cmdline()),factor_before=before_factor,
        factor_after=factor_after,matrix_before=before_matrix,matrix_after=matrix_after,
        deliberate_injection_into_private_fork_memory=True,native_arithmetic_mocked=False,
        scientific_eligibility=False))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False);assert not a.receipt.exists()
    pins={};project_sources(pins,[Path(__file__)])
    fp=Path('data/software_audits/full-weighted-parallel-fit-source-adapter-20261004-v1/qualified-common-pair-parallel.fit.plan.json')
    plan=dict(fit_plan=str(fp),fit_plan_sha256=target.sha(fp),output=str(a.output/'unused-stage'),
        scaled_variance_points=[0.,1.],pins=dict(pins),resources=dict(cpus=2,memory_gib=32,
            workers=2,worker_address_space_gib=12,reservation_capacity_gib=24,
            parent_headroom_gib=8,minimum_free_disk_gib=128),scope='Private fork-memory guard controls; no producer-stage closure or biological fit.')
    pp=a.output/'private.plan.json';write(pp,plan);target.runtime(plan)
    fit,source,bindings=target.sources(plan,pp)
    for path,digest in bindings.items():bind(pins,path,digest)
    task=None
    for index,(cohort,rows,entries) in enumerate(cohorts(source,fit)):
        _,selected,_=groups(source,fit,cohort,entries)
        if selected:task=(index,cohort,rows,entries);break
    assert task is not None
    refs=target.distinct_refs(list(target.arrays(source))+list(target.arrays(task,'task')))
    before={name:target.array_bytes(value) for name,value in refs}
    context=multiprocessing.get_context('fork');custody=context.Lock();cases=[]
    for case in ['control','source_before','source_after','task_after']:
        root=a.output/case
        process=context.Process(target=worker,args=(source,fit,dict(plan,_original_path=str(pp)),str(root),task,case,custody))
        process.start();native=psutil.Process(process.pid);created=native.create_time()
        process.join(180);assert not process.is_alive(),'Private guard child exceeded observation budget'
        assert process.exitcode==0,(case,process.exitcode)
        report=json.loads((root/'guard.json').read_text())
        assert report['worker']['pid']==process.pid and report['worker']['created']==created
        cases.append(report)
        for path in root.rglob('*'):
            if path.is_file():bind(pins,path)
    assert {name:target.array_bytes(value) for name,value in refs}==before
    verify(pins)
    result=dict(status='passed_private_timing_source_and_submitted_task_memory_guards_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),cases=cases,
        parent_source_and_submitted_task_arrays_preserved=True,source_hashes=pins,
        scientific_eligibility=False,
        scope='One complete eligible cohort with actual native timing groups, plus deliberate '
              'private source changes before admission/after arithmetic and a submitted X '
              'change after arithmetic. Actual guard callbacks reject all three and retain '
              'mutated arrays. Injection hooks do not mock native arithmetic. Parent arrays '
              'and original files remain unchanged. No repair of original corruption, full '
              'production timing, fitting or biological acceptance.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','cases']},indent=2))


if __name__=='__main__':main()
