#!/usr/bin/env python3
"""Complete weighted timing census and numeric replay with bounded fork workers."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from datetime import datetime,timezone
import fcntl
import gzip
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import resource
import shutil
import time

import numpy as np
import psutil

from full_weighted_covariance_qualification import runtime
from full_weighted_covariance_qualification_parallel_v1 import arrays,array_bytes
from full_weighted_shared_entity_fit_sources_parallel_v1 import cohorts,operators_for
from full_weighted_shared_entity_timing import groups,probe,estimate
from full_weighted_timing_contracts_v2 import validate_probes,replay,review_inputs
from full_weighted_fit_exports import atomic,finish_gzip,failure_capture
from prepare_weighted_timing_parallel_source_reference_v1 import sources,SCHEMA,PRODUCER,READER
from reference_measurement_union_sources import bind,verify
from ancestral_chain_attempt import sha

CACHE=None
SUMMARY=['logical_cases','model_setting_rows','unique_cohorts','candidate_rows',
    'eligible_candidates','timing_groups','source_status_counts','timing_status_counts',
    'measured_candidate_coverage','unmeasured_review_candidate_coverage',
    'conditional_budget_weighted_seconds']


def distinct_refs(refs):
    seen=set();result=[]
    for name,value in refs:
        if id(value) not in seen:
            seen.add(id(value));result.append((name,value))
    return result


def worker_init(source,fit,plan,root,reader,custody_lock):
    global CACHE
    limit=plan['resources']['worker_address_space_gib']*2**30
    resource.setrlimit(resource.RLIMIT_AS,(limit,limit))
    refs=distinct_refs(arrays(source));assert refs
    before={name:array_bytes(value) for name,value in refs}
    CACHE=source,fit,plan,Path(root),reader,custody_lock,refs,before


def source_guard(root,cohort,refs,before,phase):
    after={name:array_bytes(value) for name,value in refs}
    if after!=before:
        directory=root/'failures'/(cohort['cohort_id']+'-cached-source-'+phase)
        directory.mkdir(parents=True,exist_ok=False)
        atomic(directory/'failure.json',dict(cohort_id=cohort['cohort_id'],phase=phase,
            before=before,after=after,scientific_eligibility=False))
        for name,value in refs:
            if before[name]!=after[name]:
                with (directory/(hashlib.sha256(name.encode()).hexdigest()+'.npy')).open('xb') as f:np.save(f,value,allow_pickle=False)
        raise AssertionError('Cached timing source changed '+phase)


def cohort_job(task):
    source,fit,plan,root,reader,custody_lock,refs,before=CACHE
    index,cohort,rows,entries=task
    source_guard(root,cohort,refs,before,'before')
    task_refs=distinct_refs(arrays(task,'submitted-task'))
    task_before={name:array_bytes(value) for name,value in task_refs}
    assert shutil.disk_usage(root).free>=plan['resources']['minimum_free_disk_gib']*2**30
    census,selected,counts=groups(source,fit,cohort,entries);cid=cohort['cohort_id']
    folder=root/'cohorts';cp=folder/(cid+'.receipt.json')
    fp=folder/(cid+'.census.jsonl.gz');pp=folder/(cid+'.probes.json')
    stage=dict(schema=SCHEMA,plan_sha256=sha(plan['_original_path']),fit_contract=source['fit_contract'])
    if reader:
        saved=json.loads(cp.read_text());probes=json.loads(pp.read_text())
        with gzip.open(fp,'rt') as f:assert [json.loads(line) for line in f]==census
        validate_probes(probes,selected,counts,plan,fit)
        replay(source,fit,plan,rows,probes,selected,counts)
    else:
        assert not cp.exists() and not fp.exists() and not pp.exists(),'Original timing artifacts already exist'
        operators={};probes=[]
        for gid,representative in sorted(selected.items()):
            key=representative['identity']['loading_mode'],tuple(representative['source_audit']['retained_kernel_names'])
            if key not in operators:operators[key]=operators_for(source,rows,representative['source_audit'])
            try:probes.append(probe(source,fit,plan,rows,operators[key],representative,counts[gid],gid))
            except Exception as error:
                failure_capture(root/'failures',representative['identity'],probes,rows,
                    representative['diagonal'],representative['matrix'],representative['response'],
                    source['factors'][representative['identity']['tree']][rows],error)
                raise
        validate_probes(probes,selected,counts,plan,fit)
        temp=fp.with_suffix('.gz.partial')
        with gzip.open(temp,'xt') as f:
            for record in census:f.write(json.dumps(record,sort_keys=True,allow_nan=False)+'\n')
        with temp.open('rb') as f:os.fsync(f.fileno())
        finish_gzip(temp,fp);atomic(pp,probes)
        with custody_lock:
            for record in probes:review_inputs(root,record,selected[record['group_id']],source,rows)
        saved=dict(stage=stage,cohort_id=cid,census_rows=len(census),eligible_candidates=sum(counts.values()),
            timing_groups=len(probes),census_sha256=sha(fp),probes_sha256=sha(pp))
        atomic(cp,saved)
    assert saved==dict(stage=stage,cohort_id=cid,census_rows=len(census),eligible_candidates=sum(counts.values()),
        timing_groups=len(probes),census_sha256=sha(fp),probes_sha256=sha(pp))
    reviews=[]
    with custody_lock:
        for record in probes:reviews.extend(review_inputs(root,record,selected[record['group_id']],source,rows,read=True))
    source_guard(root,cohort,refs,before,'after')
    source_guard(root,cohort,task_refs,task_before,'submitted-task-after')
    manifest=dict(cohort_id=cid,census_path=str(fp.relative_to(root)),census_sha256=sha(fp),
        probes_path=str(pp.relative_to(root)),probes_sha256=sha(pp),receipt_path=str(cp.relative_to(root)),receipt_sha256=sha(cp))
    planning=estimate(probes,fit);assert planning['measured_candidate_coverage']+planning['unmeasured_review_candidate_coverage']==sum(counts.values())
    process=psutil.Process()
    result=dict(cohort_index=index,cohort_id=cid,manifest=manifest,candidate_rows=len(census),
        eligible_candidates=sum(counts.values()),timing_groups=len(probes),
        source_status_counts=dict(Counter(row['source_disposition'] for row in census)),
        timing_status_counts=dict(Counter(record['status'] for record in probes)),planning=planning,
        review_artifacts={str(path.relative_to(root)):sha(path) for path in reviews},
        worker=dict(pid=process.pid,created=process.create_time(),cmdline=process.cmdline()),
        worker_address_space_limit_bytes=resource.getrlimit(resource.RLIMIT_AS)[0],
        cached_input_array_bindings_checked=len(refs),cached_numeric_inputs_preserved=True,
        submitted_task_array_bindings_checked=len(task_refs),submitted_numeric_inputs_preserved=True,
        scientific_eligibility=False)
    atomic(root/'checkpoints'/(str(index).zfill(5)+('.reader.json' if reader else '.producer.json')),result)
    return result


def run(path,reader=False,output=None):
    path=Path(path);plan=json.loads(path.read_text());root=Path(plan['output'])
    if output is None:output=root/('readback.json' if reader else 'receipt.json')
    output=Path(output);assert not output.exists(),'Completed original timing role cannot restart'
    assert not (root/'reader_completed.json').exists(),'Completed original reader cannot restart'
    lock_path=Path(str(root)+'.parallel-timing.lock');lock_path.parent.mkdir(parents=True,exist_ok=True)
    with lock_path.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        return run_locked(path,plan,root,reader,output)


def run_locked(path,plan,root,reader,output):
    started=time.perf_counter();runtime(plan);r=plan['resources']
    assert 1<=r['workers']<=r['cpus']<=16 and r['worker_address_space_gib']==12
    assert r['workers']*12<=r['reservation_capacity_gib']
    assert r['reservation_capacity_gib']+r['parent_headroom_gib']<=r['memory_gib']
    assert shutil.disk_usage(root.parent).free>=r['minimum_free_disk_gib']*2**30
    assert not (root/'failures').exists(),'Original failed namespace cannot restart'
    fit,source,bindings=sources(plan,path);original=dict(bindings)
    stage=dict(schema=SCHEMA,plan_sha256=sha(path),fit_contract=source['fit_contract'])
    receipt=None;manifest=None
    if reader:
        rp=root/'receipt.json';receipt=json.loads(rp.read_text());bind(bindings,rp)
        assert receipt['status']==PRODUCER and receipt['plan_sha256']==sha(path)
        assert receipt['source_hashes']==original and receipt['fit_contract']==source['fit_contract']
        assert receipt['scientific_eligibility'] is False
        for name,digest in receipt['artifacts'].items():bind(bindings,root/name,digest)
        verify(bindings);manifest=json.loads((root/'cohort_manifest.json').read_text())
        assert [row['cohort_id'] for row in manifest]==[c['cohort_id'] for c in source['cohorts']]
        assert json.loads((root/'stage_plan.json').read_text())==stage
    else:
        assert psutil.virtual_memory().available>=r['memory_gib']*2**30
        root.mkdir(parents=True,exist_ok=False);(root/'cohorts').mkdir();(root/'checkpoints').mkdir()
        atomic(root/'stage_plan.json',stage)
    # Only bounded live tasks contain full X/y/audit arrays; the complete grid
    # is never materialized as an array-bearing task list in parent memory.
    stream=enumerate(cohorts(source,fit));pending={};outcomes={}
    worker_plan=dict(plan,_original_path=str(path))
    context=multiprocessing.get_context('fork');custody=context.Lock()
    with ProcessPoolExecutor(max_workers=r['workers'],mp_context=context,
            initializer=worker_init,initargs=(source,fit,worker_plan,str(root),reader,custody)) as pool:
        def submit():
            item=next(stream,None)
            if item is not None:
                index,(cohort,rows,entries)=item
                pending[pool.submit(cohort_job,(index,cohort,rows,entries))]=index
        for _ in range(r['workers']):submit()
        try:
            while pending:
                done,_=wait(pending,return_when=FIRST_COMPLETED)
                for future in done:
                    index=pending.pop(future);result=future.result()
                    assert result['cohort_index']==index and index not in outcomes
                    outcomes[index]=result
                    print('parallel_weighted_timing_cohorts',len(outcomes),'/',len(source['cohorts']),'reader',reader,flush=True)
                for _ in done:submit()
        except BaseException:
            for future in pending:future.cancel()
            raise
    assert sorted(outcomes)==list(range(len(source['cohorts'])))
    global_planning=estimate([],fit);status=Counter();tstatus=Counter();artifacts={};entries=[]
    total=eligible=timing_groups=0
    for index in sorted(outcomes):
        value=outcomes[index];entry=value['manifest'];entries.append(entry)
        if reader:assert manifest[index]==entry
        total+=value['candidate_rows'];eligible+=value['eligible_candidates'];timing_groups+=value['timing_groups']
        status.update(value['source_status_counts']);tstatus.update(value['timing_status_counts'])
        costs=value['planning']
        for field in ['measured_candidate_coverage','unmeasured_review_candidate_coverage']:
            global_planning[field]+=costs[field]
        global_planning['group_planning_costs'].extend(costs['group_planning_costs'])
        for kind in ['census','probes','receipt']:artifacts[entry[kind+'_path']]=entry[kind+'_sha256']
        artifacts.update(value['review_artifacts'])
        cp=root/'checkpoints'/(str(index).zfill(5)+'.producer.json')
        prior=json.loads(cp.read_text());assert prior['manifest']==entry and prior['cached_numeric_inputs_preserved']
        assert prior['submitted_numeric_inputs_preserved'] is True
        assert prior['worker_address_space_limit_bytes']==12*2**30
        assert prior['cached_input_array_bindings_checked']>0 and prior['submitted_task_array_bindings_checked']>0
        assert prior['cohort_index']==index and prior['cohort_id']==entry['cohort_id']
        assert prior['scientific_eligibility'] is False
        for field in ['candidate_rows','eligible_candidates','timing_groups','source_status_counts','timing_status_counts','planning','review_artifacts']:
            assert prior[field]==value[field]
        artifacts[str(cp.relative_to(root))]=sha(cp)
    global_planning['conditional_budget_weighted_seconds']=sum(
        v['conditional_budget_weighted_producer_seconds']+v['conditional_budget_weighted_reader_seconds']
        for v in global_planning['group_planning_costs'])
    assert total==source['numerical_completion']['designs']*2*len(fit['methods'])*len(fit['loading_modes'])*len(fit['trees'])*len(fit['policies'])==fit['expected']['candidate_rows']
    assert global_planning['measured_candidate_coverage']+global_planning['unmeasured_review_candidate_coverage']==eligible
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=source['numerical_completion']['settings'],
        unique_cohorts=len(entries),candidate_rows=total,eligible_candidates=eligible,timing_groups=timing_groups,
        source_status_counts=dict(status),timing_status_counts=dict(tstatus),
        measured_candidate_coverage=global_planning['measured_candidate_coverage'],
        unmeasured_review_candidate_coverage=global_planning['unmeasured_review_candidate_coverage'],
        conditional_budget_weighted_seconds=global_planning['conditional_budget_weighted_seconds'])
    if reader:
        assert json.loads((root/'conditional_planning.json').read_text())==global_planning
        assert all(receipt[k]==v for k,v in summary.items()) and receipt['production_finish_eta'] is None
    else:
        atomic(root/'cohort_manifest.json',entries);atomic(root/'conditional_planning.json',global_planning)
    for name in ['stage_plan.json','cohort_manifest.json','conditional_planning.json']:artifacts[name]=sha(root/name)
    if reader:assert receipt['artifacts']==artifacts
    verify(bindings)
    result=dict(status=READER if reader else PRODUCER,plan_sha256=sha(path),fit_contract=source['fit_contract'],
        **summary,source_hashes=bindings,artifacts=artifacts,parallel_workers=r['workers'],
        cached_numeric_inputs_preserved_for_all_cohorts=True,worker_address_space_gib=12,
        fits_computed=0,nonuniform_weighting_accepted=False,component_variance_attribution_accepted=False,
        scientific_eligibility=False,production_finish_eta=None,scope=plan['scope'])
    if reader:
        result.update(producer_receipt_sha256=sha(root/'receipt.json'),fresh_numeric_probes_replayed=timing_groups,
            independent_hardware_timing_reimplementation=False,
            reader_checkpoint_artifacts={str(cp.relative_to(root)):sha(cp) for cp in sorted((root/'checkpoints').glob('*.reader.json'))})
    else:result.update(observed_wall_seconds=time.perf_counter()-started,observed_current_rss_bytes=psutil.Process().memory_info().rss)
    atomic(output,result)
    if reader:atomic(root/'reader_completed.json',dict(output=str(output),sha256=sha(output),status=READER,plan_sha256=sha(path)))
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts','reader_checkpoint_artifacts','scope']}),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--reader',action='store_true');p.add_argument('--output',type=Path)
    a=p.parse_args();run(a.plan,a.reader,a.output)
