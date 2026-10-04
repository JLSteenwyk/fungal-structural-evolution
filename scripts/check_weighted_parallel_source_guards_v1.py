#!/usr/bin/env python3
"""Exercise cached-input guards on private fork memory before and after a cohort."""
import argparse
from datetime import datetime, timezone
import json
import multiprocessing
from pathlib import Path

import numpy as np
import psutil

from ancestral_chain_attempt import sha
import full_weighted_covariance_qualification_parallel_v1 as target
from reference_measurement_union_sources import bind, verify


def child(source,plan,root,task,case):
    root=Path(root);root.mkdir();(root/'cohorts').mkdir();(root/'checkpoints').mkdir()
    target.worker_init(source,plan,str(root),False)
    factor=source['factors'][next(iter(source['factors']))]
    before=target.array_bytes(factor)
    def mutate():
        raw=factor.view(np.uint8);raw.flat[0]=int(raw.flat[0])^1
    if case=='before':mutate()
    elif case=='after':
        original=target.sha
        def instrumented_sha(path):
            result=original(path)
            if str(path).endswith('.links.tsv.gz'):mutate()
            return result
        target.sha=instrumented_sha
    try:
        result=target.cohort_job(task)
    except AssertionError as error:
        expected='Worker cached numeric source changed '+('before cohort' if case=='before' else 'during cohort')
        assert case!='control' and str(error)==expected,(case,str(error))
        result=dict(rejected=True,error=str(error))
    else:
        assert case=='control' and result['cached_numeric_inputs_preserved']
        result=dict(rejected=False,cohort_checkpoint=result)
    after=target.array_bytes(factor)
    assert (before==after)==(case=='control')
    if case=='after':
        failures=list((root/'failures').glob('*/failure.json'));assert len(failures)==1
        assert list(failures[0].parent.glob('*.npy'))
    process=psutil.Process()
    target.write(root/'guard.json',dict(case=case,**result,worker=dict(pid=process.pid,
        created=process.create_time(),cmdline=process.cmdline()),before=before,after=after,
        injected_into_private_fork_memory=True,source_files_modified=False,
        scientific_eligibility=False))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists();a.output.mkdir(exist_ok=False)
    plan=json.loads(a.plan.read_text());target.runtime(plan)
    source,pins=target.source_inputs(plan,a.plan);bind(pins,Path(__file__))
    groups,_=target.settings_by_cohort(source);stream=target.jsonl(source['root']/'unique_designs.jsonl')
    task=(0,source['cohorts'][0],source['controls'][0],[next(stream) for _ in range(30)],groups[source['cohorts'][0]['cohort_id']])
    parent_before={name:target.array_bytes(value) for name,value in target.arrays(source)}
    cases=[]
    for case in ['control','before','after']:
        root=a.output/case
        process=multiprocessing.get_context('fork').Process(target=child,args=(source,plan,str(root),task,case))
        process.start();identity=psutil.Process(process.pid)
        created=identity.create_time();process.join(120)
        assert not process.is_alive(),'Private source-guard worker exceeded observation budget'
        assert process.exitcode==0,(case,process.exitcode)
        report=json.loads((root/'guard.json').read_text())
        assert report['worker']['pid']==process.pid and report['worker']['created']==created
        assert report['rejected']==(case!='control');cases.append(report)
        for path in root.rglob('*'):
            if path.is_file():bind(pins,path)
    assert {name:target.array_bytes(value) for name,value in target.arrays(source)}==parent_before
    verify(pins)
    result=dict(status='passed_private_fork_cached_source_mutation_guards_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),cases=cases,
        parent_numeric_source_preserved=True,original_input_files_preserved=True,
        source_hashes=pins,scientific_eligibility=False,
        scope='One unchanged actual cohort calculation and two deliberate cached-factor '
              'byte changes in private fork memory. Guards reject changes before admission '
              'and after arithmetic, retaining post-cohort mutated arrays. Source files and '
              'parent memory are preserved. Artificial controls do not identify or repair '
              'the prior original native memory corruption.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','cases']},indent=2))


if __name__=='__main__':main()
