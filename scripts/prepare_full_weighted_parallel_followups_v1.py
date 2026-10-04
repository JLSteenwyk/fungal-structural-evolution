#!/usr/bin/env python3
"""Create fresh parallel reader/closure dependency plans after exact launches exist."""
import argparse
import json
from pathlib import Path
import sys

from ancestral_chain_attempt import sha
from full_weighted_covariance_qualification import PRODUCER, READER, SUMMARY
from reference_measurement_union_sources import verify


def write(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2);f.write('\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--role',choices=['reader','closure','inventory'],required=True)
    a=p.parse_args();pp=Path('metadata/full_weighted_parallel_qualification_plan_20261004_v1.json')
    plan=json.loads(pp.read_text());verify(plan['pins'])
    producer='metadata/full_weighted_parallel_producer_launch_20261004_v1.json'
    reader='metadata/full_weighted_parallel_reader_launch_20261004_v1.json'
    closer='metadata/full_weighted_parallel_closure_launch_20261004_v1.json'
    if a.role=='inventory':
        launches=[producer,reader,closer]
        records=[json.loads(Path(path).read_text()) for path in launches]
        result=dict(status='launched_full_original_parallel_numerical_producer_with_gated_reader_and_closure',
            source_plan=str(pp),source_plan_sha256=sha(pp),launches=launches,
            original_tool_sessions=[r['original_tool_session_id'] for r in records],expected=plan['expected'],
            scientific_eligibility=False,weighted_working_model_fits_launched=False,gpu=False,
            original_serial_jobs_restarted=False,all_eight_aims_incomplete=True)
        write('metadata/full_weighted_parallel_qualification_launches_20261004_v1.json',result)
        print(json.dumps(result,indent=2));return
    if a.role=='reader':
        deps=[producer];source=pp
        command=['/usr/bin/prlimit','--as='+str(8*2**30)+':'+str(24*2**30),
            '--cpu=604800','--fsize='+str(512*2**20),'--',sys.executable,
            'scripts/full_weighted_covariance_qualification_parallel_v1.py','--plan',str(pp),'--reader']
    else:
        deps=[producer,reader]
        source=Path('metadata/full_weighted_parallel_qualification_completion_plan_20261004_v1.json')
        root=Path(plan['output'])
        write(source,dict(source_plan=str(pp),producer_receipt=str(root/'receipt.json'),
            independent_readback=str(root/'readback.json'),producer_status=PRODUCER,reader_status=READER,
            completed_status='complete_verified_full_four_control_covariance_numerical_qualification_v1',
            summary_fields=SUMMARY+['parallel_workers','worker_address_space_gib','cached_numeric_inputs_preserved_for_all_cohorts'],
            launches=deps,output=plan['completion'],scope=plan['scope'],
            pins={str(path):sha(path) for path in [pp,*map(Path,deps),Path(plan['closure_script']),
                Path('scripts/record_completed_process_handoffs_v2.py'),Path(__file__)]}))
        command=['/usr/bin/prlimit','--as='+str(24*2**30),'--cpu=7200',
            '--fsize='+str(512*2**20),'--',sys.executable,plan['closure_script'],'--plan',str(source)]
    wp=Path('metadata/full_weighted_parallel_'+a.role+'_wait_plan_20261004_v1.json')
    write(wp,dict(dependencies=deps,command=command,
        pins={str(path):sha(path) for path in [source,*map(Path,deps),
            Path('scripts/run_after_verified_dependencies_v2.py'),Path(__file__)]}))
    print(json.dumps(dict(role=a.role,wait_plan=str(wp),command=command),indent=2))


if __name__=='__main__':main()
