#!/usr/bin/env python3
"""Observe the original full reader after its producer has closed, without rerunning either."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

import psutil

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    launch=Path('metadata/full_weighted_parallel_reader_launch_20261004_v1.json')
    previous=Path('metadata/full_weighted_parallel_reader_checkpoint_20261004_goal_1306.json')
    original=json.loads(previous.read_text())
    handle=observe(str(launch),{'cpu.max':'1600000 100000','memory.max':str(200*2**30),'memory.swap.max':'0'})
    native_record=original['native_parent'];parent=psutil.Process(native_record['pid'])
    assert parent.create_time()==native_record['created'] and parent.cmdline()==native_record['cmdline']
    assert parent.rlimit(psutil.RLIMIT_AS)==(8*2**30,24*2**30)
    workers=[]
    for child in parent.children():
        created=child.create_time();command=child.cmdline()
        assert command==native_record['cmdline'] and child.rlimit(psutil.RLIMIT_AS)==(12*2**30,)*2
        assert child.create_time()==created
        workers.append(dict(pid=child.pid,created=created,status=child.status(),cmdline=command,
                            address_space_limits=child.rlimit(psutil.RLIMIT_AS),
                            cpu_seconds=child.cpu_times().user+child.cpu_times().system))
    assert len(workers)==16
    root=Path('results/phylogeny/full-four-control-covariance-parallel-qualification-20261004-v1')
    checkpoints=len(list((root/'checkpoints').glob('*.reader.json')))
    result=dict(status='verified_original_full_parallel_numerical_reader_after_producer_terminal',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_tool_session_id=75216,
        original_handle=handle,native_parent=native_record,workers=workers,
        completed_reader_checkpoints=checkpoints,expected_cohorts=4340,
        independent_readback_present=(root/'readback.json').exists(),
        full_numerical_completion_present=Path('metadata/full_weighted_parallel_qualification_completed_20261004_v1.json').exists(),
        source_hashes={str(q):sha(q) for q in [Path(__file__),launch,previous]},
        scientific_eligibility=False,biological_fits=0,gpu=False,
        scope='Exact original wrapper/native identities,16actual workers and caps; checkpoint census '
              'is progress only, not independently accepted full arithmetic or biological fitting.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['original_handle','workers','native_parent','source_hashes']},indent=2))


if __name__=='__main__':main()
