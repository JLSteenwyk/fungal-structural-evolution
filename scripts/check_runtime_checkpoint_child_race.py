#!/usr/bin/env python3
"""Verify child-observation races without relaxing original-parent identity checks."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from unittest.mock import Mock, patch

import psutil
import record_project_runtime_checkpoint_v3 as old
import record_project_runtime_checkpoint_v4 as new


def process(pid):
    result=Mock()
    result.pid=pid
    result.create_time.return_value=123.0
    result.cmdline.return_value=['worker',str(pid)]
    result.status.return_value=psutil.STATUS_RUNNING
    result.cpu_times.return_value=(2.0,3.0)
    result.children.return_value=[]
    return result


def rejected(function,error):
    try:
        function()
    except error:
        return
    raise AssertionError('Expected exception was swallowed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    parent=process(100)
    survivor=process(101)
    disappeared=process(102)
    disappeared.cmdline.side_effect=psutil.NoSuchProcess(102)
    parent.children.return_value=[disappeared,survivor]
    record=dict(pid=100,created=123.0,cmdline=['worker','100'])
    rejected(lambda:old.live_record(record,parent),psutil.NoSuchProcess)
    observed=new.live_record(record,parent)
    assert observed['children']==[dict(pid=101,created=123.0,cmdline=['worker','101'],
        status=psutil.STATUS_RUNNING,cpu_seconds_at_observation=5.0)]
    assert observed['children_unavailable_during_observation'][0]['pid']==102
    assert observed['children_unavailable_during_observation'][0]['error_type']=='NoSuchProcess'
    zombie=process(103)
    zombie.cpu_times.side_effect=psutil.ZombieProcess(103)
    parent.children.return_value=[zombie]
    assert new.live_record(record,parent)['children_unavailable_during_observation'][0]['error_type']=='ZombieProcess'
    denied=process(104)
    denied.cmdline.side_effect=psutil.AccessDenied(104)
    parent.children.return_value=[denied]
    rejected(lambda:new.live_record(record,parent),psutil.AccessDenied)
    parent.children.return_value=[]
    parent.cpu_times.side_effect=psutil.NoSuchProcess(100)
    rejected(lambda:new.live_record(record,parent),psutil.NoSuchProcess)
    reused=process(100)
    reused.create_time.return_value=124.0
    with patch.object(new.psutil,'Process',return_value=reused):
        assert new.fingerprint(record) is None
    wrong=process(100)
    wrong.cmdline.return_value=['different-worker']
    with patch.object(new.psutil,'Process',return_value=wrong):
        rejected(lambda:new.fingerprint(record),AssertionError)
    sources=[Path(__file__),Path(old.__file__),Path(new.__file__)]
    receipt=dict(status='passed_runtime_checkpoint_child_race_checks',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        checks=['original_v3_race_reproduced','disappeared_child_recorded_survivor_preserved',
            'zombie_child_recorded','child_access_denied_remains_fatal',
            'parent_disappearance_remains_fatal','reused_parent_pid_not_accepted',
            'changed_parent_command_rejected'],
        source_hashes={str(p.relative_to(Path.cwd())):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        scope='Software process-observation checks only; real scientific processes are neither stopped nor restarted.')
    with args.output.open('x') as output:
        output.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))


if __name__=='__main__':
    main()
