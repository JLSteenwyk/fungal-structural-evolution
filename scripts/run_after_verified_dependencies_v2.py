#!/usr/bin/env python3
"""Run pinned stages only after exact original dependency completion journals.

Collected transient units can expose inactive/success/0 after failures.
Require invocation-linked original process messages and completion resources;
never treat those default unit fields as proof of successful completion.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time

from run_after_verified_dependencies import live,terminal_state
from run_ortholog_pair_guide_comparison import sha


def inspect_journal(launch,rows):
    exact=[r for r in rows if r.get('_PID')==str(launch['pid']) and
           r.get('_CMDLINE')==' '.join(launch['cmdline'])]
    assert exact,('No original process journal',launch['unit'])
    invocations={r['_SYSTEMD_INVOCATION_ID'] for r in exact}
    resource=[r for r in rows if r.get('USER_INVOCATION_ID') in invocations and r.get('CPU_USAGE_NSEC')]
    assert resource,('No original completion resource record',launch['unit'])
    failed=[r for r in rows if r.get('USER_INVOCATION_ID') in invocations and
            ('Main process exited' in (r.get('MESSAGE') or '') or
             'Failed with result' in (r.get('MESSAGE') or ''))]
    if failed:raise RuntimeError('Original dependency invocation failed: '+launch['unit'])
    return dict(original_process_messages=len(exact),original_completion_resource_records=len(resource))


def completed(launch):
    assert not live(launch)
    state=terminal_state(launch['unit'])
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),state
    journal=subprocess.check_output(['journalctl','--user','-u',launch['unit'],'-o','json','--no-pager'],text=True)
    return inspect_journal(launch,[json.loads(line) for line in journal.splitlines()])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());bindings={str(args.plan):sha(args.plan),**plan['pins']}
    dependencies=[]
    for path in plan['dependencies']:
        assert path in bindings
        launch=json.loads(Path(path).read_text());assert sha(launch['plan'])==launch['plan_sha256']
        if launch['plan'] in bindings:assert bindings[launch['plan']]==launch['plan_sha256']
        bindings[launch['plan']]=launch['plan_sha256'];dependencies.append(launch)
    def verify():
        for path,digest in bindings.items():assert sha(path)==digest,('Changed dependency',path)
    verify();remaining=dependencies.copy()
    while remaining:
        pending=[]
        for launch in remaining:
            if live(launch):pending.append(launch);continue
            state=terminal_state(launch['unit'])
            if state['ActiveState'] in ['active','activating','deactivating','reloading']:
                pending.append(launch);continue
            completed(launch)
        remaining=pending
        if remaining:
            print('waiting_for_original_dependency_completion',[r['unit'] for r in remaining],flush=True);time.sleep(30)
    verify()
    for launch in dependencies:completed(launch)
    cmd=plan['command'];assert isinstance(cmd,list) and cmd and all(isinstance(x,str) for x in cmd)
    print('running_original_journal_verified_command',cmd,flush=True)
    subprocess.run(cmd,check=True);verify()


if __name__=='__main__':main()
