#!/usr/bin/env python3
"""Suspend recorded same-user process trees, then resume them with cooling limits."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import time
import psutil


def save(path, value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def live(row):
    try:
        p=psutil.Process(row['pid'])
        return p if p.create_time()==row['created'] and p.uids().real==os.getuid() else None
    except psutil.NoSuchProcess:return None


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=['pause','resume'])
    ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--state',type=Path,required=True)
    ap.add_argument('--cooling-config',type=Path)
    ap.add_argument('--cooling-state',type=Path)
    ap.add_argument('--cooling-script',type=Path)
    args=ap.parse_args();manifest=json.loads(args.manifest.read_text())
    if args.action=='pause':
        if args.state.exists():raise FileExistsError(args.state)
        state={'status':'pausing','started_at':time.time(),'processes':{},'manifest_sha256':hashlib.sha256(args.manifest.read_bytes()).hexdigest()}
        save(args.state,state)
        roots=manifest['roots']
        for _ in range(3):
            for row in roots+list(state['processes'].values()):
                p=live(row)
                if not p:continue
                # Stop the launcher before discovering its descendants.
                todo=[p]
                while todo:
                    q=todo.pop(0)
                    try:
                        if q.pid==os.getpid() or q.pid in manifest['excluded_pids']:raise RuntimeError('Control process included')
                        key=str(q.pid)+':'+str(q.create_time())
                        if key not in state['processes']:
                            state['processes'][key]={'pid':q.pid,'created':q.create_time(),'command':q.cmdline(),'was_stopped':q.status()==psutil.STATUS_STOPPED}
                            save(args.state,state)
                        q.send_signal(signal.SIGSTOP)
                        todo.extend(q.children())
                    except psutil.NoSuchProcess:pass
            time.sleep(.2)
        for row in state['processes'].values():
            p=live(row)
            if p:assert p.status() in [psutil.STATUS_STOPPED,psutil.STATUS_ZOMBIE],p.pid
        state['status']='paused';state['paused_at']=time.time();save(args.state,state)
        print('Paused',len(state['processes']),'processes',flush=True)
    else:
        state=json.loads(args.state.read_text());assert state['status']=='paused'
        assert hashlib.sha256(args.manifest.read_bytes()).hexdigest()==state['manifest_sha256']
        cooling=json.loads(args.cooling_config.read_text()) if args.cooling_config else None
        resumed=[];gone=[]
        if cooling:
            assert args.cooling_script and args.cooling_state and not args.cooling_state.exists()
            assert hashlib.sha256(args.cooling_script.read_bytes()).hexdigest()==manifest['cooling_script_sha256']
            # Reapply the CPU allowance before resuming any worker.
            for row in state['processes'].values():
                p=live(row)
                if not p:continue
                for t in p.threads():
                    try:os.sched_setaffinity(t.id,cooling['cpus'])
                    except ProcessLookupError:pass
            cooling['roots']=list(state['processes'].values());save(args.cooling_config,cooling)
        for row in state['processes'].values():
            p=live(row)
            if not p:gone.append(row['pid']);continue
            if not row['was_stopped']:p.send_signal(signal.SIGCONT);resumed.append(p.pid)
        state.update(status='resumed',resumed_at=time.time(),resumed_pids=resumed,no_longer_present_pids=gone)
        save(args.state,state);print('Resumed',len(resumed),'processes',flush=True)
        if cooling:
            os.execv(manifest['python'],[manifest['python'],str(args.cooling_script),'--config',str(args.cooling_config),'--state',str(args.cooling_state)])


if __name__=='__main__':main()
