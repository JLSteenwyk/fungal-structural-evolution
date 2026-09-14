#!/usr/bin/env python3
"""Time-limited resource control for explicitly recorded live process trees.

No global settings or analysis configuration changes. SIGTERM restores affinity
and resumes any GPU process stopped by this controller. --restore recovers from
an abrupt controller death using the persisted state. Linux, same-user jobs only.
"""
import argparse
import json
import os
from pathlib import Path
import signal
import time
import psutil


def save(path, data):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    tmp.replace(path)


def process(pid, created):
    try:
        p = psutil.Process(int(pid))
        return p if p.create_time() == created and p.uids().real == os.getuid() else None
    except psutil.NoSuchProcess:
        return None


def restore(state):
    # Resume first, even if affinity restoration encounters an error.
    gpu = state.get('stopped_gpu')
    if gpu:
        p = process(gpu['pid'], gpu['created'])
        if p:
            p.send_signal(signal.SIGCONT)
        state['stopped_gpu'] = None
    errors = []
    for row in state['processes'].values():
        p = process(row['pid'], row['created'])
        if not p:
            continue
        originals = row['threads']
        try:
            for thread in p.threads():
                # Threads spawned while limited inherit the leader's old mask
                # for restoration; pre-existing threads retain their own mask.
                old = originals.get(str(thread.id), row['original_affinity'])
                try:
                    os.sched_setaffinity(thread.id, old)
                except ProcessLookupError:
                    pass
        except (psutil.NoSuchProcess, ProcessLookupError):
            pass
        except (PermissionError, OSError) as exc:
            errors.append({'pid': row['pid'], 'error': str(exc)})
    state['restore_errors'] = errors
    state['status'] = 'restored' if not errors else 'restore_errors'
    state['ended_at'] = time.time()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--state', type=Path, required=True)
    parser.add_argument('--restore', action='store_true')
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    if args.restore:
        state = json.loads(args.state.read_text())
        live = process(state['controller_pid'], state['controller_created'])
        if live:
            raise RuntimeError('Controller still alive: send it SIGTERM to restore')
        restore(state)
        save(args.state, state)
        return
    if args.state.exists():
        raise FileExistsError(args.state)
    state = {'controller_pid': os.getpid(), 'controller_created': psutil.Process().create_time(),
             'status': 'active', 'started_at': time.time(), 'processes': {}, 'stopped_gpu': None}
    stop = False

    def finish(*_):
        nonlocal stop
        stop = True

    signal.signal(signal.SIGTERM, finish)
    signal.signal(signal.SIGINT, finish)
    started = time.monotonic()
    save(args.state, state)
    try:
        while not stop and time.monotonic() - started < config['duration_seconds']:
            targets = {}
            for root in config['roots'] + list(state['processes'].values()):
                p = process(root['pid'], root['created'])
                if not p:
                    continue
                try:
                    for child in [p] + p.children(recursive=True):
                        if child.uids().real == os.getuid() and child.pid != os.getpid():
                            targets[child.pid] = child
                except psutil.NoSuchProcess:
                    pass
            if not targets:
                break
            for pid, p in targets.items():
                try:
                    key = str(pid) + ':' + str(p.create_time())
                    if key not in state['processes']:
                        old = p.cpu_affinity()
                        # A new child may already have inherited the limit.
                        if set(old) == set(config['cpus']):
                            old = config['original_affinity']
                        state['processes'][key] = {'pid': pid, 'created': p.create_time(),
                                                  'original_affinity': old, 'threads': {}}
                    row = state['processes'][key]
                    pending = []
                    for thread in p.threads():
                        try:
                            affinity = sorted(os.sched_getaffinity(thread.id))
                            if str(thread.id) not in row['threads']:
                                row['threads'][str(thread.id)] = (row['original_affinity']
                                    if affinity == config['cpus'] else affinity)
                            pending.append(thread.id)
                        except ProcessLookupError:
                            pass
                    # Persist restoration information before changing anything.
                    save(args.state, state)
                    for tid in pending:
                        try:
                            os.sched_setaffinity(tid, config['cpus'])
                        except ProcessLookupError:
                            pass
                except (psutil.NoSuchProcess, ProcessLookupError):
                    pass
            gpu = process(config['gpu']['pid'], config['gpu']['created'])
            phase = (time.monotonic() - started) % (config['run_seconds'] + config['rest_seconds'])
            resting = phase >= config['run_seconds']
            if gpu and resting and not state['stopped_gpu']:
                if gpu.status() == psutil.STATUS_STOPPED:
                    raise RuntimeError('GPU process independently stopped; refusing to take ownership')
                state['stopped_gpu'] = config['gpu']
                save(args.state, state)
                gpu.send_signal(signal.SIGSTOP)
            elif state['stopped_gpu'] and (not resting or not gpu):
                if gpu:
                    gpu.send_signal(signal.SIGCONT)
                state['stopped_gpu'] = None
            state['last_update'] = time.time()
            save(args.state, state)
            time.sleep(2)
    finally:
        restore(state)
        save(args.state, state)


if __name__ == '__main__':
    main()
