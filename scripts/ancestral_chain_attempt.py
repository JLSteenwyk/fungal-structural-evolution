#!/usr/bin/env python3
"""Run isolated chain attempts; exit zero is not posterior qualification.

Callers supply an absolute executable, immutable input hashes, seed and model
settings in config. Commands may contain {attempt} for the new output path.
Incomplete attempts lacking process identity require manual review.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time

import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def check_inputs(config):
    assert Path(config['command'][0]).is_absolute()
    assert config['command'][0] in config['pins'], 'Executable must be pinned'
    assert config['timeout_seconds'] > 0
    assert config['pins']
    for name, digest in config['pins'].items():
        assert Path(name).is_absolute(), name
        assert sha(name) == digest, name


def live_group(group):
    live = []
    for process in psutil.process_iter(['pid', 'status']):
        try:
            if os.getpgid(process.pid) == group and process.status() != psutil.STATUS_ZOMBIE:
                live.append(process.pid)
        except (ProcessLookupError, psutil.NoSuchProcess):
            pass
    return live


def run_attempt(root, config):
    """Reuse checked exit-zero output or create a separate attempt under lock.

    Caller must validate biological outputs separately. No native checkpoint,
    sample concatenation, deletion or automatic adoption of orphaned outputs.
    """
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    check_inputs(config)
    encoded = json.dumps(config, sort_keys=True, separators=(',', ':'), allow_nan=False)
    digest = hashlib.sha256(encoded.encode()).hexdigest()
    with (root / 'attempt.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        binding = root / 'configuration.json'
        if binding.exists():
            assert json.loads(binding.read_text()) == config, 'Configuration changed'
        else:
            write_json(binding, config)
        previous = sorted(root.glob('attempt-[0-9]*'))
        reusable = None
        for folder in previous:
            receipt = folder / 'receipt.json'
            identity = folder / 'process.json'
            if not identity.exists():
                raise RuntimeError('Missing process identity; review incomplete attempt: ' + str(folder))
            info = json.loads(identity.read_text())
            assert not live_group(info['pgid']), 'An earlier process group is still live'
            if psutil.pid_exists(info['pid']):
                process = psutil.Process(info['pid'])
                assert process.create_time() != info['created'] or process.status() == psutil.STATUS_ZOMBIE
            if receipt.exists():
                result = json.loads(receipt.read_text())
                assert result['configuration_sha256'] == digest
                for name, expected in result['artifacts'].items():
                    assert sha(folder / name) == expected, name
                if result['exit_code'] == 0 and result['status'] == 'exited_zero_pending_scientific_validation':
                    reusable = receipt
        if reusable is not None:
            return reusable
        folder = root / ('attempt-%04d' % (len(previous) + 1))
        folder.mkdir(exist_ok=False)
        command = [part.replace('{attempt}', str(folder)) for part in config['command']]
        write_json(folder / 'command.json', command)
        started = time.monotonic()
        with (folder / 'stdout.log').open('w') as stdout, (folder / 'stderr.log').open('w') as stderr:
            process = subprocess.Popen(command, cwd=folder, stdout=stdout, stderr=stderr,
                                       start_new_session=True, pass_fds=(lock.fileno(),))
            identity = psutil.Process(process.pid)
            write_json(folder / 'process.json', dict(pid=process.pid, pgid=process.pid,
                       created=identity.create_time(), command=command))
            timed_out = False
            try:
                code = process.wait(timeout=config['timeout_seconds'])
            except subprocess.TimeoutExpired:
                timed_out = True
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                code = process.wait()
        assert not live_group(process.pid), 'Child descendants still active; receipt withheld'
        check_inputs(config)
        result = dict(configuration_sha256=digest, exit_code=code,
                      status='timeout' if timed_out else ('exited_zero_pending_scientific_validation' if code == 0 else 'failed'),
                      elapsed_seconds=time.monotonic() - started,
                      artifacts={str(p.relative_to(folder)): sha(p) for p in folder.rglob('*') if p.is_file()},
                      scope='Whole new attempt; never continued or concatenated with earlier samples.')
        receipt = folder / 'receipt.json'
        write_json(receipt, result)
        return receipt
