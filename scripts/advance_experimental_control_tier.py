#!/usr/bin/env python3
"""Launch the next control tier only after the pinned predecessor completes."""
import argparse
import fcntl
import json
from pathlib import Path
import subprocess
import sys
import time
import psutil
from run_experimental_control_tier import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    a = p.parse_args()
    config_hash = sha(a.config)
    c = json.loads(a.config.read_text())
    status = Path(c['status_path'])
    lock = status.with_suffix('.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if status.exists():
        raise FileExistsError('Review prior controller state before retrying')
    def verify():
        if sha(a.config) != config_hash:
            raise ValueError('Controller configuration changed')
        for name, digest in c['pins'].items():
            if sha(name) != digest:
                raise ValueError('Changed handoff dependency: ' + name)
    def save(state, **extra):
        temp = status.with_suffix('.partial')
        temp.write_text(json.dumps(dict(status=state, updated_at=time.time(),
            config_sha256=config_hash, **extra), indent=2) + '\n')
        temp.replace(status)
    try:
        verify()
        save('waiting_for_pinned_predecessor')
        while True:
            try:
                process = psutil.Process(c['predecessor_pid'])
                live = (process.create_time() == c['predecessor_created']
                        and process.status() != psutil.STATUS_ZOMBIE)
            except psutil.NoSuchProcess:
                live = False
            if not live:
                break
            time.sleep(20)
        verify()
        previous = json.loads(Path(c['predecessor_receipt']).read_text())
        if (previous['status'] != 'complete_prediction_audit_conversion'
                or previous['plan_sha256'] != c['predecessor_plan_sha256']):
            raise ValueError('Predecessor did not complete prediction, audit and conversion')
        save('launching_next_tier', predecessor_receipt_sha256=sha(c['predecessor_receipt']))
        subprocess.run([sys.executable, 'scripts/run_experimental_control_tier.py',
                        '--plan', c['next_plan']], check=True)
        save('complete_next_tier')
    except Exception as exc:
        save('failed_requires_review', error=repr(exc))
        raise


if __name__ == '__main__':
    main()
