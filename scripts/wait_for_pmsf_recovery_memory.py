"""Wait for the reviewed PMSF memory prerequisite, then run checkpoint recovery."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    digest = sha(args.plan)
    plan = json.loads(args.plan.read_text())
    config = json.loads((Path(plan['original_output']) / 'config.json').read_text())
    lock = Path(plan['control'] + '.memory-wait.lock')
    lock.parent.mkdir(parents=True, exist_ok=True)
    handle = lock.open('a')
    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    while True:
        if sha(args.plan) != digest:
            raise ValueError('Recovery plan changed while waiting')
        for name, checksum in plan['pins'].items():
            if sha(name) != checksum:
                raise ValueError('Recovery dependency changed: ' + name)
        if Path(plan['original_output'], 'receipt.json').exists():
            raise RuntimeError('Another execution produced a receipt; inspect before proceeding')
        for process in psutil.process_iter(['cmdline']):
            if process.info['cmdline'] == config['command']:
                raise RuntimeError('Native PMSF is already active')
        available = psutil.virtual_memory().available / 2**30
        required = plan['resources']['minimum_available_memory_gib']
        if available >= required:
            break
        print(f'waiting_for_memory available_gib={available:.2f} required_gib={required}', flush=True)
        time.sleep(30)
    # The existing recovery runner rechecks hashes, checkpoint, memory, disk,
    # duplicate native processes and the global PMSF lock before launching.
    subprocess.run([sys.executable, 'scripts/recover_species_pmsf_checkpoint.py',
                    '--plan', str(args.plan)], check=True)


if __name__ == '__main__':
    main()
