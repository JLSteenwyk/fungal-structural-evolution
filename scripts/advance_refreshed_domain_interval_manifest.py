"""Wait for full registry validation, then prepare and audit every domain interval."""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
import psutil
from readback_whole_proteome_catalog import sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--launch', type=Path, required=True)
    args = parser.parse_args()
    pins = {p: sha(p) for p in (args.plan, args.launch)}
    plan = json.loads(args.plan.read_text())
    launch = json.loads(args.launch.read_text())
    assert sha(plan['producer_plan']) == launch['plan_sha256']
    while True:
        for path, expected in pins.items():
            assert sha(path) == expected
        try:
            proc = psutil.Process(launch['pid'])
            live = proc.create_time() == launch['created'] and proc.status() != psutil.STATUS_ZOMBIE
            if live:
                assert proc.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            live = False
        text = subprocess.check_output(['systemctl', '--user', 'show', launch['unit'],
                                        '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True)
        state = dict(line.split('=', 1) for line in text.splitlines())
        if state['ActiveState'] == 'failed':
            raise RuntimeError('Domain registry failed: ' + str(state))
        if not live and state['ActiveState'] == 'inactive':
            assert state['Result'] == 'success' and state['ExecMainStatus'] == '0'
            break
        print('Waiting for identified domain registry to finish successfully', flush=True)
        time.sleep(20)
    for name, digest in plan['pins'].items():
        assert sha(name) == digest
    subprocess.run([sys.executable, 'scripts/prepare_domain_extraction_manifest_v2.py',
                    '--registry', plan['registry'], '--audit', plan['registry_audit'],
                    '--output', plan['output']], check=True)
    subprocess.run([sys.executable, 'scripts/readback_domain_extraction_manifest_v2.py',
                    '--manifest', plan['output'], '--database', str(Path(plan['registry']) / 'structure_domains.sqlite'),
                    '--output', plan['readback_output']], check=True)


if __name__ == '__main__':
    main()
