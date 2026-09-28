"""Wait for the identified full atom audit, then build and verify the refreshed domain search database."""
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
    assert sha(plan['audit_plan']) == launch['plan_sha256']
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
            raise RuntimeError('Domain coordinate audit failed: ' + str(state))
        if not live and state['ActiveState'] == 'inactive':
            assert state['Result'] == 'success' and state['ExecMainStatus'] == '0'
            break
        print('Waiting for identified full domain atom audit to finish successfully', flush=True)
        time.sleep(20)
    for name, digest in plan['pins'].items():
        assert sha(name) == digest
    subprocess.run([sys.executable, 'scripts/advance_domain_search_database.py',
                    '--plan', str(args.plan)], check=True)


if __name__ == '__main__':
    main()
