#!/usr/bin/env python3
"""Run candidate-length screens after each scalar quartet report is complete."""
import argparse
import json
from pathlib import Path
import subprocess
import time
import psutil
from ancestral_chain_attempt import sha, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    digest = sha(args.plan)
    def verify():
        assert sha(args.plan) == digest
        for path, expected in plan['pins'].items():
            assert sha(path) == expected, path
    verify()
    launch = json.loads(Path(plan['scalar_launch']).read_text())
    scalar_plan = json.loads(Path(plan['scalar_plan']).read_text())
    assert launch['plan_sha256'] == sha(plan['scalar_plan'])
    producer = json.loads(Path(plan['producer_plan']).read_text())
    jobs = json.loads(Path(producer['jobs']).read_text())
    groups = sorted({job['config']['model_input_identity'] for job in jobs})
    assert len(groups) == 405
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    complete = {}
    while True:
        verify()
        state = dict(line.split('=', 1) for line in subprocess.check_output([
            'systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'MainPID',
            '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        terminal = state['ActiveState'] in ['inactive', 'failed']
        if not terminal:
            assert state['ActiveState'] == 'active' and int(state['MainPID']) == launch['pid']
            try:
                process = psutil.Process(launch['pid'])
                assert process.create_time() == launch['created'] and process.cmdline() == launch['cmdline']
            except psutil.NoSuchProcess:
                time.sleep(1)
                continue
        for group in groups:
            if group in complete:
                continue
            source = Path(scalar_plan['output']) / group / 'receipt.json'
            if not source.exists():
                continue
            report = json.loads(source.read_text())
            assert report['status'] == 'both_burnin_scalar_screens_complete'
            assert set(report['outputs']) == {'250', '500'}
            for path, expected in report['evidence'].items():
                assert sha(path) == expected
            for output in report['outputs'].values():
                assert sha(output['path']) == output['sha256']
            target = out / group
            subprocess.run([plan['runner_python'], plan['script'], '--producer-plan', plan['producer_plan'],
                '--group', group, '--output', str(target), '--python', plan['diagnostic_python']],
                check=True, stdout=subprocess.DEVNULL)
            rp = target / 'receipt.json'
            result = json.loads(rp.read_text())
            assert result['status'] == 'candidate_length_screens_complete_not_posterior_qualification'
            assert result['model_input_identity'] == group
            for name, expected in result['artifacts'].items():
                assert sha(target / name) == expected
            assert {cutoff: value['retained_samples_per_chain'] for cutoff, value in result['outputs'].items()} == {'250': 75, '500': 50}
            complete[group] = dict(status='length_screens_complete_pending_review',
                receipt=str(rp), receipt_sha256=sha(rp), scalar_receipt=str(source), scalar_receipt_sha256=sha(source))
            print('Length-screen quartets', len(complete), '/405', group, flush=True)
        if terminal:
            break
        time.sleep(30)
    write_json(out / 'receipt.json', dict(status='terminal_length_group_accounting_requires_review',
        scalar_service_terminal_state=state, plan_sha256=digest,
        groups={group: complete.get(group, dict(status='unresolved_no_complete_scalar_quartet')) for group in groups},
        complete_quartets=len(complete), unresolved_quartets=405-len(complete),
        scope='Length mixing screens only; no state/homology convergence or ancestral structure qualification.'))


if __name__ == '__main__':
    main()
