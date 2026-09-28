#!/usr/bin/env python3
"""Audit every paired FastML output separately, preserving original strict failures."""
import argparse
import json
from pathlib import Path
import subprocess
import time
import psutil
from ancestral_chain_attempt import sha, write_json
from readback_fastml_variant_roundoff import readback


def producer_state(launch):
    state = dict(x.split('=', 1) for x in subprocess.check_output(
        ['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState',
         '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    if state['ActiveState'] in ['active', 'activating']:
        try:
            p = psutil.Process(launch['pid'])
            assert p.create_time() == launch['created'] and p.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            # The service can become terminal between state and process reads.
            return None
        return None
    assert state['ActiveState'] in ['inactive', 'failed']
    return state


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    for p, h in plan['pins'].items():
        assert sha(p) == h
    producer = json.loads(Path(plan['producer_plan']).read_text())
    launch = json.loads(Path(plan['producer_launch']).read_text())
    assert launch['plan_sha256'] == sha(plan['producer_plan'])
    for p, h in producer['pins'].items():
        assert sha(p) == h
    root = Path(producer['output'])
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    completed = {}
    while True:
        terminal = producer_state(launch)
        for item in producer['jobs']:
            if item['id'] in completed:
                continue
            source = root / item['id'] / 'readback.json'
            if not source.exists():
                continue
            original = json.loads(source.read_text())
            assert original['job'] == item['job'] and original['variant'] == item['variant']
            evidence = {str(source): sha(source)}
            if not item['job']['character_count']:
                assert original['status'] == 'no_coded_characters'
                result = dict(status='no_coded_characters')
            else:
                rp = Path(original['attempt_receipt'])
                assert rp.resolve().is_relative_to((root / item['id'] / 'attempts').resolve())
                assert sha(rp) == original['attempt_receipt_sha256']
                attempt = json.loads(rp.read_text())
                config = rp.parent.parent / 'configuration.json'
                assert json.loads(config.read_text()) == item['config']
                command = [x.replace('{attempt}', str(rp.parent)) for x in item['config']['command']]
                assert json.loads((rp.parent / 'command.json').read_text()) == command
                for name, digest in attempt['artifacts'].items():
                    assert sha(rp.parent / name) == digest
                evidence.update({str(rp): sha(rp), str(config): sha(config)})
                if attempt['exit_code'] != 0:
                    result = dict(status='inference_failed', exit_code=attempt['exit_code'])
                else:
                    try:
                        result = readback(rp.parent, item['job'])
                    except Exception as error:
                        result = dict(status='roundoff_readback_failed', error=repr(error))
            target = out / (item['id'] + '.json')
            result.update(job_id=item['job']['job_id'], variant=item['variant'], evidence=evidence,
                          original_readback_status=original['status'])
            write_json(target, result)
            completed[item['id']] = dict(path=str(target), sha256=sha(target), status=result['status'])
            print(item['id'], result['status'], flush=True)
        if terminal is not None:
            break
        time.sleep(30)
    for p, h in plan['pins'].items():
        assert sha(p) == h
    for p, h in producer['pins'].items():
        assert sha(p) == h
    write_json(out / 'receipt.json', dict(status='terminal_paired_readback_accounting_requires_review',
        plan_sha256=sha(args.plan), producer_terminal_state=terminal, readbacks=completed,
        unresolved=[j['id'] for j in producer['jobs'] if j['id'] not in completed],
        scope='Independent replay with explicit boundary excursions; original raw values and strict '
              'failure records preserved. All failed/unresolved entries retained. No optimizer '
              'or posterior ensemble qualification.'))


if __name__ == '__main__':
    main()
