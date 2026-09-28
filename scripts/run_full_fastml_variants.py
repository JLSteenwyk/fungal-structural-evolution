#!/usr/bin/env python3
"""Run all paired FastML variants with recoverable attempts and per-fit replay."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import json
from pathlib import Path
import shutil
from ancestral_chain_attempt import run_attempt, sha, write_json
from readback_fastml_variant import readback
from replay_fastml_indel_probabilities import analytic_check


def execute(item, root):
    destination = root / item['id']
    destination.mkdir(parents=True, exist_ok=True)
    final = destination / 'readback.json'
    if final.exists():
        prior = json.loads(final.read_text())
        if prior['status'] == 'no_coded_characters':
            assert item['job']['character_count'] == 0 and prior['job'] == item['job']
            return str(final)
        # run_attempt verifies configuration, attempt identity and every artifact
        # before any successful result is reused; replay is recomputed below.
    if not item['job']['character_count']:
        write_json(final, dict(status='no_coded_characters', job=item['job'], variant=item['variant']))
        return str(final)
    rp = run_attempt(destination / 'attempts', item['config'])
    receipt = json.loads(rp.read_text())
    result = dict(status='execution_failed', exit_code=receipt['exit_code'])
    if receipt['exit_code'] == 0:
        try:
            result = readback(rp.parent, item['job'])
        except Exception as error:
            result = dict(status='readback_failed', error=repr(error))
    result.update(job=item['job'], variant=item['variant'], attempt_receipt=str(rp),
                  attempt_receipt_sha256=sha(rp))
    write_json(final, result)
    return str(final)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    for p, h in plan['pins'].items():
        assert sha(p) == h
    assert len(plan['jobs']) == 312 and len({j['id'] for j in plan['jobs']}) == 312
    assert sum(j['job']['character_count'] > 0 for j in plan['jobs']) == 306
    assert shutil.disk_usage('.').free > plan['minimum_free_disk_bytes']
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        binding = out / 'plan_sha256.txt'
        if binding.exists():
            assert binding.read_text().strip() == sha(args.plan)
        else:
            binding.write_text(sha(args.plan) + '\n')
        analytic_check()
        completed = []
        with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
            tasks = {pool.submit(execute, item, out): item['id'] for item in plan['jobs']}
            for future in as_completed(tasks):
                path = future.result()
                completed.append(path)
                print(tasks[future], json.loads(Path(path).read_text())['status'], flush=True)
        assert len(completed) == 312
        for p, h in plan['pins'].items():
            assert sha(p) == h
        write_json(out / 'receipt.json', dict(status='all_paired_inputs_accounted_requires_scientific_review',
            plan_sha256=sha(args.plan), readbacks={p: sha(p) for p in sorted(completed)},
            analytic_enumeration='passed', scope='312 variant/input records, including six empty cases. '
            'Failures retained explicitly. Completion does not certify optimization or model adequacy.'))


if __name__ == '__main__':
    main()
