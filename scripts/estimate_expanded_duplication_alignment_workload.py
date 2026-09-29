#!/usr/bin/env python3
"""Bin expanded masked inputs against historical task wall times, without reuse."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import subprocess
import time
import psutil
from screen_duplication_alignment_reuse import load, sha


def length_bin(length):
    assert isinstance(length, int) and length >= 3
    return next((str(n) for n in (250, 500, 1000, 2000, 4000) if length <= n), '>4000')


def summarize(counts, timing, workers):
    rows = []
    covered = uncovered = unavailable = 0
    contributions = []
    for (mask, status, bin_name), count in sorted(counts.items()):
        old = timing.get((mask, bin_name)) if status == 'ready' else None
        row = dict(mask=mask, status=status, length_bin=bin_name, directed_dispositions=count)
        if status != 'ready':
            unavailable += count
            row['cost_status'] = 'no_native_call'
        elif old is None:
            uncovered += count
            row['cost_status'] = 'no_historical_length_bin'
        else:
            covered += count
            mean = old['summed_elapsed_seconds'] / old['records']
            contributions.append(count * mean)
            row.update(cost_status='historical_bin_mean_scenario', historical_records=old['records'],
                       historical_mean_seconds=mean, summed_task_seconds=count * mean)
        rows.append(row)
    seconds = math.fsum(contributions)
    return dict(strata=rows, covered_native_calls=covered, uncovered_native_calls=uncovered,
                unavailable_dispositions=unavailable,
                covered_only_scenarios=[dict(workers=w, cost_multiplier=m,
                    idealized_wall_hours=seconds * m / w / 3600)
                    for w in workers for m in (1, 4)],
                covered_summed_task_hours=seconds / 3600)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    launch = json.loads(Path(plan['audit_launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            process = psutil.Process(launch['pid'])
            if abs(process.create_time() - launch['created']) > .01 or process.status() == psutil.STATUS_ZOMBIE:
                break
            assert process.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(line.split('=', 1) for line in subprocess.check_output([
        'systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p',
        'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    root = Path(plan['inputs'])
    proof = json.loads(Path(plan['audit_proof']).read_text())
    assert proof['status'] == 'passed_full_duplication_alignment_input_readback'
    assert proof['source_receipt_sha256'] == sha(root / 'receipt.json')
    bindings[plan['audit_proof']] = sha(plan['audit_proof'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_duplication_alignment_input_materialization'
    models, pairs, queue_bindings = load(Path(plan['queue']))
    bindings.update(queue_bindings)
    assert receipt['queue_receipt_sha256'] == sha(Path(plan['queue']) / 'receipt.json')
    bindings[str(root / 'receipt.json')] = sha(root / 'receipt.json')
    bindings[str(root / 'inputs.jsonl')] = receipt['artifacts']['inputs.jsonl']
    verify()
    inputs = {}
    input_counts = Counter()
    for line in (root / 'inputs.jsonl').open():
        row = json.loads(line)
        key = row['model_id'], row['version'], row['mask']
        assert key not in inputs and key[:2] in models and key[2] in ('full', 'plddt70')
        assert row['source_sha256'] == models[key[:2]]['sha256']
        assert row['status'] in ('ready', 'too_few_retained_residues', 'source_rejected')
        if row['status'] == 'ready':
            assert row['retained_residues'] == len(row['sequence']) >= 3
        inputs[key] = row['status'], row.get('retained_residues', 0)
        input_counts[row['mask'] + ':' + row['status']] += 1
    assert len(inputs) == 2 * len(models) == receipt['input_dispositions']
    assert dict(input_counts) == receipt['counts']
    costs = json.loads(Path(plan['costs']).read_text())
    assert costs['status'] == 'complete_recorded_duplication_alignment_cost_census'
    timing = {(r['mask'], r['maximum_length_upper_bin']): r for r in costs['strata'] if r['status'] == 'aligned'}
    counts = Counter()
    for ends in pairs.values():
        for mask in ('full', 'plddt70'):
            rows = [inputs[(*end, mask)] for end in ends]
            ready = all(r[0] == 'ready' for r in rows)
            counts[(mask, 'ready' if ready else 'unavailable',
                    length_bin(max(r[1] for r in rows)) if ready else 'not_aligned')] += 2
    assert sum(counts.values()) == 4 * len(pairs)
    result = summarize(counts, timing, plan['worker_scenarios'])
    result.update(status='complete_expanded_duplication_alignment_workload_scenarios',
                  pairs=len(pairs), directed_dispositions=sum(counts.values()), source_hashes=bindings,
                  scope='All expanded pairs, both input orders and masks, assuming fresh alignment of every ready input. Scenarios apply historical mean task wall time within mask and longer-input length bins; uncovered bins are explicit and excluded from estimated hours. No result reuse, native alignment, calibrated ETA, memory bound, or queue/IO overhead estimate. Bin membership alone does not guarantee runtime transferability; sequence divergence and both lengths also matter.')
    verify()
    output = Path(plan['output'])
    output.mkdir(parents=True, exist_ok=False)
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
