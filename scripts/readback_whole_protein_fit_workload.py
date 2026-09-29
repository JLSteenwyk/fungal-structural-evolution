#!/usr/bin/env python3
"""Check workload coverage and arithmetic using grouped, compensated sums."""
import json
import math
import subprocess
import time
from collections import Counter
from pathlib import Path

import psutil
from screen_duplication_alignment_reuse import sha


def reconstruct(recipes, references):
    trees = sorted({r['tree'] for r in references})
    refs = {(r['records_upper_bound'], r['coefficients'], r['tree']): r for r in references}
    assert len(refs) == len(references)
    counts, workload, sizes, seen = Counter(), Counter(), Counter(), set()
    maxima = [0, 0]
    for row in recipes:
        assert row['fit_input_id'] not in seen
        seen.add(row['fit_input_id'])
        spec = row['specification']
        n, p = spec['records'], len(spec['columns']) - 1
        assert n == row['records'] and n > p > 0
        maxima = [max(maxima[0], n), max(maxima[1], p)]
        bucket = min(b for b in (250, 500, 1000, 2000, 5000, 10000, 20000, 1000000000) if b >= n)
        for tree in trees:
            key = bucket, p, tree
            ref = refs.get(key)
            if ref is None:
                status = 'missing_size_stratum'
            elif n < ref['records_min'] or n > ref['records_max']:
                status = 'outside_observed_record_range'
            else:
                status = 'covered_size_stratum'
                sizes[key] += 1
            counts[status] += 1
            workload[(row['variant'], row['outcome'], tree, status)] += 1
    totals = {
        'mean_seconds': math.fsum(v * refs[k]['wall_seconds_mean'] for k, v in sizes.items()),
        'p95_seconds': math.fsum(v * refs[k]['wall_seconds_quantiles']['p95'] for k, v in sizes.items()),
        'serialized_bytes': math.fsum(v * refs[k]['fit_bytes_mean'] for k, v in sizes.items()),
    }
    return seen, counts, workload, totals, maxima


def main():
    pp = Path('metadata/whole_protein_fit_workload_readback_plan_20260929.json')
    plan = json.loads(pp.read_text())
    bindings = {str(pp): sha(pp), **plan['pins']}
    launch = json.loads(Path(plan['launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p = psutil.Process(launch['pid'])
            if abs(p.create_time() - launch['created']) > .01 or p.status() == psutil.STATUS_ZOMBIE:
                break
            assert p.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(line.split('=', 1) for line in subprocess.check_output(
        ['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    producer = json.loads(Path(launch['plan']).read_text())
    result_path = Path(producer['output']) / 'workload.json'
    result = json.loads(result_path.read_text())
    assert result['status'] == 'complete_whole_protein_workload_planning_reference'
    bindings.update(result['source_hashes'])
    bindings[str(result_path)] = sha(result_path)
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    root = Path(producer['inventory'])
    receipt = json.loads((root / 'receipt.json').read_text())
    refs = json.loads(Path(producer['cost_strata']).read_text())
    with (root / 'unique_input_recipes.jsonl').open() as handle:
        seen, counts, workload, totals, maxima = reconstruct((json.loads(line) for line in handle), refs)
    assert len(seen) == receipt['unique_inputs'] == result['unique_inputs']
    assert sum(counts.values()) == receipt['unique_tree_fits'] == result['tree_fits']
    assert dict(counts) == result['coverage_counts']
    expected = [dict(variant=k[0], outcome=k[1], tree=k[2], reference_status=k[3], fits=v) for k, v in sorted(workload.items())]
    assert expected == result['workload']
    assert maxima == [result['max_records'], result['max_coefficients']]
    assert result['full_workload_cost_covered'] == (counts['covered_size_stratum'] == sum(counts.values()))
    def close(a, b):
        assert math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-8), (a, b)
    for name, value in totals.items():
        close(value, result['covered_fit_totals'][name])
    scenarios = result['covered_only_scenarios']
    assert [(s['workers'], s['cost_multiplier']) for s in scenarios] == [(w, m) for w in (8, 16, 32) for m in (1, 4)]
    for s in scenarios:
        close(s['covered_mean_hours'], totals['mean_seconds'] / 3600 / s['workers'] * s['cost_multiplier'])
        close(s['covered_p95_cost_hours'], totals['p95_seconds'] / 3600 / s['workers'] * s['cost_multiplier'])
    verify()
    proof = dict(status='passed_whole_protein_workload_arithmetic_readback', unique_inputs=len(seen), tree_fits=sum(counts.values()), coverage_counts=dict(counts), producer_terminal_state=state, source_result_sha256=sha(result_path), checker_sha256=sha(__file__), scope='Grouped compensated sums independently check coverage and arithmetic. Empirical timing transfer and ideal parallel scaling remain uncalibrated; no full ETA for uncovered inputs, memory bound or fit launch.')
    with Path(plan['proof']).open('x') as handle:
        json.dump(proof, handle, indent=2)
        handle.write('\n')
    print(json.dumps(proof), flush=True)


if __name__ == '__main__':
    main()
