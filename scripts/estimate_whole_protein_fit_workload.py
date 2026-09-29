#!/usr/bin/env python3
"""Transfer empirical costs only to covered model-size strata after input audit."""
import json
import time
import subprocess
from pathlib import Path
from collections import Counter
import psutil
from screen_duplication_alignment_reuse import sha

BOUNDS = (250, 500, 1000, 2000, 5000, 10000, 20000, 1000000000)


def reference_status(n, columns, tree, references):
    bucket = next(b for b in BOUNDS if n <= b)
    ref = references.get((bucket, columns, tree))
    if ref is None:
        return 'missing_size_stratum', None
    if not ref['records_min'] <= n <= ref['records_max']:
        return 'outside_observed_record_range', None
    return 'covered_size_stratum', ref


def main():
    pp = Path('metadata/whole_protein_fit_workload_plan_20260929.json')
    plan = json.loads(pp.read_text())
    bindings = {str(pp): sha(pp), **plan['pins']}
    launch = json.loads(Path(plan['audit_launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p = psutil.Process(launch['pid'])
            if abs(p.create_time() - launch['created']) > .01 or p.status() == psutil.STATUS_ZOMBIE:
                break
            assert p.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state == dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    root = Path(plan['inventory'])
    receipt = json.loads((root / 'receipt.json').read_text())
    proof = json.loads(Path(plan['audit_proof']).read_text())
    assert proof['status'] == 'passed_full_whole_protein_input_inventory_readback'
    assert proof['source_receipt_sha256'] == sha(root / 'receipt.json')
    bindings[str(root / 'receipt.json')] = sha(root / 'receipt.json')
    bindings[plan['audit_proof']] = sha(plan['audit_proof'])
    bindings.update({str(root / n): h for n,h in receipt['artifacts'].items()})
    def verify():
        for path,digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    references = {(r['records_upper_bound'],r['coefficients'],r['tree']):r for r in json.loads(Path(plan['cost_strata']).read_text())}
    trees = sorted({k[2] for k in references})
    assert len(trees) == receipt['tree_alternatives'] == 5
    counts, workload, seen = Counter(), Counter(), set()
    sums = dict(mean_seconds=0., p95_seconds=0., serialized_bytes=0.)
    max_records, max_columns = 0, 0
    for line in (root / 'unique_input_recipes.jsonl').open():
        recipe = json.loads(line)
        assert recipe['fit_input_id'] not in seen
        seen.add(recipe['fit_input_id'])
        spec = recipe['specification']
        n, columns = spec['records'], len(spec['columns']) - 1
        assert n == recipe['records'] and n > columns > 0
        max_records, max_columns = max(max_records,n), max(max_columns,columns)
        for tree in trees:
            status, ref = reference_status(n,columns,tree,references)
            counts[status] += 1
            workload[(recipe['variant'],recipe['outcome'],tree,status)] += 1
            if ref:
                sums['mean_seconds'] += ref['wall_seconds_mean']
                sums['p95_seconds'] += ref['wall_seconds_quantiles']['p95']
                sums['serialized_bytes'] += ref['fit_bytes_mean']
    assert len(seen) == receipt['unique_inputs'] and sum(counts.values()) == receipt['unique_tree_fits']
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    scenarios = [dict(workers=w, cost_multiplier=m, covered_mean_hours=sums['mean_seconds']*m/w/3600, covered_p95_cost_hours=sums['p95_seconds']*m/w/3600) for w in (8,16,32) for m in (1,4)]
    result = dict(status='complete_whole_protein_workload_planning_reference', unique_inputs=len(seen), tree_fits=sum(counts.values()), coverage_counts=dict(counts), max_records=max_records,max_coefficients=max_columns,covered_fit_totals=sums,covered_only_scenarios=scenarios,full_workload_cost_covered=(counts['covered_size_stratum']==sum(counts.values())),workload=[dict(variant=k[0],outcome=k[1],tree=k[2],reference_status=k[3],fits=v) for k,v in sorted(workload.items())],source_hashes=bindings,script_sha256=sha(__file__),scope='Planning scenarios only, not ETAs or confidence intervals. Uncovered inputs excluded from cost sums and retained explicitly. Representative variant labels can hide exact cross-variant deduplication. Domain-fit completion-order timings include contention; covariance grouping and response differences uncalibrated. Ideal parallel scaling assumed, no setup/audit/refinement overhead. Serialized fit bytes exclude input caches; memory not estimated. No fitting launched.')
    verify()
    (out / 'workload.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','workload']}),flush=True)


if __name__ == '__main__':
    main()
