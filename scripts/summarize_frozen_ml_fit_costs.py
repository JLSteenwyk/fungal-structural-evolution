#!/usr/bin/env python3
"""Summarize full frozen fit-prefix timings for subsequent resource planning."""
import hashlib
import json
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from screen_duplication_alignment_reuse import sha


def main():
    plan_path = Path('metadata/frozen_ml_fit_costs_plan_20260929.json')
    plan = json.loads(plan_path.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    root = Path(plan['census'])
    receipt = json.loads((root / 'receipt.json').read_text())
    prefix = root / 'manifest_prefix.jsonl'
    assert sha(prefix) == receipt['artifacts']['manifest_prefix.jsonl']
    groups, statuses, seen = defaultdict(list), Counter(), set()
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    fields = ['fit_input_id','tree','polynomial_degree','records','coefficients','wall_seconds','fit_bytes','status']
    import csv
    with (out / 'fit_costs.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for line in prefix.open():
            row = json.loads(line)
            key = row['fit_input_id'], row['tree']
            assert key not in seen
            seen.add(key)
            path = Path(row['path'])
            data = path.read_bytes()
            assert hashlib.sha256(data).hexdigest() == row['sha256'], str(path)
            saved = json.loads(data)
            p = saved['payload']
            assert saved['fit_input_id'] == key[0] and saved['tree'] == key[1]
            assert p['status'] == row['status']
            seconds = p['wall_seconds']
            assert np.isfinite(seconds) and seconds >= 0
            n, columns = int(p['records']), len(p['beta'])
            assert n > columns > 0
            record = dict(fit_input_id=key[0], tree=key[1], polynomial_degree=saved['polynomial_degree'], records=n, coefficients=columns, wall_seconds=seconds, fit_bytes=len(data), status=p['status'])
            writer.writerow(record)
            bucket = next(b for b in [250,500,1000,2000,5000,10000,20000,1000000000] if n <= b)
            groups[(bucket, columns, key[1])].append((seconds,len(data),n))
            statuses[p['status']] += 1
            if len(seen) % 10000 == 0:
                print('Summarized frozen fit costs', len(seen), flush=True)
    assert len(seen) == receipt['prefix_rows'] and dict(statuses) == receipt['manifest_status_counts']
    summary = []
    for key, values in sorted(groups.items()):
        a = np.asarray(values)
        summary.append(dict(records_upper_bound=key[0], coefficients=key[1], tree=key[2], fits=len(a), records_min=int(a[:,2].min()), records_max=int(a[:,2].max()), wall_seconds_total=float(a[:,0].sum()), wall_seconds_mean=float(a[:,0].mean()), wall_seconds_quantiles=dict(zip(['p50','p90','p95','p99','maximum'], map(float,np.quantile(a[:,0],[.5,.9,.95,.99,1])))), fit_bytes_mean=float(a[:,1].mean()), fit_bytes_total=int(a[:,1].sum())))
    (out / 'cost_strata.json').write_text(json.dumps(summary,indent=2)+'\n')
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    result = dict(status='complete_frozen_ml_fit_cost_census', fits=len(seen), strata=len(summary), status_counts=dict(statuses), source_hashes=plan['pins'], plan_sha256=sha(plan_path), script_sha256=sha(__file__), artifacts={p.name:sha(p) for p in out.iterdir()}, scope='All completed fits in the previously frozen prefix, including review flags. Wall times include concurrent machine contention and are not CPU times. Completion-order sample, not a random benchmark. Whole-protein response/design and covariance grouping differences can change costs; these strata alone are not an ETA, memory bound, calibrated prediction or authorization to launch fits.')
    (out / 'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
