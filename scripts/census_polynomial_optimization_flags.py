#!/usr/bin/env python3
"""Classify flagged completed fits in an immutable prefix of an active manifest."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan_hash = sha(args.plan)
    plan = json.loads(args.plan.read_text())
    source = Path(plan['output']) / 'fit_manifest.jsonl'
    raw = source.read_bytes()
    prefix = raw[:raw.rfind(b'\n') + 1]
    assert prefix
    args.output.mkdir(exist_ok=False)
    snapshot = args.output / 'manifest_prefix.jsonl'
    snapshot.write_bytes(prefix)
    counts, reasons, combinations, degrees = Counter(), Counter(), Counter(), Counter()
    rows, seen, bindings = [], set(), {}
    for line in prefix.splitlines():
        record = json.loads(line)
        key = (record['fit_input_id'], record['tree'])
        assert key not in seen
        seen.add(key)
        counts[record['status']] += 1
        if record['status'] == 'ml_candidate_passed_numerical_optimization_checks':
            continue
        path = Path(record['path'])
        assert sha(path) == record['sha256']
        bindings[str(path)] = record['sha256']
        saved = json.loads(path.read_text())
        assert saved['plan_sha256'] == plan_hash and (saved['fit_input_id'], saved['tree']) == key
        payload = saved['payload']
        encoded = json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False)
        assert hashlib.sha256(encoded.encode()).hexdigest() == saved['payload_sha256']
        assert payload['status'] == record['status']
        checks = payload.get('checks', {})
        failed = sorted(k for k, value in checks.items() if (value if k == 'upper_bound_contact' else not value))
        if record['status'] == 'fit_error_requires_review':
            failed.append('fit_error:' + payload['error_type'])
        assert failed, record
        reasons.update(failed)
        combinations['|'.join(failed)] += 1
        degrees[str(saved['polynomial_degree'])] += 1
        candidates = payload.get('candidates', [])
        best = min((c['objective'] for c in candidates), default=None)
        successful = [c['objective'] for c in candidates if c['success']]
        gap = min(successful) - best if successful else None
        rows.append(dict(fit_input_id=key[0], tree=key[1], polynomial_degree=saved['polynomial_degree'],
                         status=record['status'], reasons='|'.join(failed),
                         maximum_absolute_projected_gradient=max(map(abs, payload.get('projected_gradient', [])), default=None),
                         successful_candidate_objective_gap=gap,
                         candidates=len(candidates), successful_candidates=len(successful),
                         source_path=str(path), source_sha256=record['sha256']))
    table = args.output / 'flagged_fits.tsv'
    assert rows
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    with table.open() as handle:
        serialized = list(csv.DictReader(handle, delimiter='\t'))
    assert len(serialized) == len(rows)
    for actual, expected in zip(serialized, rows):
        assert actual == {k: '' if v is None else str(v) for k, v in expected.items()}
    with source.open('rb') as handle:
        assert handle.read(len(prefix)) == prefix
    assert sha(args.plan) == plan_hash
    for path, digest in bindings.items():
        assert sha(path) == digest
    result = dict(status='complete_fixed_prefix_flag_census_not_full_grid_audit',
                  script_sha256=sha(__file__), source_plan=str(args.plan), source_plan_sha256=plan_hash,
                  manifest_path=str(source), prefix_bytes=len(prefix), prefix_rows=len(seen),
                  manifest_status_counts=dict(counts), flagged_fits=len(rows),
                  flag_reason_counts=dict(reasons), flag_combinations=dict(combinations),
                  flagged_polynomial_degrees=dict(degrees),
                  successful_candidate_gap_at_most_1e_minus_7=sum(r['successful_candidate_objective_gap'] is not None and r['successful_candidate_objective_gap'] <= 1e-7 for r in rows),
                  source_hashes=bindings, artifacts={p.name: sha(p) for p in [snapshot, table]},
                  scope='Every non-passing fit in the captured completed-manifest prefix is hash- and payload-verified and classified. Passing fit counts are manifest counts only. Near-equal successful candidates do not resolve flags; no fit statuses, tolerances or source results changed. Full-grid audit and optimization review remain required.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
