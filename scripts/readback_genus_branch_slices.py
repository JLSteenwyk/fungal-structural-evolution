#!/usr/bin/env python3
"""Read back every conditional branch-slice point from immutable HyPhy logs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def equal(a, b):
    if not math.isfinite(float(a)) or not math.isclose(float(a), float(b), rel_tol=1e-11, abs_tol=1e-9):
        raise ValueError(f'Numerical mismatch: {a}, {b}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['slices', 'fits', 'normalized', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    a = parser.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    receipt = json.loads((a.slices / 'receipt.json').read_text())
    for path, key in [(a.fits / 'receipt.json', 'source_fit_receipt_sha256'),
                      (a.normalized / 'receipt.json', 'normalization_receipt_sha256')]:
        if sha(path) != receipt[key]:
            raise ValueError('Source receipt changed')
    for name, digest in receipt['artifacts'].items():
        if sha(a.slices / name) != digest:
            raise ValueError('Slice table changed')
    normalized_receipt = json.loads((a.normalized / 'receipt.json').read_text())
    if sha(a.normalized / 'normalized_branches.tsv') != normalized_receipt['artifacts']['normalized_branches.tsv']:
        raise ValueError('Normalization table changed')
    branches = {}
    for row in rows(a.normalized / 'normalized_branches.tsv'):
        branches.setdefault(row['case_id'], []).append(row)
    summaries = rows(a.slices / 'case_review.tsv')
    points = rows(a.slices / 'likelihood_slices.tsv')
    grouped = {}
    for point in points:
        grouped.setdefault(point['case_id'], []).append(point)
    checks = {r['case_id']: r for r in receipt['checks']}
    if len(summaries) != receipt['cases'] or len(checks) != len(summaries) or len(points) != receipt['slice_evaluations']:
        raise ValueError('Wrong row counts')
    if set(branches) != set(grouped) or set(checks) != set(branches) or {r['case_id'] for r in summaries} != set(checks):
        raise ValueError('Case grid mismatch')
    largest_improvement = -math.inf
    for summary in summaries:
        case = summary['case_id']; check = checks[case]
        for path, key in [(a.fits / case / 'receipt.json', 'source_receipt_sha256'),
                          (a.slices / case / 'slice.bf', 'slice_script_sha256'),
                          (a.slices / case / 'slice.log', 'slice_log_sha256')]:
            if sha(path) != check[key]:
                raise ValueError('Raw artifact changed')
        target = min(branches[case], key=lambda r: (-float(r['dS_upstream_equal_alternative_convention']), r['node']))
        if summary['target_node'] != target['node']:
            raise ValueError('Wrong target branch')
        ds = float(target['dS_upstream_equal_alternative_convention'])
        equal(summary['target_dS'], ds)
        equal(summary['target_dN'], target['dN_upstream_equal_alternative_convention'])
        baseline = float(json.loads((a.fits / case / 'receipt.json').read_text())['log_likelihood'])
        log = (a.slices / case / 'slice.log').read_text()
        raw = re.findall(r'^SLICE\t(\S+)\t(\S+)\t(\S+)', log, re.M)
        expected_factors = [.1, .25, .5, 1., 2., 4., 10.]
        if len(raw) != 7 or [float(r[0]) for r in raw] != expected_factors or len(grouped[case]) != 7:
            raise ValueError('Wrong raw slice grid')
        saved = {float(r['branch_multiplier']): r for r in grouped[case]}
        if set(saved) != set(expected_factors):
            raise ValueError('Wrong table slice grid')
        deltas = {}
        for factor, parameter, ll in raw:
            factor, parameter, ll = float(factor), float(parameter), float(ll)
            row = saved[factor]
            if row['node'] != target['node']:
                raise ValueError('Wrong point branch')
            equal(row['branch_parameter'], parameter)
            equal(parameter, float(summary['target_branch_parameter']) * factor)
            equal(row['normalized_dS'], ds * factor)
            equal(row['log_likelihood'], ll)
            equal(row['delta_log_likelihood_from_saved_fit'], ll - baseline)
            deltas[factor] = ll - baseline
        restored = re.findall(r'^RESTORED=(\S+)', log, re.M)
        if len(restored) != 1:
            raise ValueError('Missing restored baseline')
        equal(restored[0], baseline); equal(deltas[1.], 0)
        for factor, field in [(.5, 'delta_ll_at_half'), (2., 'delta_ll_at_double'), (10., 'delta_ll_at_tenfold')]:
            equal(summary[field], deltas[factor])
        equal(summary['maximum_slice_improvement'], max(deltas.values()))
        largest_improvement = max(largest_improvement, max(deltas.values()))
    result = {'status': 'passed_full_raw_log_and_table_readback', 'cases': len(summaries),
              'points': len(points), 'maximum_conditional_improvement': largest_improvement,
              'source_receipt_sha256': sha(a.slices / 'receipt.json'),
              'script_sha256': sha(Path(__file__)),
              'scope': 'Every raw log, point, target selection, distance scaling and summary checked. No new likelihood evaluation or optimization; does not establish model adequacy.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
