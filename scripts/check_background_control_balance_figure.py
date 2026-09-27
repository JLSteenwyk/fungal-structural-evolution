#!/usr/bin/env python3
"""Verify all figure points against the independently audited balance table."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.receipt.read_text())
    for path, expected in (receipt['source_bindings'] | receipt['artifacts']).items():
        assert digest(path) == expected, path
    source = next(p for p in receipt['source_bindings'] if p.endswith('/covariate_balance.tsv'))
    with open(source) as handle:
        rows = {(r['guide'], r['feature']): r for r in csv.DictReader(handle, delimiter='\t')
                if r['scenario_id'] == 'S45' and r['policy'] == 'alignment_evalue'}
    assert len(rows) == 16
    fields = {'matched_balance': 'standardized_mean_difference',
              'target_selection_shift': 'selection_shift_in_baseline_sd'}
    table = next(p for p in receipt['artifacts'] if p.endswith('.tsv'))
    seen = set()
    with open(table) as handle:
        for point in csv.DictReader(handle, delimiter='\t'):
            key = point['guide'], point['feature']
            identity = (*key, point['panel'])
            assert identity not in seen
            seen.add(identity)
            row = rows[key]
            assert point['pairs'] == row['pairs']
            assert point['baseline_targets'] == row['baseline_targets']
            assert math.isfinite(float(point['value']))
            assert float(point['value']) == float(row[fields[point['panel']]])
            if point['panel'] == 'matched_balance':
                pooled = math.sqrt((float(row['target_sd'])**2 + float(row['control_sd'])**2) / 2)
                expected = (float(row['target_mean']) - float(row['control_mean'])) / pooled
                assert math.isclose(float(point['value']), expected, rel_tol=1e-9, abs_tol=1e-12)
    assert seen == {(*key, panel) for key in rows for panel in fields}
    args.output.write_text(json.dumps({
        'status': 'passed_all_background_control_balance_figure_points',
        'figure_receipt_sha256': digest(args.receipt), 'points': len(seen),
        'scope': 'All 32 values and sample counts checked against the fully audited source table; matched pooled-SD differences also recalculated from moments. Selection shifts inherit the independently verified original-target SD calculation. Export hashes verified. No structural effect inference.'
    }, indent=2) + '\n')
    print(args.output.read_text())


if __name__ == '__main__':
    main()
