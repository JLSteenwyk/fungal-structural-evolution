#!/usr/bin/env python3
"""Independently check all figure rows against the complete candidate table."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    receipt = json.loads(Path(plan['receipt']).read_text())
    assert receipt['plan_sha256'] == sha(args.plan)
    for path, digest in {**plan['pins'], **receipt['artifacts']}.items():
        assert sha(path) == digest, path
    source_path = Path(plan['source']) / 'candidate_control_summary.tsv'
    source_receipt = json.loads((Path(plan['source']) / 'receipt.json').read_text())
    assert source_receipt['artifacts'][source_path.name] == sha(source_path)
    source = rows(source_path)
    plotted = rows(Path(plan['output_prefix']).with_suffix('.tsv'))
    key_fields = ['family', 'gene_a', 'gene_b', 'pfam_accession']
    def key(row):
        return tuple(row[field] for field in key_fields)
    lookup = {key(row): row for row in source}
    assert len(lookup) == len(source) == len(plotted) == 48
    assert {key(row) for row in plotted} == set(lookup)
    endpoints = {
        'original_min': 'all_alternatives_min_rmsd_ar_minus_br',
        'original_max': 'all_alternatives_max_rmsd_ar_minus_br',
        'control_min': 'control_contrast_min',
        'control_max': 'control_contrast_max',
    }
    classifications = Counter()
    for number, row in enumerate(plotted, 1):
        original = lookup[key(row)]
        assert int(row['display_row']) == number
        for field in row.keys() - endpoints.keys() - {'display_row'}:
            assert row[field] == original[field], (number, field)
        for target, field in endpoints.items():
            assert float(row[target]) == float(original[field]), (number, target)
        low, high = float(row['control_min']), float(row['control_max'])
        assert low <= high
        if low > .1:
            category = 'a_farther_from_reference'
        elif high < -.1:
            category = 'b_farther_from_reference'
        elif low < -.1 and high > .1:
            category = 'direction_flip_beyond_margin'
        else:
            category = 'touches_or_enters_margin_band'
        assert row['control_direction_margin_0_1'] == category
        classifications[category] += 1
    assert sum(int(row['control_retains_original_direction']) for row in plotted) == 44
    assert classifications['direction_flip_beyond_margin'] == 3
    assert classifications['touches_or_enters_margin_band'] == 1
    result = dict(status='passed_full_figure_table_readback', candidate_rows=48,
                  endpoints_checked=192, classifications=dict(classifications),
                  plan_sha256=sha(args.plan), producer_receipt_sha256=sha(plan['receipt']),
                  source_table_sha256=sha(source_path), checker_sha256=sha(__file__),
                  scope='All candidate identities, copied fields, endpoints, classification thresholds and artifact hashes. Rendered layout requires visual inspection.')
    assert not args.output.exists()
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
