#!/usr/bin/env python3
"""Reconstruct each boundary pair with CSV dictionaries, independently of pandas joins."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', required=True, type=Path)
    ap.add_argument('--result', required=True, type=Path)
    args = ap.parse_args()
    r = json.loads((args.result / 'receipt.json').read_text())
    for p, digest in r['source_hashes'].items():
        if sha(p) != digest:
            raise ValueError('Source hash differs')
    data = args.result / 'boundary_sensitivity.tsv.gz'
    if sha(data) != r['artifacts'][data.name]:
        raise ValueError('Output hash differs')
    intervals = {}
    with (args.manifest / 'intervals.tsv').open() as f:
        for x in csv.DictReader(f, delimiter='\t'):
            key = x['interval_id']
            if key in intervals:
                raise ValueError('Duplicate interval')
            intervals[key] = (x['model_key'], int(x['start']), int(x['end']))
    pairs = {}
    with gzip.open(args.manifest / 'boundary_links.tsv.gz', 'rt') as f:
        for x in csv.DictReader(f, delimiter='\t'):
            model, start, end = intervals[x['interval_id']]
            if model != x['model_key']:
                raise ValueError('Wrong model')
            bounds = pairs.setdefault((model, x['hit_id']), {})
            if x['boundary'] in bounds:
                raise ValueError('Duplicate boundary')
            bounds[x['boundary']] = (start, end)
    seen = set(); models = set(); identical = large = fraction = 0
    with gzip.open(data, 'rt') as f:
        for x in csv.DictReader(f, delimiter='\t'):
            key = x['model_key'], x['hit_id']
            if key in seen or key not in pairs:
                raise ValueError('Extra/repeated pair')
            seen.add(key); models.add(key[0]); bounds = pairs[key]
            if set(bounds) != {'alignment', 'envelope'}:
                raise ValueError('Missing boundary')
            a, b = bounds['alignment']; c, d = bounds['envelope']
            expected = {'alignment_start': a, 'alignment_end': b, 'alignment_residues': b-a+1,
                        'envelope_start': c, 'envelope_end': d, 'envelope_residues': d-c+1,
                        'n_terminal_extension': a-c, 'c_terminal_extension': d-b,
                        'added_residues': (d-c)-(b-a)}
            if not c <= a <= b <= d or any(int(x[k]) != v for k, v in expected.items()):
                raise ValueError('Reconstructed interval differs')
            extra = expected['added_residues']; length = b-a+1
            if not math.isclose(float(x['added_fraction_of_alignment']), extra/length, abs_tol=1e-14):
                raise ValueError('Fraction differs')
            identical += extra == 0; large += extra >= 10; fraction += 5*extra >= length
    checks = {'candidate_model_hit_pairs': len(seen), 'source_models': len(models),
              'identical_boundary_pairs': identical, 'different_boundary_pairs': len(seen)-identical,
              'pairs_with_at_least_10_added_residues': large,
              'pairs_with_at_least_20_percent_added_residues': fraction}
    if seen != set(pairs) or any(r[k] != v for k, v in checks.items()):
        raise ValueError('Incomplete scope or summary mismatch')
    receipt = {'status': 'passed_all_boundary_pair_rows_and_threshold_counts', **checks,
               'producer_receipt_sha256': sha(args.result / 'receipt.json'), 'script_sha256': sha(__file__),
               'scope': 'All emitted identities, endpoints, lengths, extensions, fractions and threshold counts independently reconstructed with CSV dictionaries from original manifest. Quantiles not independently recomputed. Biological boundaries and structural effects not validated.'}
    (args.result / 'readback.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
