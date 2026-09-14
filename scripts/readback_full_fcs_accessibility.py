#!/usr/bin/env python3
"""Independently compare every merged exposure field with filtered baseline rows."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--merged', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    receipt = json.loads((a.merged / 'receipt.json').read_text())
    for key, path in receipt['source_paths'].items():
        folder = Path(path)
        if digest(folder / 'receipt.json') != receipt['source_receipts'][key]:
            raise ValueError('Source receipt changed')
        source = json.loads((folder / 'receipt.json').read_text())
        for filename, expected in source['artifacts'].items():
            if digest(folder / filename) != expected:
                raise ValueError('Source artifact changed')
    for filename, expected in receipt['artifacts'].items():
        if digest(a.merged / filename) != expected:
            raise ValueError('Merged artifact changed')
    with (Path(receipt['source_paths']['inputs']) / 'omitted_marker_observations.tsv').open() as handle:
        excluded = {(r['marker'], r['taxon_id']) for r in csv.DictReader(handle, delimiter='\t')}
    original = Path(receipt['source_paths']['baseline']) / 'normalized_paired_sites.tsv.gz'
    count = 0
    with gzip.open(original, 'rt') as src, gzip.open(a.merged / 'normalized_paired_sites.tsv.gz', 'rt') as dst:
        before, after = csv.reader(src, delimiter='\t'), csv.reader(dst, delimiter='\t')
        fields = next(before)
        if fields != next(after):
            raise ValueError('Header differs')
        m, t = fields.index('marker'), fields.index('taxon_id')
        retained = (row for row in before if (row[m], row[t]) not in excluded)
        for left, right in itertools.zip_longest(retained, after):
            if left != right:
                raise ValueError('Filtered baseline and merged rows differ')
            count += 1
    if count != receipt['rows']:
        raise ValueError('Merged count differs')
    result = {'status': 'passed_full_filtered_baseline_field_readback', 'rows': count,
              'fields_per_row': len(fields), 'field_values_checked': count * len(fields),
              'merged_receipt_sha256': digest(a.merged / 'receipt.json'),
              'script_sha256': digest(Path(__file__)),
              'scope': 'All saved source/merged hashes checked. Every retained output field and row order matches baseline after the specified whole marker/taxon omissions. Does not independently recalculate solvent area or validate sequence provenance.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
