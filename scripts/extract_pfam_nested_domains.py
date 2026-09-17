#!/usr/bin/env python3
"""Preserve every repeated Pfam NE tag for subsequent overlap review."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--metadata', type=Path, required=True)
    ap.add_argument('--release', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    release = json.loads(a.release.read_text())
    source = next(r for r in release['files'] if r['uncompressed_file'] == a.metadata.name)
    assert sha(a.metadata) == source['uncompressed_sha256']
    records = []
    row = defaultdict(list)
    raw_ne_lines = 0
    with a.metadata.open() as handle:
        for line in handle:
            if line.startswith('#=GF '):
                _, key, value = line.rstrip().split(maxsplit=2)
                row[key].append(value.strip())
                raw_ne_lines += key == 'NE'
            elif line.strip() == '//':
                assert len(row['ID']) == len(row['AC']) == 1
                records.append(dict(row))
                row = defaultdict(list)
    assert not row and len(records) == release['families']
    names = {r['ID'][0]: r for r in records}
    assert len(names) == len(records)
    output = []
    for record in records:
        for occurrence, nested in enumerate(record.get('NE', []), 1):
            assert nested in names, nested
            target = names[nested]
            output.append(dict(containing_accession=record['AC'][0], containing_name=record['ID'][0],
                               nested_accession=target['AC'][0], nested_name=nested,
                               containing_clan=';'.join(record.get('CL', [])),
                               nested_clan=';'.join(target.get('CL', [])), source_NE_occurrence=occurrence))
    assert len(output) == raw_ne_lines
    pairs = {(r['containing_accession'], r['nested_accession']) for r in output}
    a.output.mkdir(parents=True)
    table = a.output / 'nested_domains.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(output)
    result = dict(status='complete_release_nested_domain_metadata', pfam_families=len(records),
                  nested_relationships=len(pairs), source_NE_occurrences=len(output),
                  repeated_pair_occurrences=len(output)-len(pairs),
                  containing_families=sum(bool(r.get('NE')) for r in records),
                  families_with_multiple_NE_tags=sum(len(r.get('NE', [])) > 1 for r in records),
                  maximum_NE_tags=max(len(r.get('NE', [])) for r in records),
                  source_sha256=sha(a.metadata), release_receipt_sha256=sha(a.release),
                  script_sha256=sha(Path(__file__)), artifacts={table.name: sha(table)},
                  scope='All source NE occurrences preserved and names resolved within the same release. These are curated model relationships, not observed nesting in project proteins or automatic permission to accept all overlaps.')
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
