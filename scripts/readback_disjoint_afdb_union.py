#!/usr/bin/env python3
"""Read back every union row against its unchanged source cohorts."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); plan = json.loads(a.plan.read_text()); root = Path(plan['output'])
    pins = {str(a.plan): sha(a.plan), **plan['pins']}
    for name in ['mapping', 'encodings']:
        directory = root / name; receipt = json.loads((directory / 'receipt.json').read_text())
        if receipt['plan_sha256'] != sha(a.plan):
            raise ValueError('Different union plan')
        pins[str(directory / 'receipt.json')] = sha(directory / 'receipt.json')
        pins.update({str(directory / f): h for f, h in receipt['artifacts'].items()})
    def verify():
        for path, digest in pins.items():
            if sha(Path(path)) != digest:
                raise ValueError('Changed input: ' + path)
    verify()
    for section, filename, key in [('mapping', 'marker_structure_links.tsv', lambda r: (r['marker'], r['taxon_id'])),
                                    ('encodings', 'model_summary.tsv', lambda r: r['model_name'])]:
        expected = {}
        for cohort in plan['cohorts']:
            source = Path(cohort['mapping' if section == 'mapping' else 'encodings'])
            for row in rows(source / filename):
                if key(row) in expected:
                    raise ValueError('Repeated source row')
                expected[key(row)] = row
        actual = rows(root / section / filename)
        if len(actual) != len(expected) or {key(r): r for r in actual} != expected:
            raise ValueError('Union table differs')
    expected_models = []
    for cohort in plan['cohorts']:
        expected_models.extend(json.loads((Path(cohort['mapping']) / 'model_provenance.json').read_text()))
    model_key = lambda m: (m['model_id'], m['version'])
    actual_models = json.loads((root / 'mapping/model_provenance.json').read_text())
    if sorted(actual_models, key=model_key) != sorted(expected_models, key=model_key):
        raise ValueError('Union model provenance differs')
    count = 0
    with gzip.open(root / 'mapping/matrix_to_structure_residues.tsv.gz', 'rb') as target:
        header = next(target)
        for cohort in plan['cohorts']:
            with gzip.open(Path(cohort['mapping']) / 'matrix_to_structure_residues.tsv.gz', 'rb') as source:
                if next(source) != header:
                    raise ValueError('Residue table header differs')
                for line in source:
                    if next(target, None) != line:
                        raise ValueError('Union residue row changed, omitted, or inserted')
                    count += 1
        if next(target, None) is not None or count != plan['residue_links']:
            raise ValueError('Union residue extent differs')
    verify()
    result = dict(status='passed_full_disjoint_afdb_union_row_readback', models=len(actual_models),
                  marker_links=plan['marker_links'], residue_links=count, source_sha256=pins,
                  script_sha256=sha(Path(__file__)),
                  scope='All model provenance, marker-link and encoding-summary rows equal the source union; every decompressed residue row is byte-identical to its source. Does not rerun source coordinate or confidence audits.')
    with a.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}))


if __name__ == '__main__':
    main()
