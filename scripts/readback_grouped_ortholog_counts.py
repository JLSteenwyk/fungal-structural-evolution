#!/usr/bin/env python3
"""Check completed grouped-ortholog audit totals and snapshot bindings.

Does not repeat the raw-row identity scan or validate orthology semantics.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def counts(path, key, fields):
    found = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(reader.fieldnames == fields, 'Unexpected count-table schema')
        for row in reader:
            require(None not in row and all(v is not None for v in row.values()),
                    'Malformed count-table row')
            require(row[key] and row[key] not in found, 'Empty or duplicated count key')
            for name in ['rows', 'directed_incidences']:
                require(row[name].isdigit(), 'Nonnegative integer count required')
                row[name] = int(row[name])
            require(row['directed_incidences'] >= row['rows'], 'Invalid incidence count')
            found[row[key]] = row
    return found


def check(plan_path):
    plan = json.loads(plan_path.read_text())
    root = Path(plan['output'])
    receipt_path = root / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    require(receipt['status'] == 'complete_grouped_ortholog_protein_family_identity_audit',
            'Completed identity audit required')
    require(receipt['plan_sha256'] == digest(plan_path), 'Changed execution plan')
    producer = 'scripts/audit_grouped_ortholog_identities.py'
    require(receipt['script_sha256'] == plan['pins'][producer], 'Producer binding differs')
    require(plan['supplement_receipt'] in plan['pins'], 'Unpinned source snapshot')
    pins = dict(plan['pins'])
    pins[str(plan_path)] = digest(plan_path)
    pins[str(receipt_path)] = digest(receipt_path)
    snapshot = json.loads(Path(plan['supplement_receipt']).read_text())
    require(snapshot['status'] == 'complete_separate_small_family_ortholog_supplement',
            'Completed source snapshot required')
    for path, expected in snapshot['input_hashes'].items():
        require(path not in pins or pins[path] == expected, 'Conflicting source pins')
        pins[path] = expected
    for path, expected in pins.items():
        require(digest(path) == expected, 'Changed input: ' + path)
    require(set(receipt['artifacts']) == {'table_counts.tsv', 'family_counts.tsv'},
            'Unexpected audit artifact set')
    for name, expected in receipt['artifacts'].items():
        require(digest(root / name) == expected, 'Changed count artifact: ' + name)
        pins[str(root / name)] = expected
    tables = counts(root / 'table_counts.tsv', 'taxon',
                    ['taxon', 'rows', 'directed_incidences', 'sha256'])
    families = counts(root / 'family_counts.tsv', 'family',
                      ['family', 'rows', 'directed_incidences'])
    native = {Path(path).stem: sha for path, sha in snapshot['native_ortholog_hashes'].items()}
    require(len(native) == len(snapshot['native_ortholog_hashes']), 'Ambiguous native table names')
    require(set(tables) == set(native), 'Native table coverage differs')
    require(all(row['sha256'] == native[taxon] for taxon, row in tables.items()),
            'Recorded native table hash differs from snapshot')
    for name in ['rows', 'directed_incidences']:
        require(sum(row[name] for row in tables.values()) == receipt[name]
                == sum(row[name] for row in families.values()), 'Count totals differ: ' + name)
    require(receipt['rows'] == snapshot['native_rows_scanned'], 'Snapshot row count differs')
    require(receipt['taxa'] == len(tables), 'Taxon count differs')
    require(receipt['families_with_rows'] == len(families), 'Family count differs')
    sequence_files = [Path(path) for path in snapshot['input_hashes']
                      if Path(path).name == 'SequenceIDs.txt']
    require(len(sequence_files) == 1, 'Ambiguous source protein file')
    with sequence_files[0].open() as handle:
        proteins = sum(1 for line in handle if line.strip())
    require(receipt['proteins'] == proteins, 'Protein count differs')
    return dict(status='complete_grouped_ortholog_count_and_binding_readback',
                taxa=len(tables), proteins=proteins, families_with_rows=len(families),
                rows=receipt['rows'], directed_incidences=receipt['directed_incidences'],
                source_sha256=pins, script_sha256=digest(__file__),
                scope='Recomputed totals from both complete count tables, exact native table coverage and recorded hash equality to the pinned snapshot, rehashed plan/producer/source identities/count artifacts. Raw ortholog tables were not rescanned; per-row identities, per-family assignments, reciprocal equality, cross-row duplicates, missing pairs and biological orthology were not independently tested by this readback.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = check(args.plan)
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ['status', 'taxa', 'rows', 'directed_incidences']}))


if __name__ == '__main__':
    main()
