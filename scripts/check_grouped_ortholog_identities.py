#!/usr/bin/env python3
"""Reproduce grouped-ortholog audit fixture and corruption checks."""
import argparse
import copy
import csv
import json
import subprocess
import sys
from pathlib import Path

from assess_small_family_output_exposure import sha
from audit_grouped_ortholog_identities import inspect_row, labels


def write_json(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    snapshot_path = Path('metadata/native_small_family_supplement_sorted_receipt.json')
    snapshot = json.loads(snapshot_path.read_text())
    source = Path('results/orthology/native-small-family-output-fixture-v2/Source/WorkingDirectory')
    for path, digest in {**snapshot['input_hashes'], **snapshot['native_ortholog_hashes']}.items():
        if sha(path) != digest:
            raise ValueError('Fixture changed: ' + path)
    lookup = labels(source)
    path = Path(sorted(snapshot['native_ortholog_hashes'])[0])
    with path.open('rb') as handle:
        next(handle)
        raw = next(handle)
    row = next(csv.reader([raw.decode()], delimiter='\t'))
    taxon = path.stem
    variants = {
        'unknown_family': ['OG9999999', *row[1:]],
        'same_species': [row[0], taxon, *row[2:]],
        'unknown_species': [row[0], '__missing__', *row[2:]],
        'unknown_protein': [row[0], row[1], '__missing__', row[3]],
        'repeated_protein': [row[0], row[1], row[2] + ', ' + row[2], row[3]],
        'malformed': row[:3],
    }
    rejected = []
    for name, variant in variants.items():
        try:
            inspect_row(('\t'.join(variant) + '\n').encode(), taxon, lookup)
        except ValueError:
            rejected.append(name)
        else:
            raise AssertionError('Accepted corrupted row: ' + name)
    producer = 'scripts/audit_grouped_ortholog_identities.py'
    plan = dict(output=str(args.output / 'baseline'), supplement_receipt=str(snapshot_path),
                pins={p: sha(p) for p in [str(snapshot_path), producer,
                                        'scripts/assess_small_family_output_exposure.py']},
                resources=dict(minimum_free_disk_gib=0))

    def run(name, selected, expected_error=None):
        plan_path = args.output / (name + '.plan.json')
        write_json(plan_path, selected)
        result = subprocess.run([sys.executable, producer, '--plan', str(plan_path)],
                                capture_output=True, text=True, timeout=60)
        (args.output / (name + '.stdout.txt')).write_text(result.stdout)
        (args.output / (name + '.stderr.txt')).write_text(result.stderr)
        if expected_error is None:
            if result.returncode:
                raise AssertionError(result.stderr)
        elif result.returncode == 0 or expected_error not in result.stderr:
            raise AssertionError('Incorrect failure for ' + name)
        return result

    run('baseline', plan)
    receipt = json.loads((args.output / 'baseline/receipt.json').read_text())
    if tuple(receipt[k] for k in ['taxa', 'proteins', 'rows', 'directed_incidences']) != (4, 18, 20, 22):
        raise AssertionError('Incorrect native fixture totals')
    for name, digest in receipt['artifacts'].items():
        if sha(args.output / 'baseline' / name) != digest:
            raise AssertionError('Changed output')
    # Rebuild per-table counts by explicit pair enumeration on this tiny fixture.
    expected = {}
    for filename in snapshot['native_ortholog_hashes']:
        table = Path(filename)
        nrows = incidences = 0
        with table.open() as handle:
            reader = csv.reader(handle, delimiter='\t')
            next(reader)
            for family, target, left, right in reader:
                nrows += 1
                incidences += len([(a, b) for a in left.split(', ') for b in right.split(', ')])
        expected[table.stem] = (nrows, incidences)
    with (args.output / 'baseline/table_counts.tsv').open() as handle:
        actual = {r['taxon']: (int(r['rows']), int(r['directed_incidences']))
                  for r in csv.DictReader(handle, delimiter='\t')}
    if actual != expected:
        raise AssertionError('Per-table fixture counts differ')
    changed = copy.deepcopy(plan)
    changed['output'] = str(args.output / 'wrong-pin-output')
    changed['pins'][str(snapshot_path)] = '0' * 64
    run('wrong-pin', changed, 'Changed plan input')
    if Path(changed['output']).exists():
        raise AssertionError('Invalid provenance created output')
    run('refuse-overwrite', plan, 'FileExistsError')
    result = dict(status='passed_reproducible_grouped_identity_checks',
                  rejected_row_cases=rejected, cli_checks=['baseline', 'wrong-pin', 'refuse-overwrite'],
                  expected_native_fixture=dict(taxa=4, proteins=18, rows=20, directed_incidences=22),
                  source_sha256={str(snapshot_path): sha(snapshot_path), producer: sha(producer)},
                  script_sha256=sha(__file__),
                  artifacts={str(p.relative_to(args.output)): sha(p) for p in args.output.rglob('*') if p.is_file()},
                  scope='Software fixture checks only. No full production completion or ortholog-pair semantics claim.')
    write_json(args.output / 'receipt.json', result)
    print(json.dumps({k: v for k, v in result.items() if k != 'artifacts'}, indent=2))


if __name__ == '__main__':
    main()
