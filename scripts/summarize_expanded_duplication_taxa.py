#!/usr/bin/env python3
"""Update all taxon coverage rows while preserving audited cohort membership."""
import argparse
import csv
import gzip
import json
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from readback_expanded_duplication_coverage import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--readback', type=Path, required=True)
    p.add_argument('--membership', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    audit = json.loads(args.readback.read_text())
    assert audit['status'] == 'passed_full_expanded_duplication_event_and_sql_join_readback'
    old_receipt_path = args.membership / 'receipt.json'
    old_receipt = json.loads(old_receipt_path.read_text())
    assert old_receipt['status'] == 'complete_duplication_coverage_reconciliation_membership'
    table = args.membership / 'taxon_coverage_with_membership.tsv'
    sources = dict(audit['sources'])
    sources.update(old_receipt['source_sha256'])
    sources.update({str(args.readback): digest(args.readback), str(old_receipt_path): digest(old_receipt_path), str(table): old_receipt['artifacts'][table.name]})
    for path, expected in sources.items():
        assert digest(path) == expected, path
    events = [Path(s) for s in audit['sources'] if s.endswith('/event_coverage.tsv.gz')]
    assert len(events) == 1
    with table.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        rows = list(reader)
    keys = {(r['guide'], r['taxon']) for r in rows}
    assert len(keys) == len(rows)
    counts = defaultdict(Counter)
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE events (guide TEXT, taxon TEXT, label TEXT, same INTEGER, old_label TEXT)')
    batch = []
    n = 0
    with gzip.open(events[0], 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            key = (r['guide'], r['taxon_id'])
            assert key in keys
            c = counts[key]
            c['terminal_singleton_events'] += 1
            c[r['new_coverage_class']] += 1
            c['both_same_model'] += int(r['new_same_model'])
            batch.append((*key, r['new_coverage_class'], int(r['new_same_model']), r['coverage_class']))
            n += 1
            if len(batch) == 10000:
                db.executemany('INSERT INTO events VALUES (?,?,?,?,?)', batch)
                batch.clear()
        db.executemany('INSERT INTO events VALUES (?,?,?,?,?)', batch)
    assert n == audit['event_rows']
    numeric = ['terminal_singleton_events', 'neither_model', 'one_model', 'both_models', 'both_same_model']
    sql = {tuple(r[:2]): dict(zip(numeric, r[2:])) for r in db.execute("SELECT guide,taxon,COUNT(*),SUM(label='neither_model'),SUM(label='one_model'),SUM(label='both_models'),SUM(same) FROM events GROUP BY guide,taxon")}
    old_sql = {tuple(r[:2]): dict(zip(numeric[:4], r[2:])) for r in db.execute("SELECT guide,taxon,COUNT(*),SUM(old_label='neither_model'),SUM(old_label='one_model'),SUM(old_label='both_models') FROM events GROUP BY guide,taxon")}
    changes = []
    for r in rows:
        key = (r['guide'], r['taxon'])
        observed = {k: counts[key][k] for k in numeric}
        assert observed == sql.get(key, dict.fromkeys(numeric, 0))
        assert {k: int(r[k]) for k in numeric[:4]} == old_sql.get(key, dict.fromkeys(numeric[:4], 0))
        assert bool(observed['terminal_singleton_events']) == bool(int(r['in_reconciliation']))
        changes.append(dict(guide=r['guide'], taxon=r['taxon'], species_name=r['species_name'], study_role=r['study_role'], lineage=r['lineage'], in_reconciliation=r['in_reconciliation'], old_both_models=r['both_models'], new_both_models=observed['both_models'], net_both_model_change=observed['both_models']-int(r['both_models'])))
        r.update(observed)
        r['both_model_fraction'] = r['both_models']/r['terminal_singleton_events'] if r['terminal_singleton_events'] else ''
    args.output.mkdir(parents=True, exist_ok=False)
    for name, columns, data in [('taxon_coverage_with_membership.tsv', fields, rows), ('taxon_coverage_change.tsv', list(changes[0]), changes)]:
        path = args.output / name
        with path.open('w') as f:
            w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', lineterminator='\n');w.writeheader();w.writerows(data)
        with path.open() as f:
            actual = list(csv.DictReader(f, delimiter='\t'))
        assert actual == [{k: str(r[k]) for k in columns} for r in data]
    for path, expected in sources.items():
        assert digest(path) == expected, path
    receipt = dict(status='complete_duplication_coverage_reconciliation_membership',
                   source_sha256=sources, script_sha256=digest(__file__), event_rows=n, rows=len(rows),
                   artifacts={p.name: digest(p) for p in args.output.iterdir()},
                   scope='Expanded catalog coverage with original cohort, metadata and denominators preserved. All group counts checked by SQL against streaming counters; serialized rows read back exactly. Uses audited event/model joins; availability only, no confidence qualification or evolutionary test.')
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(dict(rows=len(rows), events=n, newly_covered_taxa_by_guide={g: sum(r['guide']==g and int(r['old_both_models'])==0 and r['new_both_models']>0 for r in changes) for g in ['profile','mafft']})), flush=True)


if __name__ == '__main__':
    main()
