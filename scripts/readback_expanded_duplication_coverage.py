#!/usr/bin/env python3
"""Check every refreshed event against its source and independent SQL model joins."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import sqlite3
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    plan = json.loads(args.plan.read_text())
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['plan_sha256'] == digest(args.plan)
    pins = dict(plan['pins'])
    pins[str(args.plan)] = digest(args.plan)
    pins[str(root / 'receipt.json')] = digest(root / 'receipt.json')
    for name, expected in receipt['artifacts'].items():
        pins[str(root / name)] = expected
    for path, expected in pins.items():
        assert digest(path) == expected, path
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE models (gene TEXT PRIMARY KEY, model TEXT NOT NULL)')
    with open(plan['links']) as stream:
        db.executemany('INSERT INTO models VALUES (?,?)',
                       ((r['taxon_id'] + '_' + r['protein_id'], r['model_id'])
                        for r in csv.DictReader(stream, delimiter='\t')))
    assert db.execute('SELECT COUNT(*) FROM models').fetchone()[0] == receipt['catalog_links']
    db.execute('CREATE TABLE events (guide TEXT, family TEXT, taxon TEXT, a TEXT, b TEXT, old_class TEXT, new_a TEXT, new_b TEXT, new_class TEXT, same INTEGER)')
    n = 0
    with gzip.open(plan['events'], 'rt') as original, gzip.open(root / 'event_coverage.tsv.gz', 'rt') as refreshed:
        old = csv.DictReader(original, delimiter='\t')
        new = csv.DictReader(refreshed, delimiter='\t')
        batch = []
        for a, b in itertools.zip_longest(old, new):
            assert a is not None and b is not None
            assert all(b[k] == v for k, v in a.items()), 'Changed event identity/source data'
            batch.append((b['guide'], b['family'], b['taxon_id'], b['gene_a'], b['gene_b'], b['coverage_class'], b['new_model_a'], b['new_model_b'], b['new_coverage_class'], int(b['new_same_model'])))
            n += 1
            if len(batch) == 10000:
                db.executemany('INSERT INTO events VALUES (?,?,?,?,?,?,?,?,?,?)', batch)
                batch.clear()
        db.executemany('INSERT INTO events VALUES (?,?,?,?,?,?,?,?,?,?)', batch)
    assert n == receipt['event_rows']
    bad = db.execute('''SELECT COUNT(*) FROM events e
        LEFT JOIN models a ON e.a=a.gene LEFT JOIN models b ON e.b=b.gene
        WHERE e.new_a != COALESCE(a.model,'') OR e.new_b != COALESCE(b.model,'')
        OR e.new_class != CASE WHEN a.model IS NOT NULL AND b.model IS NOT NULL THEN 'both_models'
                              WHEN a.model IS NOT NULL OR b.model IS NOT NULL THEN 'one_model' ELSE 'neither_model' END
        OR e.same != CASE WHEN a.model IS NOT NULL AND a.model=b.model THEN 1 ELSE 0 END''').fetchone()[0]
    assert bad == 0, bad
    summaries = {}
    for guide, in db.execute('SELECT DISTINCT guide FROM events'):
        counts = {'events': db.execute('SELECT COUNT(*) FROM events WHERE guide=?', (guide,)).fetchone()[0]}
        for prefix, col in [('old', 'old_class'), ('new', 'new_class')]:
            for label in ['neither_model', 'one_model', 'both_models']:
                counts[prefix + '_' + label] = db.execute('SELECT COUNT(*) FROM events WHERE guide=? AND ' + col + '=?', (guide, label)).fetchone()[0]
        for key, expression in [
            ('new_same_model', 'SUM(same)'),
            ('newly_both', "SUM(old_class!='both_models' AND new_class='both_models')"),
            ('lost_both', "SUM(old_class='both_models' AND new_class!='both_models')"),
            ('taxa_with_both', "COUNT(DISTINCT CASE WHEN new_class='both_models' THEN taxon END)"),
            ('families_with_both', "COUNT(DISTINCT CASE WHEN new_class='both_models' THEN family END)")]:
            counts[key] = db.execute('SELECT ' + expression + ' FROM events WHERE guide=?', (guide,)).fetchone()[0]
        summaries[guide] = counts
    assert summaries == receipt['guides']
    for path, expected in pins.items():
        assert digest(path) == expected, path
    result = dict(status='passed_full_expanded_duplication_event_and_sql_join_readback',
                  event_rows=n, guides=summaries, sources=pins,
                  auditor_sha256=digest(__file__),
                  scope='All event fields preserved; all model assignments and coverage classes independently joined with SQL; every aggregate recomputed. Uses previously audited event selection and catalog identity, does not repeat native event extraction or coordinate validation. No new confidence qualification or evolutionary inference.')
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
