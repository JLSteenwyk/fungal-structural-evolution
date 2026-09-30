#!/usr/bin/env python3
"""Independently check every pair/source reuse disposition with SQL joins."""
import argparse
import csv
import hashlib
import json
import sqlite3
from collections import Counter
from pathlib import Path

from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    root = Path(plan['output'])
    rp = root / 'receipt.json'
    receipt = json.loads(rp.read_text())
    table = root / 'pair_reuse_candidates.tsv'
    bindings = {str(args.plan): sha(args.plan), **plan['pins'], str(rp): sha(rp), str(table): sha(table)}

    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed source: ' + path)

    verify()
    if receipt['status'] != 'complete_full_pair_union_catalog_reuse_screen_not_authorization':
        raise ValueError('Incomplete producer')
    if receipt['source_hashes'] != {str(args.plan): sha(args.plan), **plan['pins']} or receipt['artifacts'][table.name] != sha(table):
        raise ValueError('Receipt lineage differs')
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE models(era,label,id,version INTEGER,coordinate,sequence,length INTEGER,PRIMARY KEY(era,label,id,version))')
    db.execute('CREATE TABLE pairs(era,label,key,a,av INTEGER,b,bv INTEGER,PRIMARY KEY(era,label,key))')
    sources = {'old': {}, 'new': {}}
    for era in ['old', 'new']:
        for spec in plan[era]:
            label, source = spec['label'], Path(spec['root'])
            if label in sources[era]:
                raise ValueError('Repeated source label')
            r = json.loads((source / 'receipt.json').read_text())
            if r['status'] != spec['status']:
                raise ValueError('Incomplete source')
            for name in ['receipt.json', 'models.jsonl', 'model_pairs.tsv']:
                if str(source / name) not in bindings:
                    raise ValueError('Unpinned source')
                if name != 'receipt.json' and sha(source / name) != r['artifacts'][name]:
                    raise ValueError('Changed artifact')
            n = 0
            with (source / 'models.jsonl').open() as handle:
                for line in handle:
                    m = json.loads(line)
                    if type(m['version']) is not int or m['length'] <= 0:
                        raise ValueError('Invalid native model descriptor')
                    db.execute('INSERT INTO models VALUES(?,?,?,?,?,?,?)', (era, label, m['model_id'], m['version'], m['sha256'], m['sequence_sha256'], m['length']))
                    n += 1
            p = 0
            with (source / 'model_pairs.tsv').open() as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    a, b = (row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))
                    if a == b:
                        raise ValueError('Identical endpoints')
                    # Explicit tuple comparison preserves numeric version order.
                    left, right = (a, b) if a < b else (b, a)
                    canonical = json.dumps([list(left), list(right)], separators=(',', ':'))
                    if hashlib.sha256(canonical.encode()).hexdigest() != row['pair_key']:
                        raise ValueError('Invalid pair digest')
                    db.execute('INSERT INTO pairs VALUES(?,?,?,?,?,?,?)', (era, label, row['pair_key'], *left, *right))
                    p += 1
            if n != r[spec['model_count_field']] or p != r[spec['pair_count_field']]:
                raise ValueError('Native counts differ')
            sources[era][label] = dict(models=n, pairs=p)
    missing = db.execute('SELECT COUNT(*) FROM pairs p LEFT JOIN models a ON p.era=a.era AND p.label=a.label AND p.a=a.id AND p.av=a.version LEFT JOIN models b ON p.era=b.era AND p.label=b.label AND p.b=b.id AND p.bv=b.version WHERE a.id IS NULL OR b.id IS NULL').fetchone()[0]
    if missing:
        raise ValueError('Pair endpoint absent from native models')
    db.execute('CREATE TABLE descriptors AS SELECT p.*,a.coordinate AS ac,a.sequence AS ass,a.length AS al,b.coordinate AS bc,b.sequence AS bss,b.length AS bl FROM pairs p JOIN models a ON p.era=a.era AND p.label=a.label AND p.a=a.id AND p.av=a.version JOIN models b ON p.era=b.era AND p.label=b.label AND p.b=b.id AND p.bv=b.version')
    db.execute('CREATE INDEX descriptor_keys ON descriptors(era,key)')
    db.execute('CREATE TABLE reviewed(key PRIMARY KEY)')
    counts, per_source, checked = Counter(), Counter(), 0
    with table.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['pair_key']
            db.execute('INSERT INTO reviewed VALUES(?)', (key,))
            native = db.execute("SELECT label,a,av,b,bv,ac,ass,al,bc,bss,bl FROM descriptors WHERE era='new' AND key=? ORDER BY label", (key,)).fetchall()
            if not native or any(r[1:] != native[0][1:] for r in native):
                raise ValueError('Missing or conflicting new source descriptors')
            fields = (row['model_a'], int(row['version_a']), row['model_b'], int(row['version_b']))
            if fields != native[0][1:5]:
                raise ValueError('Exported endpoints differ')
            previous = db.execute("SELECT label,a,av,b,bv,ac,ass,al,bc,bss,bl FROM descriptors WHERE era='old' AND key=? ORDER BY label", (key,)).fetchall()
            if any(r[1:5] != fields for r in previous):
                raise ValueError('Previous endpoint collision')
            matching = [r[0] for r in previous if r[5:] == native[0][5:]]
            changed = [r[0] for r in previous if r[5:] != native[0][5:]]
            status = ('matching_catalog_sources_pending_input_and_result_checks' if matching else
                      'changed_catalog_sources_require_alignment' if previous else 'new_pair_requires_alignment')
            expected = [[r[0] for r in native], matching, changed, status]
            actual = [json.loads(row['new_sources']), json.loads(row['matching_old_sources']), json.loads(row['changed_old_sources']), row['disposition']]
            if actual != expected:
                raise ValueError('False reuse disposition or missing source')
            counts[status] += 1
            for label in expected[0]:
                per_source[label + ':' + status] += 1
            checked += 1
    new_union = db.execute("SELECT COUNT(DISTINCT key) FROM pairs WHERE era='new'").fetchone()[0]
    old_union = db.execute("SELECT COUNT(DISTINCT key) FROM pairs WHERE era='old'").fetchone()[0]
    old_absent = db.execute("SELECT COUNT(DISTINCT key) FROM pairs WHERE era='old' AND key NOT IN (SELECT key FROM pairs WHERE era='new')").fetchone()[0]
    if checked != new_union or dict(counts) != receipt['counts'] or dict(per_source) != receipt['per_new_source_counts']:
        raise ValueError('Incomplete saved union or summaries')
    for field, value in [('new_union_pairs', new_union), ('old_union_pairs', old_union), ('old_union_pairs_absent', old_absent), ('old_sources', sources['old']), ('new_sources', sources['new'])]:
        if receipt[field] != value:
            raise ValueError('False source summary: ' + field)
    verify()
    result = dict(status='passed_full_pair_union_reuse_candidate_sql_readback',
                  producer_receipt_sha256=sha(rp), new_union_pairs=new_union, old_union_pairs=old_union,
                  counts=dict(counts), per_new_source_counts=dict(per_source), source_hashes=bindings,
                  scientific_eligibility=False,
                  scope='Every native model, pair and source membership loaded independently into SQL; full distinct pair union and all per-endpoint coordinate/sequence/length matching dispositions reconstructed. Rehashed false dispositions and incomplete unions rejected. This is a catalog-level planning check, not alignment reuse authorization or numerical/biological qualification.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(status=result['status'], new_union_pairs=new_union, counts=dict(counts))), flush=True)


if __name__ == '__main__':
    main()
