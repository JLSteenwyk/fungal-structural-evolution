#!/usr/bin/env python3
"""Join audited domain clusters to Pfam annotations, preserving boundary alternatives."""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    started = time.monotonic()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)
    def verify():
        if sha(args.plan) != plan_hash:
            raise ValueError('Plan changed')
        for p, digest in plan['pins'].items():
            if sha(p) != digest:
                raise ValueError('Input changed: ' + p)
    verify()
    audit = json.loads(Path(plan['boundary_audit']).read_text())
    registry_audit = json.loads(Path(plan['registry_audit']).read_text())
    if audit['status'] != 'passed_full_domain_boundary_cluster_output_readback':
        raise ValueError('Full boundary audit required')
    if registry_audit['status'] != 'passed_full_structure_domain_registry_readback':
        raise ValueError('Full registry audit required')
    if registry_audit['database_sha256'] != plan['pins'][plan['registry']]:
        raise ValueError('Registry audit binding differs')
    for key in ['pairs', 'intervals']:
        if audit['source_hashes'][plan[key]] != plan['pins'][plan[key]]:
            raise ValueError('Boundary audit binding differs')
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    registry = sqlite3.connect('file:' + plan['registry'] + '?mode=ro', uri=True)
    registry.row_factory = sqlite3.Row
    intervals = {}
    with open(plan['intervals']) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['interval_id'] in intervals:
                raise ValueError('Duplicate interval')
            intervals[r['interval_id']] = (r['model_key'], int(r['start']), int(r['end']))
    db = sqlite3.connect(out / 'cluster_pfam.sqlite')
    db.execute('CREATE TABLE links(boundary TEXT, representative TEXT, interval_id TEXT, model TEXT, hit_id TEXT, pfam_accession TEXT, pfam_clan TEXT, PRIMARY KEY(boundary,model,hit_id))')
    # These sets are built directly from source rows, independently of SQL aggregation.
    stats = defaultdict(lambda: [set(), set(), set()])
    clusters = set(); seen_intervals = set(); pairs = 0
    with gzip.open(plan['pairs'], 'rt') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            model, hit = row['model_key'], row['hit_id']
            segment = registry.execute('SELECT * FROM segments WHERE model_key=? AND hit_id=?', (model, hit)).fetchone()
            if segment is None or segment['pfam_type'] != 'Domain':
                raise ValueError('Missing/non-domain annotation')
            pfam, clan = segment['pfam_accession'], segment['pfam_clan']
            if not pfam or clan is None:
                raise ValueError('Missing annotation identity')
            for boundary in ['alignment', 'envelope']:
                interval, rep = row[boundary + '_interval'], row[boundary + '_cluster']
                if not rep or intervals[interval] != (model, segment[boundary + '_start'], segment[boundary + '_end']):
                    raise ValueError('Boundary coordinates or cluster differ')
                db.execute('INSERT INTO links VALUES(?,?,?,?,?,?,?)', (boundary, rep, interval, model, hit, pfam, clan))
                key = boundary, rep, pfam, clan
                stats[key][0].add(model)
                stats[key][1].add((model, hit))
                stats[key][2].add(interval)
                clusters.add(rep); seen_intervals.add(interval)
            pairs += 1
            if pairs % 100000 == 0:
                db.commit(); print(pairs, 'pairs linked', flush=True)
    db.commit()
    if pairs != audit['candidate_pairs'] or len(intervals) != audit['intervals'] or seen_intervals != set(intervals):
        raise ValueError('Incomplete paired-boundary scope')
    header = ['boundary', 'representative', 'pfam_accession', 'pfam_clan', 'models', 'hits', 'intervals']
    table = out / 'cluster_pfam.tsv'
    groups = 0
    with table.open('w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(header)
        query = '''SELECT boundary,representative,pfam_accession,pfam_clan,
                   COUNT(DISTINCT model),COUNT(*),COUNT(DISTINCT interval_id)
                   FROM links GROUP BY boundary,representative,pfam_accession,pfam_clan
                   ORDER BY boundary,representative,pfam_accession,pfam_clan'''
        for row in db.execute(query):
            key = tuple(row[:4]); expected = tuple(len(v) for v in stats.pop(key))
            if tuple(row[4:]) != expected:
                raise ValueError('SQL/source set aggregation differs')
            w.writerow(row); groups += 1
    if stats:
        raise ValueError('Missing annotation groups')
    # Independently check every serialized row against SQL; no sampling.
    db.execute('CREATE INDEX group_idx ON links(boundary,representative,pfam_accession,pfam_clan)')
    read = 0; multi = defaultdict(int); per_cluster = defaultdict(set)
    with table.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            key = tuple(row[k] for k in header[:4])
            found = db.execute('SELECT COUNT(DISTINCT model),COUNT(*),COUNT(DISTINCT interval_id) FROM links WHERE boundary=? AND representative=? AND pfam_accession=? AND pfam_clan=?', key).fetchone()
            if found != tuple(int(row[k]) for k in header[4:]):
                raise ValueError('Serialized count differs')
            per_cluster[key[:2]].add(key[2]); read += 1
    if read != groups:
        raise ValueError('Serialized scope differs')
    for (boundary, rep), accessions in per_cluster.items():
        if len(accessions) > 1:
            multi[boundary] += 1
    db.commit(); db.close(); registry.close(); verify()
    result = dict(status='complete_domain_cluster_pfam_source_set_sql_agreement',
                  plan_sha256=plan_hash, pairs=pairs, intervals=len(intervals), clusters=len(clusters),
                  boundary_links=2*pairs, annotation_groups=groups,
                  represented_clusters_by_boundary={b: sum(k[0] == b for k in per_cluster) for b in ['alignment','envelope']},
                  multiple_pfam_clusters_by_boundary={b: multi[b] for b in ['alignment','envelope']},
                  elapsed_seconds=time.monotonic()-started,
                  artifacts={p.name: sha(p) for p in [table, out/'cluster_pfam.sqlite']},
                  scope='Full selected Domain/Pfam annotation join. Alignment and envelope alternatives remain separate; counts distinguish source models, model/hit pairs and intervals. SQL aggregates agree with source-derived Python sets and all serialized counts are checked. Shared source annotation and cluster partition; no new annotation searches, independent clustering, confidence qualification, established homology, function or evolutionary events.')
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
