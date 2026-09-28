#!/usr/bin/env python3
"""Inventory model/version pairs for resource planning, without launching alignments."""
import argparse
import csv
import gzip
import hashlib
import json
import sqlite3
from collections import Counter
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--readback', type=Path, required=True)
    ap.add_argument('--old-queue-evidence', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    audit = json.loads(a.readback.read_text())
    assert audit['status'] == 'passed_full_expanded_duplication_event_and_sql_join_readback'
    old = json.loads(a.old_queue_evidence.read_text())
    rp = Path(old['receipt_path'])
    assert sha(rp) == old['receipt_sha256']
    receipt = json.loads(rp.read_text())
    assert receipt == old['receipt']
    pair_path = rp.parent / 'model_pairs.tsv'
    sources = dict(audit['sources'])
    sources.update({str(a.readback): sha(a.readback), str(a.old_queue_evidence): sha(a.old_queue_evidence),
                    str(rp): old['receipt_sha256'], str(pair_path): receipt['artifacts'][pair_path.name]})
    for path, expected in sources.items():
        assert sha(path) == expected, path
    event_path, = [Path(p) for p in sources if p.endswith('/event_coverage.tsv.gz')]
    link_path, = [Path(p) for p in sources if p.endswith('/protein_model_links.tsv')]
    links = {}
    with link_path.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            key = r['taxon_id'] + '_' + r['protein_id']
            assert key not in links
            links[key] = (r['model_id'], r['version'])
    old_pairs = set()
    with pair_path.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            pair = tuple(sorted([(r['model_a'], r['version_a']), (r['model_b'], r['version_b'])]))
            assert pair not in old_pairs and pair[0] != pair[1]
            old_pairs.add(pair)
    assert len(old_pairs) == receipt['unique_distinct_model_pairs']
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE events (guide TEXT, ma TEXT, va TEXT, mb TEXT, vb TEXT)')
    counts = Counter()
    pair_counts = Counter()
    n = 0
    with gzip.open(event_path, 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            n += 1
            if r['new_coverage_class'] != 'both_models':
                continue
            ma, mb = links[r['gene_a']], links[r['gene_b']]
            assert (ma[0], mb[0]) == (r['new_model_a'], r['new_model_b'])
            pair = tuple(sorted([ma, mb]))
            counts[(r['guide'], 'both')] += 1
            counts[(r['guide'], 'same')] += int(ma == mb)
            db.execute('INSERT INTO events VALUES (?,?,?,?,?)', (r['guide'], *pair[0], *pair[1]))
            if ma != mb:
                pair_counts[pair] += 1
    assert n == audit['event_rows']
    for guide, c in audit['guides'].items():
        assert counts[(guide, 'both')] == c['new_both_models']
        assert counts[(guide, 'same')] == c['new_same_model']
    # A separate SQL aggregation checks every pair and its multiplicity.
    sql_counts = {((ma, va), (mb, vb)): k for ma, va, mb, vb, k in db.execute(
        'SELECT ma,va,mb,vb,COUNT(*) FROM events WHERE ma != mb OR va != vb GROUP BY ma,va,mb,vb')}
    assert dict(pair_counts) == sql_counts
    pairs = set(pair_counts)
    common = pairs & old_pairs
    new = pairs - old_pairs
    a.output.mkdir(parents=True)
    table = a.output / 'model_version_pairs.tsv'
    with table.open('w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['model_a', 'version_a', 'model_b', 'version_b', 'event_links', 'in_old_queue'])
        for pair in sorted(pairs):
            w.writerow([*pair[0], *pair[1], pair_counts[pair], int(pair in old_pairs)])
    observed = {}
    with table.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            pair = ((r['model_a'], r['version_a']), (r['model_b'], r['version_b']))
            assert pair not in observed and int(r['in_old_queue']) == int(pair in old_pairs)
            observed[pair] = int(r['event_links'])
    assert observed == sql_counts
    for path, expected in sources.items():
        assert sha(path) == expected, path
    result = dict(status='complete_model_version_pair_resource_inventory', source_sha256=sources,
                  script_sha256=sha(__file__), event_rows=n, distinct_pair_event_links=sum(pair_counts.values()),
                  same_model_event_links=sum(v for (g, k), v in counts.items() if k == 'same'),
                  old_distinct_pairs=len(old_pairs), expanded_distinct_pairs=len(pairs),
                  model_version_pairs_shared_with_old_queue=len(common), new_model_version_pairs=len(new),
                  old_pairs_not_in_expanded_set=len(old_pairs - pairs),
                  new_directed_comparisons_one_mask=2 * len(new),
                  new_directed_comparisons_two_masks=4 * len(new),
                  artifacts={table.name: sha(table)},
                  validation='Every pair multiplicity checked against SQL GROUP BY and every serialized row checked.',
                  scope='Resource planning only. Expanded tree-node correspondence, coordinate/confidence validation, masks and alignment eligibility remain required. Matching model IDs/versions do not establish identical coordinate bytes or permit automatic reuse of old alignments. All event associations remain in the source table; pair reuse does not imply independent events. No new structural comparisons launched.')
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'scope', 'validation')}, indent=2))


if __name__ == '__main__':
    main()
