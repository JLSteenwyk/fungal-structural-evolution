#!/usr/bin/env python3
"""Compare every database row to immutable source tables and check all taxon joins."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    plan = json.loads(a.plan.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['plan_sha256'] == sha(a.plan)
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
    conn = sqlite3.connect('file:' + str((root / 'domains.sqlite').resolve()) + '?mode=ro', uri=True)
    fields = [r[1] for r in conn.execute('PRAGMA table_info(hits)')][2:]
    actual = iter(conn.execute('SELECT * FROM hits ORDER BY hit_number'))
    counts = Counter()
    n = 0
    for source in json.loads((root / 'source_shards.json').read_text()):
        path = Path(source['path'])
        assert sha(path) == source['sha256']
        opener = gzip.open if path.suffix == '.gz' else open
        lines = 0
        with opener(path, 'rt') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if source['source_shard'] == 'marker':
                    row['search_partition'] = 'marker-inputs-v1'
                n += 1
                lines += 1
                expected = (n, source['source_shard']) + tuple(row[k] for k in fields)
                assert next(actual, None) == expected, (source['source_shard'], lines)
                counts[row['sequence_id']] += 1
        assert lines == source['rows']
    assert next(actual, None) is None and n == receipt['total_hit_rows']
    summaries = defaultdict(lambda: [0, 0, 0])
    rows = conn.execute('SELECT taxon_id,protein_id,sequence_id,query_source FROM proteins ORDER BY rowid')
    with (Path(plan['full_inputs']) / 'protein_links.tsv').open() as handle:
        for source, stored in itertools.zip_longest(csv.DictReader(handle, delimiter='\t'), rows):
            assert source is not None and stored is not None
            assert tuple(source[k] for k in ['taxon_id','protein_id','sequence_id','query_source']) == stored
            summary = summaries[source['taxon_id']]
            summary[0] += 1
            summary[1] += counts[source['sequence_id']] > 0
            summary[2] += counts[source['sequence_id']]
    sql = '''SELECT p.taxon_id,COUNT(*),SUM(CASE WHEN c.n IS NULL THEN 0 ELSE 1 END),SUM(COALESCE(c.n,0))
             FROM proteins p LEFT JOIN (SELECT sequence_id,COUNT(*) AS n FROM hits GROUP BY sequence_id) c
             ON p.sequence_id=c.sequence_id GROUP BY p.taxon_id'''
    joined = {r[0]: list(r[1:]) for r in conn.execute(sql)}
    assert joined == dict(summaries)
    assert len(joined) == 526 and sum(r[0] for r in joined.values()) == 5815847
    assert conn.execute('PRAGMA integrity_check').fetchall() == [('ok',)]
    conn.close()
    table = a.output.with_name('taxon_domain_detection.tsv')
    if table.exists():
        raise FileExistsError(table)
    with table.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['taxon_id','proteins','proteins_with_raw_hits','proteins_without_GA_hits','protein_linked_raw_hit_rows'])
        for taxon, (total, detected, hits) in sorted(joined.items()):
            writer.writerow([taxon,total,detected,total-detected,hits])
    result = dict(status='passed_complete_domain_database_source_readback', proteins=5815847,
                  taxa=526, hit_rows=n, proteins_with_raw_hits=sum(r[1] for r in joined.values()),
                  protein_linked_raw_hit_rows=sum(r[2] for r in joined.values()),
                  database_receipt_sha256=sha(root / 'receipt.json'), script_sha256=sha(Path(__file__)),
                  artifacts={table.name: sha(table)},
                  scope='Every hit field and protein link compared exactly to source tables; all 526 taxon join totals checked against independently counted source hits. No architecture resolution or proof of domain absence.')
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
