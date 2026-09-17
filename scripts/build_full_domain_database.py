#!/usr/bin/env python3
"""Join complete Pfam partitions and representative-protein links without filtering hits."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    a = ap.parse_args()
    plan = load(a.plan)
    for path, digest in plan['pins'].items():
        if sha(path) != digest:
            raise ValueError('Changed pinned input: ' + path)
    output = Path(plan['output'])
    if output.exists():
        raise FileExistsError(output)
    if shutil.disk_usage(output.parent).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
        raise ValueError('Insufficient free space')
    full = Path(plan['full_inputs'])
    marker = Path(plan['marker_annotations'])
    catalog = Path(plan['catalog'])
    fr, mr, cr = load(full / 'receipt.json'), load(marker / 'receipt.json'), load(catalog / 'receipt.json')
    links = full / 'protein_links.tsv'
    marker_hits = marker / 'raw_annotated_hits.tsv'
    assert sha(links) == fr['artifacts'][links.name]
    assert sha(marker_hits) == mr['artifacts'][marker_hits.name]
    assert sha(catalog / 'annotation_shards.tsv') == cr['artifacts']['annotation_shards.tsv']
    with (catalog / 'annotation_shards.tsv').open() as handle:
        shards = list(csv.DictReader(handle, delimiter='\t'))
    assert len(shards) == cr['chunks'] == 64
    assert len({r['chunk'] for r in shards}) == len(shards)
    for row in shards:
        assert sha(row['annotation_path']) == row['annotation_sha256']
        assert sha(row['shard_receipt']) == row['shard_receipt_sha256']
    output.mkdir(parents=True)
    dbpath = output / 'domains.sqlite'
    conn = sqlite3.connect(dbpath)
    conn.execute('PRAGMA cache_size=-262144')
    conn.execute('PRAGMA synchronous=NORMAL')
    conn.execute('CREATE TABLE proteins(taxon_id TEXT, protein_id TEXT, sequence_id TEXT, query_source TEXT, PRIMARY KEY(taxon_id,protein_id))')
    with links.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        batch = []
        for row in reader:
            batch.append(tuple(row[k] for k in ['taxon_id', 'protein_id', 'sequence_id', 'query_source']))
            if len(batch) == 10000:
                conn.executemany('INSERT INTO proteins VALUES(?,?,?,?)', batch)
                batch.clear()
        conn.executemany('INSERT INTO proteins VALUES(?,?,?,?)', batch)
    conn.execute('CREATE INDEX proteins_by_sequence ON proteins(sequence_id)')
    conn.execute('CREATE TABLE queries(sequence_id TEXT PRIMARY KEY, search_partition TEXT NOT NULL)')
    conn.execute('INSERT INTO queries SELECT DISTINCT sequence_id,query_source FROM proteins')
    # Keep the four marker-only queries explicit; they have no representative-protein link.
    mlinks = Path(plan['marker_inputs']) / 'protein_links.tsv'
    mir = load(Path(plan['marker_inputs']) / 'receipt.json')
    assert sha(mlinks) == mir['artifacts'][mlinks.name]
    marker_ids = set()
    with mlinks.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            marker_ids.add(row['sequence_id'])
    assert len(marker_ids) == mir['unique_sequences']
    conn.executemany('INSERT OR IGNORE INTO queries VALUES(?,?)', [(x, 'marker-inputs-v1') for x in sorted(marker_ids)])
    assert all(conn.execute('SELECT search_partition FROM queries WHERE sequence_id=?', (x,)).fetchone()[0] == 'marker-inputs-v1' for x in marker_ids)
    conn.commit()
    print('Protein links and query universe imported', flush=True)
    with gzip.open(shards[0]['annotation_path'], 'rt') as handle:
        fields = next(csv.reader(handle, delimiter='\t'))
    assert all(k.replace('_', '').isalnum() for k in fields)
    assert 'search_partition' in fields and 'hit_id' in fields
    conn.execute('CREATE TABLE hits(hit_number INTEGER PRIMARY KEY,source_shard TEXT NOT NULL,' + ','.join('"'+k+'" TEXT NOT NULL' for k in fields) + ')')
    sql = 'INSERT INTO hits(source_shard,' + ','.join('"'+k+'"' for k in fields) + ') VALUES(' + ','.join('?' for _ in range(len(fields)+1)) + ')'
    manifest = []
    sources = [(r['chunk'], Path(r['annotation_path']), 'additional_full_proteome', int(r['hit_rows'])) for r in shards]
    sources.append(('marker', marker_hits, 'marker-inputs-v1', mr['raw_hits']))
    for label, path, partition, expected in sources:
        opener = gzip.open if path.suffix == '.gz' else open
        count = 0
        with opener(path, 'rt') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            assert set(reader.fieldnames) == set(fields) - ({'search_partition'} if label == 'marker' else set())
            batch = []
            for row in reader:
                if label == 'marker':
                    row['search_partition'] = partition
                assert row['search_partition'] == partition
                batch.append((label,) + tuple(row[k] for k in fields))
                count += 1
                if len(batch) == 10000:
                    conn.executemany(sql, batch)
                    batch.clear()
            conn.executemany(sql, batch)
        assert count == expected
        conn.commit()
        manifest.append(dict(source_shard=label, path=str(path), sha256=sha(path), rows=count, search_partition=partition))
        print(label, count, flush=True)
        if dbpath.stat().st_size > plan['resources']['output_allowance_gib'] * 2**30:
            raise ValueError('Database output allowance exceeded')
    conn.execute('CREATE UNIQUE INDEX unique_source_hit ON hits(search_partition,hit_id)')
    conn.execute('CREATE INDEX hits_by_sequence ON hits(sequence_id)')
    assert conn.execute('SELECT COUNT(*) FROM hits h LEFT JOIN queries q ON q.sequence_id=h.sequence_id WHERE q.sequence_id IS NULL OR q.search_partition!=h.search_partition').fetchone()[0] == 0
    assert conn.execute('SELECT COUNT(*) FROM proteins').fetchone()[0] == 5815847
    assert conn.execute('SELECT COUNT(DISTINCT taxon_id) FROM proteins').fetchone()[0] == 526
    assert conn.execute('SELECT COUNT(*) FROM queries').fetchone()[0] == 5713603
    assert conn.execute('SELECT COUNT(*) FROM hits').fetchone()[0] == 8103610
    conn.execute('CREATE VIEW protein_hits AS SELECT p.taxon_id,p.protein_id,h.* FROM proteins p JOIN hits h USING(sequence_id)')
    conn.execute('CREATE VIEW protein_detection AS SELECT p.*, CASE WHEN EXISTS(SELECT 1 FROM hits h WHERE h.sequence_id=p.sequence_id) THEN "raw_hits_require_architecture_review" ELSE "no_GA_hit_not_proven_absence" END AS detection_status FROM proteins p')
    assert conn.execute('PRAGMA integrity_check').fetchall() == [('ok',)]
    conn.commit()
    conn.close()
    (output / 'source_shards.json').write_text(json.dumps(manifest, indent=2) + '\n')
    result = dict(status='complete_domain_partition_database_requires_independent_readback',
                  proteins=5815847, taxa=526, queries=5713603, marker_only_queries=4,
                  additional_hit_rows=7939960, marker_hit_rows=163650, total_hit_rows=8103610,
                  plan_sha256=sha(a.plan), script_sha256=sha(Path(__file__)),
                  artifacts={p.name: sha(p) for p in [dbpath, output / 'source_shards.json']},
                  scope='Every source hit and full-proteome protein link retained. Text values preserve source scores, coordinates and partition-specific E-values exactly. Views expose protein joins and explicitly qualified no-hit states; no overlap resolution or biological absence inference.')
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
