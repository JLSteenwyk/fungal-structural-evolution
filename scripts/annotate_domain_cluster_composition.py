#!/usr/bin/env python3
"""Preserve domain-cluster links to proteins, taxa and both family partitions."""
import argparse
import csv
import hashlib
import json
import shutil
import sqlite3
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8388608), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan_hash = sha(args.plan)
    p = json.loads(args.plan.read_text())

    def verify():
        if sha(args.plan) != plan_hash:
            raise ValueError('Plan changed')
        for name, digest in p['pins'].items():
            if sha(name) != digest:
                raise ValueError('Changed source: '+name)

    verify()
    out = Path(p['output'])
    if out.exists():
        raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free < 100*2**30:
        raise RuntimeError('Insufficient disk reserve')
    out.mkdir()
    c = sqlite3.connect(out/'domain_cluster_composition.sqlite', uri=True)
    c.execute('PRAGMA cache_size=-262144')
    c.execute('ATTACH DATABASE ? AS source', ('file:'+str(Path(p['source_database']).resolve())+'?mode=ro',))
    c.executescript('''
      CREATE TABLE intervals(interval_id TEXT PRIMARY KEY,model TEXT NOT NULL,start INTEGER,end INTEGER);
      CREATE TABLE members(interval_id TEXT PRIMARY KEY,representative TEXT NOT NULL);
    ''')
    with Path(p['intervals']).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            start, end = int(row['start']), int(row['end'])
            if start < 1 or end-start+1 != int(row['residues']):
                raise ValueError('Invalid interval')
            c.execute('INSERT INTO intervals VALUES(?,?,?,?)', (row['interval_id'], row['model_key'], start, end))
    with Path(p['members']).open() as handle:
        c.executemany('INSERT INTO members VALUES(?,?)', ((m,r) for r,m in csv.reader(handle, delimiter='\t')))
    count = lambda sql: c.execute(sql).fetchone()[0]
    if (count('SELECT COUNT(*) FROM intervals'), count('SELECT COUNT(*) FROM members')) != (1078592,1078592):
        raise ValueError('Wrong interval universe')
    if count('SELECT COUNT(*) FROM intervals i LEFT JOIN members m USING(interval_id) WHERE m.interval_id IS NULL'):
        raise ValueError('Interval without membership')
    if count('SELECT COUNT(*) FROM members m LEFT JOIN members r ON m.representative=r.interval_id WHERE r.interval_id IS NULL OR r.interval_id!=r.representative'):
        raise ValueError('Invalid representative')
    c.execute('CREATE INDEX intervals_model ON intervals(model)')
    c.execute('CREATE INDEX members_representative ON members(representative)')
    c.execute('CREATE TABLE proteins(native_gene_id TEXT PRIMARY KEY,taxon_id TEXT,protein_id TEXT,model TEXT)')
    c.execute('''INSERT INTO proteins SELECT l.native_gene_id,l.taxon_id,l.protein_id,l.model
      FROM source.links l WHERE EXISTS(SELECT 1 FROM intervals i WHERE i.model=l.model)''')
    c.execute('CREATE INDEX proteins_model ON proteins(model)')
    if count('SELECT COUNT(DISTINCT model) FROM proteins') != 427255:
        raise ValueError('Incomplete source-model join')
    c.execute('CREATE TABLE family_links(guide TEXT,native_gene_id TEXT,family TEXT,PRIMARY KEY(guide,native_gene_id))')
    c.execute('''INSERT INTO family_links SELECT f.guide,f.native_gene_id,f.family
      FROM source.family_links f JOIN proteins p USING(native_gene_id)''')
    proteins = count('SELECT COUNT(*) FROM proteins')
    for guide in ['profile','mafft']:
        if c.execute('SELECT COUNT(*) FROM family_links WHERE guide=?',(guide,)).fetchone()[0] != proteins:
            raise ValueError('Incomplete family join')
    c.execute('CREATE TABLE cluster_proteins(representative TEXT,native_gene_id TEXT,PRIMARY KEY(representative,native_gene_id))')
    c.execute('''INSERT INTO cluster_proteins SELECT DISTINCT m.representative,p.native_gene_id
      FROM members m JOIN intervals i USING(interval_id) JOIN proteins p USING(model)''')
    c.execute('CREATE TABLE cluster_intervals AS SELECT representative,COUNT(*) AS intervals,COUNT(DISTINCT model) AS models FROM members JOIN intervals USING(interval_id) GROUP BY representative')
    c.execute('CREATE UNIQUE INDEX cluster_intervals_rep ON cluster_intervals(representative)')
    c.commit()
    print('All source joins complete; aggregating clusters', flush=True)
    clusters = 0
    with (out/'domain_cluster_composition.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['representative','intervals','models','proteins','taxa','profile_families','mafft_families'])
        query = '''SELECT cp.representative,ci.intervals,ci.models,COUNT(*),COUNT(DISTINCT p.taxon_id),COUNT(DISTINCT f.family),COUNT(DISTINCT g.family)
          FROM cluster_proteins cp JOIN proteins p USING(native_gene_id)
          JOIN cluster_intervals ci USING(representative)
          JOIN family_links f ON f.native_gene_id=p.native_gene_id AND f.guide='profile'
          JOIN family_links g ON g.native_gene_id=p.native_gene_id AND g.guide='mafft'
          GROUP BY cp.representative ORDER BY cp.representative'''
        for row in c.execute(query):
            writer.writerow(row)
            clusters += 1
    if clusters != 70537 or count('SELECT SUM(intervals) FROM cluster_intervals') != 1078592:
        raise ValueError('Cluster totals differ')
    links = count('SELECT COUNT(*) FROM cluster_proteins')
    if c.execute('PRAGMA integrity_check').fetchall() != [('ok',)]:
        raise ValueError('Database integrity failed')
    c.close()
    verify()
    receipt = dict(status='complete_domain_cluster_composition_pending_independent_readback',
                   plan_sha256=plan_hash,intervals=1078592,models=427255,proteins=proteins,
                   clusters=clusters,cluster_protein_links=links,
                   artifacts={f.name:sha(f) for f in out.iterdir()},
                   scope='Both family partitions and all taxon/protein links retained. Alternative boundaries remain linked to their original interval identities. Counts are descriptive; not homology, orthology, independent observations, gains/losses or evolutionary events. Independent source and aggregation readback pending.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__ == '__main__':
    main()
