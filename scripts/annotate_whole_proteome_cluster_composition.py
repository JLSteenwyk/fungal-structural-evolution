#!/usr/bin/env python3
"""Join the full structural partition to audited protein and family identities."""
import argparse,csv,hashlib,json,shutil,sqlite3
from pathlib import Path

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8388608),b''): h.update(block)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    a=ap.parse_args();p=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Plan changed')
        for f,h in p['pins'].items():
            if sha(f)!=h:raise ValueError('Changed input: '+f)
    verify();out=Path(p['output'])
    if out.exists():raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free<100*2**30:raise RuntimeError('Disk reserve')
    out.mkdir();c=sqlite3.connect(out/'cluster_composition.sqlite',uri=True)
    c.execute('PRAGMA cache_size=-262144')
    for alias,key in [('source','structures'),('family','families')]:
        c.execute('ATTACH DATABASE ? AS '+alias,('file:'+str(Path(p[key]).resolve())+'?mode=ro',))
    c.execute('CREATE TABLE members(model TEXT PRIMARY KEY,representative TEXT NOT NULL)')
    with Path(p['members']).open() as f:
        c.executemany('INSERT INTO members VALUES(?,?)',((m,r) for r,m in csv.reader(f,delimiter='\t')))
    c.execute('CREATE INDEX members_representative ON members(representative)')
    n=c.execute('SELECT COUNT(*) FROM members').fetchone()[0]
    if n!=1290278:raise ValueError('Wrong model universe')
    if c.execute('SELECT COUNT(*) FROM members m LEFT JOIN members r ON m.representative=r.model WHERE r.model IS NULL OR r.model!=r.representative').fetchone()[0]:raise ValueError('Invalid representative')
    c.execute('CREATE TABLE links(native_gene_id TEXT PRIMARY KEY,taxon_id TEXT,protein_id TEXT,model TEXT,representative TEXT)')
    c.execute('''INSERT INTO links SELECT s.native_gene_id,s.taxon_id,s.protein_id,m.model,m.representative
      FROM source.structures s JOIN members m ON m.model=s.model_id||'-v'||s.version''')
    links=c.execute('SELECT COUNT(*) FROM links').fetchone()[0]
    if links!=1319513 or links!=c.execute('SELECT COUNT(*) FROM source.structures').fetchone()[0]:raise ValueError('Unmapped proteins')
    if c.execute('SELECT COUNT(DISTINCT model) FROM links').fetchone()[0]!=n:raise ValueError('Unmapped models')
    c.execute('CREATE INDEX links_representative ON links(representative)')
    c.execute('CREATE TABLE family_links(guide TEXT,native_gene_id TEXT,family TEXT,PRIMARY KEY(guide,native_gene_id))')
    for guide in ['profile','mafft']:
        c.execute('INSERT INTO family_links SELECT a.guide,a.native_gene_id,a.family FROM family.assignments a JOIN links l USING(native_gene_id) WHERE a.guide=?',(guide,))
        if c.execute('SELECT COUNT(*) FROM family_links WHERE guide=?',(guide,)).fetchone()[0]!=links:raise ValueError('Incomplete family join')
    c.commit()
    with (out/'cluster_composition.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['representative','models','proteins','taxa','profile_families','mafft_families'])
        count=0;total_models=total_proteins=0
        for row in c.execute('''SELECT l.representative,COUNT(DISTINCT l.model),COUNT(*),COUNT(DISTINCT l.taxon_id),COUNT(DISTINCT p.family),COUNT(DISTINCT m.family)
          FROM links l JOIN family_links p ON p.native_gene_id=l.native_gene_id AND p.guide='profile'
          JOIN family_links m ON m.native_gene_id=l.native_gene_id AND m.guide='mafft'
          GROUP BY l.representative ORDER BY l.representative'''):
            w.writerow(row);count+=1;total_models+=row[1];total_proteins+=row[2]
    if (count,total_models,total_proteins)!=(400127,n,links):raise ValueError('Composition totals differ')
    if c.execute('PRAGMA integrity_check').fetchall()!=[('ok',)]:raise ValueError('Database integrity')
    c.close();verify()
    r=dict(status='complete_cluster_composition_pending_independent_readback',plan_sha256=ph,models=n,protein_links=links,clusters=count,artifacts={f.name:sha(f) for f in out.iterdir()},scope='Complete model/protein/taxon joins and both family partitions. Counts describe represented proteins, not independent observations or verified species. No confidence qualification, homology, orthology, novelty, enrichment or evolutionary event inference. Independent source-join and aggregate readback remains required.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':main()
