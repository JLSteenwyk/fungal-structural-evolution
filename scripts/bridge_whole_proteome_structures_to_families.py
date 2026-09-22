#!/usr/bin/env python3
"""Attach complete structural coverage to both full gene-family partitions."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sqlite3
import time
import psutil
from catalog_whole_proteome_structures import sha


def family_rows(connection,guide):
    return connection.execute('''
      SELECT a.family, COUNT(*) AS proteins, COUNT(DISTINCT p.taxon_id) AS taxa,
             SUM(s.native_gene_id IS NOT NULL) AS structure_proteins,
             COUNT(DISTINCT CASE WHEN s.native_gene_id IS NOT NULL THEN p.taxon_id END) AS structure_taxa,
             COUNT(DISTINCT s.sequence_sha256) AS structure_sequences,
             COUNT(DISTINCT s.model_path) AS structure_models
      FROM bridge.assignments a INDEXED BY assignments_family
      JOIN bridge.proteins p ON p.native_gene_id=a.native_gene_id
      LEFT JOIN structures s ON s.native_gene_id=p.native_gene_id
      WHERE a.guide=? GROUP BY a.family ORDER BY a.family''',(guide,))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan_sha=sha(args.plan)
    def verify():
        if sha(args.plan)!=plan_sha:raise ValueError('Changed plan')
        for path,digest in plan['pins'].items():
            if sha(path)!=digest:raise ValueError('Changed pinned dependency')
    verify();dependency=plan['predecessor']
    while True:
        try:
            p=psutil.Process(dependency['pid'])
            live=p.create_time()==dependency['create_time'] and p.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify();catalog=Path(plan['catalog']);receipt=json.loads((catalog/'receipt.json').read_text())
    audit=json.loads(Path(plan['readback']).read_text())
    if (audit['status']!='passed_full_proteome_sequence_and_model_selection_readback'
            or audit['producer_receipt_sha256']!=sha(catalog/'receipt.json')
            or audit['plan_sha256']!=sha(plan['readback_plan'])):raise ValueError('Catalog readback incomplete')
    for name,digest in receipt['artifacts'].items():
        if sha(catalog/name)!=digest:raise ValueError('Changed catalog artifact')
    bridge_receipt=json.loads(Path(plan['bridge_receipt']).read_text())
    bridge_audit=json.loads(Path(plan['bridge_readback']).read_text())
    if bridge_audit['status']!=plan['bridge_readback_status']:raise ValueError('Family bridge not audited')
    if bridge_audit['producer_receipt_sha256']!=sha(plan['bridge_receipt']):raise ValueError('Wrong family bridge readback')
    if sha(plan['bridge'])!=bridge_receipt['bridge_sha256']:raise ValueError('Changed family bridge')
    out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
    out.mkdir();db=out/'structure_family_bridge.sqlite'
    c=sqlite3.connect(db,uri=True)
    c.execute('PRAGMA cache_size=-262144');c.execute('PRAGMA temp_store=MEMORY')
    c.execute('ATTACH DATABASE ? AS bridge',('file:'+str(Path(plan['bridge']).resolve())+'?mode=ro',))
    c.execute('CREATE TABLE raw_links(taxon_id TEXT,protein_id TEXT,sequence_sha256 TEXT,model_id TEXT,version INTEGER,model_path TEXT,PRIMARY KEY(taxon_id,protein_id))')
    with (catalog/'protein_model_links.tsv').open() as handle:
        rows=csv.DictReader(handle,delimiter='\t')
        c.executemany('INSERT INTO raw_links VALUES(?,?,?,?,?,?)',((r['taxon_id'],r['protein_id'],r['sequence_sha256'],r['model_id'],int(r['version']),r['model_path']) for r in rows))
    n=c.execute('SELECT COUNT(*) FROM raw_links').fetchone()[0]
    if n!=audit['protein_links']:raise ValueError('Catalog link count differs')
    c.execute('''CREATE TABLE structures(native_gene_id TEXT PRIMARY KEY,taxon_id TEXT,protein_id TEXT,
                 sequence_sha256 TEXT,model_id TEXT,version INTEGER,model_path TEXT,UNIQUE(taxon_id,protein_id))''')
    c.execute('''INSERT INTO structures SELECT p.native_gene_id,l.* FROM raw_links l
                 JOIN bridge.proteins p ON p.taxon_id=l.taxon_id AND p.protein_id=l.protein_id
                 WHERE p.sequence_id='S'||l.sequence_sha256''')
    if c.execute('SELECT COUNT(*) FROM structures').fetchone()[0]!=n:raise ValueError('Structure links do not match family protein sequences')
    c.execute('DROP TABLE raw_links');c.commit()
    summaries=[]
    fields=['guide','family','proteins','taxa','structure_proteins','structure_taxa','structure_sequences','structure_models']
    with (out/'family_structure_coverage.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n');writer.writerow(fields)
        for expected in bridge_receipt['guides']:
            guide=expected['guide'];families=proteins=linked=multi_taxon=fully_covered=0
            for row in family_rows(c,guide):
                family,count,taxa,covered,covered_taxa,sequences,models=row
                if not 0<=covered<=count or not 0<=covered_taxa<=taxa or not models<=sequences<=covered:
                    raise ValueError('Impossible family coverage')
                writer.writerow([guide,*row]);families+=1;proteins+=count;linked+=covered
                multi_taxon+=covered_taxa>=2;fully_covered+=covered==count
            if (families,proteins,linked)!=(expected['families'],expected['proteins'],n):
                raise ValueError('Incomplete family partition or coverage')
            summaries.append(dict(guide=guide,families=families,proteins=proteins,structure_proteins=linked,
                                  families_with_models_in_multiple_taxa=multi_taxon,families_with_all_proteins_modeled=fully_covered))
            print(summaries[-1],flush=True)
    if c.execute('PRAGMA quick_check').fetchall()!=[('ok',)]:raise ValueError('Output database integrity failed')
    c.close();verify()
    result=dict(status='complete_structure_family_coverage_bridge_pending_independent_readback',
                protein_links=n,guides=summaries,plan_sha256=plan_sha,catalog_receipt_sha256=sha(catalog/'receipt.json'),
                catalog_readback_sha256=sha(plan['readback']),artifacts={p.name:sha(p) for p in out.iterdir()},
                scope='Exact taxon/protein/sequence join to both full family partitions. Includes families without structures. Unqualified model availability; family membership is not reconciled orthology. No duplication, domain-event, homology or structural acceleration inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
