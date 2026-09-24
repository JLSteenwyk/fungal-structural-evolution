#!/usr/bin/env python3
"""Report domain-cluster representation for every sampled taxon, including zeros."""
import argparse
import csv
import hashlib
import json
import sqlite3
from collections import defaultdict
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8388608),b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for key in ('composition','audit','manifest','output'):
        ap.add_argument('--'+key,type=Path,required=True)
    a = ap.parse_args()
    rp = a.composition/'receipt.json'
    receipt = json.loads(rp.read_text())
    audit = json.loads(a.audit.read_text())
    dbpath = a.composition/'domain_cluster_composition.sqlite'
    if audit['status'] != 'passed_full_domain_cluster_composition_source_readback' or audit['producer_receipt_sha256'] != sha(rp):
        raise ValueError('Missing bound full audit')
    if sha(dbpath) != receipt['artifacts'][dbpath.name]:
        raise ValueError('Database changed')
    pins = {str(p):sha(p) for p in (rp,a.audit,dbpath,a.manifest)}
    with a.manifest.open() as f:
        rows = list(csv.DictReader(f,delimiter='\t'))
    manifest = {r['taxon_id']:r for r in rows}
    if len(rows) != 526 or len(manifest) != 526:
        raise ValueError('Unexpected sampling universe')
    c = sqlite3.connect('file:'+str(dbpath.resolve())+'?mode=ro',uri=True)
    c.execute('PRAGMA cache_size=-131072')
    genes = {}
    models = defaultdict(set)
    proteins = defaultdict(set)
    for gene,taxon,model in c.execute('SELECT native_gene_id,taxon_id,model FROM proteins'):
        if gene in genes or taxon not in manifest:
            raise ValueError('Invalid protein identity')
        genes[gene] = (taxon,model)
        proteins[taxon].add(gene)
        models[taxon].add(model)
    model_intervals = defaultdict(set)
    for interval,model in c.execute('SELECT interval_id,model FROM intervals'):
        model_intervals[model].add(interval)
    intervals = {taxon:set().union(*(model_intervals[m] for m in ms)) for taxon,ms in models.items()}
    clusters = defaultdict(set)
    for rep,gene in c.execute('SELECT representative,native_gene_id FROM cluster_proteins'):
        clusters[genes[gene][0]].add(rep)
    families = {g:defaultdict(set) for g in ('profile','mafft')}
    for guide,gene,family in c.execute('SELECT guide,native_gene_id,family FROM family_links'):
        families[guide][genes[gene][0]].add(family)
    expected = {t:(len(models[t]),len(proteins[t]),len(intervals.get(t,set())),len(clusters[t]),len(families['profile'][t]),len(families['mafft'][t])) for t in manifest}
    # SQL independently aggregates each relation, avoiding multiplicative joins.
    sql = {t:[0]*6 for t in manifest}
    for t,n,m in c.execute('SELECT taxon_id,COUNT(*),COUNT(DISTINCT model) FROM proteins GROUP BY taxon_id'):
        sql[t][0:2] = [m,n]
    for t,n in c.execute('SELECT p.taxon_id,COUNT(DISTINCT i.interval_id) FROM proteins p JOIN intervals i USING(model) GROUP BY p.taxon_id'):
        sql[t][2] = n
    for t,n in c.execute('SELECT p.taxon_id,COUNT(DISTINCT cp.representative) FROM cluster_proteins cp JOIN proteins p USING(native_gene_id) GROUP BY p.taxon_id'):
        sql[t][3] = n
    for t,g,n in c.execute('SELECT p.taxon_id,f.guide,COUNT(DISTINCT f.family) FROM family_links f JOIN proteins p USING(native_gene_id) GROUP BY p.taxon_id,f.guide'):
        sql[t][4 if g=='profile' else 5] = n
    if any(tuple(sql[t]) != expected[t] for t in manifest):
        raise ValueError('SQL and set aggregation differ')
    if len(genes) != audit['proteins']:
        raise ValueError('Protein universe differs')
    c.close()
    a.output.mkdir(parents=True,exist_ok=False)
    table = a.output/'taxon_coverage.tsv'
    with table.open('w') as f:
        w = csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['taxon_id','species_name','study_role','lineage','models','proteins','intervals','clusters','profile_families','mafft_families'])
        for t in sorted(manifest):
            row = manifest[t]
            w.writerow([t,row['species_name'],row['study_role'],row['lineage'],*expected[t]])
    for path,digest in pins.items():
        if sha(path) != digest:
            raise ValueError('Source changed')
    result = dict(status='complete_all_taxon_domain_coverage_with_full_sql_set_agreement',taxa=526,
                  represented_taxa=sum(v[1]>0 for v in expected.values()),unrepresented_taxa=sum(v[1]==0 for v in expected.values()),
                  proteins=len(genes),source_hashes=pins,script_sha256=sha(__file__),artifacts={table.name:sha(table)},
                  scope='Full manifest retained, including zero representation. Counts describe the selected AlphaFold/Pfam domain catalog; zero does not establish gene/domain absence. Alternative boundaries and shared models are not independent observations. No confidence qualification or evolutionary-event inference.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
