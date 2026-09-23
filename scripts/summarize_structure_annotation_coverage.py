#!/usr/bin/env python3
"""Measure full-catalog structural availability by taxon and Pfam hit status."""
import argparse
from collections import Counter,defaultdict
import csv
import json
from pathlib import Path
import sqlite3
from catalog_whole_proteome_structures import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Plan changed')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Source changed: '+p)
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);modeled={}
    with Path(plan['links']).open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            key=(r['taxon_id'],r['protein_id'])
            if key in modeled:raise ValueError('Duplicate catalog protein')
            modeled[key]='S'+r['sequence_sha256']
    catalog_count=len(modeled);counts=defaultdict(Counter)
    con=sqlite3.connect('file:'+str(Path(plan['architectures']).resolve())+'?mode=ro',uri=True)
    n=0
    for taxon,protein,sequence,raw in con.execute('SELECT p.taxon_id,p.protein_id,p.sequence_id,q.raw_hits FROM proteins p JOIN queries q USING(sequence_id)'):
        n+=1;selected=modeled.pop((taxon,protein),None)
        if selected is not None and selected!=sequence:raise ValueError('Catalog/annotation sequence mismatch')
        key=(taxon,'with_pfam_hit' if raw else 'without_pfam_hit');counts[key]['proteins']+=1;counts[key]['modeled']+=selected is not None
    if modeled or n!=plan['expected_proteins'] or catalog_count!=plan['expected_links']:raise ValueError('Incomplete scope')
    # Independent SQL grouping of the source total; separate verified structure
    # bridge supplies model presence, rather than reusing the catalog dictionary.
    con.execute('ATTACH DATABASE ? AS atlas',('file:'+str(Path(plan['structure_bridge']).resolve())+'?mode=ro',))
    sql={}
    for taxon,hit,total,models in con.execute('SELECT p.taxon_id,(q.raw_hits>0),count(*),sum(s.native_gene_id IS NOT NULL) FROM proteins p JOIN queries q USING(sequence_id) LEFT JOIN atlas.structures s ON p.taxon_id=s.taxon_id AND p.protein_id=s.protein_id GROUP BY p.taxon_id,(q.raw_hits>0)'):
        sql[(taxon,'with_pfam_hit' if hit else 'without_pfam_hit')]={'proteins':total,'modeled':models}
    if dict(counts)!=sql:raise ValueError('Independent SQL/source-bridge aggregation differs')
    rows=[];pooled=defaultdict(Counter)
    for (taxon,hit),c in sorted(counts.items()):
        rows.append(dict(taxon_id=taxon,annotation_status=hit,proteins=c['proteins'],modeled=c['modeled'],without_model=c['proteins']-c['modeled'],modeled_fraction=c['modeled']/c['proteins']))
        pooled[hit].update(c)
    taxon_counts=defaultdict(Counter)
    for (taxon,_),c in counts.items():taxon_counts[taxon].update(c)
    with Path(plan['taxon_coverage']).open() as f:
        coverage=list(csv.DictReader(f,delimiter='\t'))
    if len(coverage)!=len(taxon_counts):raise ValueError('Taxon grid differs')
    for r in coverage:
        c=taxon_counts[r['taxon_id']]
        if c['proteins']!=int(r['representative_proteins']) or c['modeled']!=int(r['proteins_with_model']):raise ValueError('Original catalog taxon counts differ')
    table=out/'taxon_annotation_coverage.tsv'
    with table.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    con.close();verify()
    receipt={'status':'complete_full_structure_annotation_coverage_with_independent_aggregation','proteins':n,'modeled_protein_links':catalog_count,'taxa':len(taxon_counts),'strata':len(rows),'pooled':{k:dict(v,modeled_fraction=v['modeled']/v['proteins']) for k,v in pooled.items()},'plan_sha256':ph,'artifacts':{'taxon_annotation_coverage.tsv':sha(table)},'scope':'Descriptive frozen AlphaFold availability by taxon and any raw Pfam gathering-threshold hit. Complete protein-level join and independent SQL/structure-bridge aggregation; original taxon totals checked. Models not confidence-qualified; proteins and taxa not independent replicates. No causal, novelty, loss or statistical significance claim; excludes ESMFold additions.'}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)


if __name__=='__main__':main()
