#!/usr/bin/env python3
"""Read back the full interval manifest against an independent SQL union."""
import csv
import hashlib
import json
from pathlib import Path
import sqlite3
from catalog_whole_proteome_structures import sha


def main():
    root=Path('results/domains/domain-extraction-manifest-20260923-v1');rp=root/'receipt.json';r=json.loads(rp.read_text())
    db=Path('results/domains/whole-proteome-structure-domain-registry-20260923-v1/structure_domains.sqlite')
    if sha(db)!=r['registry_database_sha256']:raise ValueError('Changed source registry')
    for name,h in r['artifacts'].items():
        if sha(root/name)!=h:raise ValueError('Changed manifest output')
    con=sqlite3.connect('file:'+str(db.resolve())+'?mode=ro',uri=True)
    expected=set(con.execute('''WITH eligible AS (SELECT DISTINCT model_key,hit_id FROM policy_hits WHERE domain_interval_candidate=1), spans AS (
SELECT s.model_key,alignment_start AS start,alignment_end AS end FROM segments s JOIN eligible e USING(model_key,hit_id)
UNION SELECT s.model_key,envelope_start,envelope_end FROM segments s JOIN eligible e USING(model_key,hit_id))
SELECT spans.model_key,m.sequence_sha256,m.path,start,end,m.length FROM spans JOIN models m USING(model_key)'''))
    models=set();count=0;residues=0;maximum=0;ids=set()
    with (root/'intervals.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            start,end,n=map(int,(row['start'],row['end'],row['residues']))
            iid='D'+hashlib.sha256(f"{row['full_sequence_sha256']}:{start}:{end}".encode()).hexdigest()
            if iid!=row['interval_id'] or iid in ids or n!=end-start+1:raise ValueError('Interval identity or length differs')
            ids.add(iid);models.add(row['model_key'])
            length=con.execute('SELECT length FROM models WHERE model_key=?',(row['model_key'],)).fetchone()[0]
            identity=(row['model_key'],row['full_sequence_sha256'],row['model_path'],start,end,length)
            if identity not in expected or not 1<=start<=end<=length:raise ValueError('Interval absent, repeated or out of bounds')
            expected.remove(identity);count+=1;residues+=n;maximum=max(maximum,n)
    if expected or count!=r['unique_intervals'] or len(models)!=r['source_models'] or residues!=r['unique_interval_residues'] or maximum!=r['maximum_interval_residues']:raise ValueError('Scope or dimensions differ')
    con.close()
    result={'status':'passed_full_domain_interval_manifest_sql_union_readback','intervals':count,'models':len(models),'residues':residues,'maximum_interval_residues':maximum,'source_receipt_sha256':sha(rp),'script_sha256':sha(__file__),'scope':'Independent SQL union of alignment/envelope spans compared with every emitted interval, exact source model/sequence/path, inclusive bounds and interval hash. Complements producer boundary-link round trip; does not extract or validate coordinates.'}
    target=Path('metadata/domain_extraction_manifest_readback.json')
    if target.exists():raise FileExistsError(target)
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
