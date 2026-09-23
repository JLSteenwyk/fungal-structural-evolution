#!/usr/bin/env python3
"""Prepare all policy-qualified domain boundaries with reversible associations."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sqlite3
from catalog_whole_proteome_structures import sha


def main():
    source=Path('results/domains/whole-proteome-structure-domain-registry-20260923-v1')
    rp=source/'receipt.json';ap=Path('results/domains/whole-proteome-structure-domain-registry-readback-20260923-v1/receipt.json')
    receipt=json.loads(rp.read_text());audit=json.loads(ap.read_text());database=source/'structure_domains.sqlite'
    if audit['status']!='passed_full_structure_domain_registry_readback' or audit['producer_receipt_sha256']!=sha(rp) or sha(database)!=audit['database_sha256']:raise ValueError('Unverified interval registry')
    out=Path('results/domains/domain-extraction-manifest-20260923-v1')
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);con=sqlite3.connect('file:'+str(database.resolve())+'?mode=ro',uri=True)
    query='''SELECT s.model_key,s.hit_id,s.alignment_start,s.alignment_end,s.envelope_start,s.envelope_end,m.sequence_sha256,m.path,m.length
FROM segments s JOIN models m USING(model_key)
WHERE EXISTS (SELECT 1 FROM policy_hits p WHERE p.model_key=s.model_key AND p.hit_id=s.hit_id AND p.domain_interval_candidate=1)
ORDER BY s.model_key,s.hit_id'''
    intervals={};models=set();hits=0;boundaries=0
    # One interval per exact full-sequence hash and inclusive coordinate span.
    # Boundary-policy/hit links remain explicit even when coordinates coincide.
    with (out/'boundary_links.tsv.gz').open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
            with io.TextIOWrapper(gz,encoding='utf-8',newline='') as f:
                w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['interval_id','model_key','hit_id','boundary'])
                for key,hit,astart,aend,estart,eend,sequence,path,length in con.execute(query):
                    models.add(key);hits+=1
                    for boundary,start,end in [('alignment',astart,aend),('envelope',estart,eend)]:
                        if not 1<=start<=end<=length:raise ValueError('Interval outside protein')
                        identity=f'{sequence}:{start}:{end}';iid='D'+hashlib.sha256(identity.encode()).hexdigest()
                        record=(iid,key,sequence,path,start,end,end-start+1)
                        if iid in intervals and intervals[iid]!=record:raise ValueError('Ambiguous interval identity')
                        intervals[iid]=record;w.writerow([iid,key,hit,boundary]);boundaries+=1
    table=out/'intervals.tsv'
    with table.open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['interval_id','model_key','full_sequence_sha256','model_path','start','end','residues']);w.writerows(intervals[k] for k in sorted(intervals))
    expected=con.execute('SELECT count(*) FROM (SELECT DISTINCT model_key,hit_id FROM policy_hits WHERE domain_interval_candidate=1)').fetchone()[0]
    if expected!=hits or boundaries!=2*hits:raise ValueError('Candidate/boundary scope differs')
    # Read back every emitted association and its exact source bounds.
    emitted={}
    with gzip.open(out/'boundary_links.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            label=(r['model_key'],r['hit_id'],r['boundary'])
            if label in emitted or r['interval_id'] not in intervals:raise ValueError('Duplicate or unknown interval association')
            emitted[label]=r['interval_id']
    for key,hit,astart,aend,estart,eend,sequence,path,length in con.execute(query):
        for boundary,start,end in [('alignment',astart,aend),('envelope',estart,eend)]:
            iid=emitted.pop((key,hit,boundary));record=intervals[iid]
            if record[1:]!=(key,sequence,path,start,end,end-start+1):raise ValueError('Boundary readback differs')
    if emitted:raise ValueError('Unexpected emitted association')
    con.close()
    result={'status':'complete_all_candidate_domain_extraction_manifest','source_models':len(models),'candidate_model_hit_pairs':hits,'boundary_links':boundaries,'unique_intervals':len(intervals),'unique_interval_residues':sum(r[-1] for r in intervals.values()),'maximum_interval_residues':max(r[-1] for r in intervals.values()),'registry_receipt_sha256':sha(rp),'registry_readback_sha256':sha(ap),'registry_database_sha256':audit['database_sha256'],'script_sha256':sha(__file__),'artifacts':{p.name:sha(p) for p in out.iterdir()},'scope':'Union of all four policy-qualified Domain hits, retaining alignment and envelope boundaries and exact model/hit associations. Identical source-sequence intervals deduplicated reversibly. Full boundary-link round trip checked. Coordinate extraction, residue confidence, PAE and independent structural validation not performed.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
