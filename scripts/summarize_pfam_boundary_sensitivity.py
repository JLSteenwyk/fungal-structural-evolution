#!/usr/bin/env python3
"""Summarize full domain-boundary sensitivity by Pfam, with original-registry SQL checks."""
import argparse
import json
from pathlib import Path
import sqlite3
import pandas as pd
from catalog_whole_proteome_structures import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registry',type=Path,required=True)
    p.add_argument('--registry-readback',type=Path,required=True)
    p.add_argument('--boundaries',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    rr=json.loads((a.registry/'receipt.json').read_text()); audit=json.loads(a.registry_readback.read_text())
    db=a.registry/'structure_domains.sqlite'; database_hash=sha(db)
    if audit['status']!='passed_full_structure_domain_registry_readback' or audit['database_sha256']!=database_hash or rr['database_sha256']!=database_hash or audit['producer_receipt_sha256']!=sha(a.registry/'receipt.json'):
        raise ValueError('Registry audit binding differs')
    br=json.loads((a.boundaries/'receipt.json').read_text()); ba=json.loads((a.boundaries/'readback.json').read_text())
    if ba['status']!='passed_all_boundary_pair_rows_and_threshold_counts' or ba['producer_receipt_sha256']!=sha(a.boundaries/'receipt.json'):
        raise ValueError('Boundary audit binding differs')
    table=a.boundaries/'boundary_sensitivity.tsv.gz'
    if sha(table)!=br['artifacts'][table.name]:raise ValueError('Changed boundary table')
    c=sqlite3.connect('file:'+str(db.resolve())+'?mode=ro',uri=True)
    candidates='SELECT DISTINCT model_key,hit_id FROM policy_hits WHERE domain_interval_candidate=1'
    annotations=pd.read_sql_query('SELECT s.* FROM segments s JOIN ('+candidates+') q USING(model_key,hit_id)',c)
    b=pd.read_csv(table,sep='\t')
    joined=b.merge(annotations,on=['model_key','hit_id'],how='outer',validate='one_to_one',suffixes=('','_registry'),indicator=True)
    if not joined['_merge'].eq('both').all() or not joined.pfam_type.eq('Domain').all():raise ValueError('Candidate identity grid differs')
    for k in ['alignment_start','alignment_end','envelope_start','envelope_end']:
        if not joined[k].eq(joined[k+'_registry']).all():raise ValueError('Boundary endpoint mismatch')
    joined['identical']=joined.added_residues.eq(0).astype(int)
    joined['added_ge10']=joined.added_residues.ge(10).astype(int)
    joined['added_ge20pct']=(5*joined.added_residues>=joined.alignment_residues).astype(int)
    summary=joined.groupby('pfam_accession').agg(candidate_pairs=('hit_id','size'),models=('model_key','nunique'),identical_pairs=('identical','sum'),added_ge10=('added_ge10','sum'),added_ge20pct=('added_ge20pct','sum'),total_added_residues=('added_residues','sum'),maximum_added_residues=('added_residues','max'))
    # Separate aggregation directly from original endpoints, without using the joined result or its flags.
    query='''SELECT s.pfam_accession, COUNT(*) candidate_pairs, COUNT(DISTINCT s.model_key) models,
    SUM(s.alignment_start=s.envelope_start AND s.alignment_end=s.envelope_end) identical_pairs,
    SUM((s.envelope_end-s.envelope_start)-(s.alignment_end-s.alignment_start)>=10) added_ge10,
    SUM(5*((s.envelope_end-s.envelope_start)-(s.alignment_end-s.alignment_start))>=s.alignment_end-s.alignment_start+1) added_ge20pct,
    SUM((s.envelope_end-s.envelope_start)-(s.alignment_end-s.alignment_start)) total_added_residues,
    MAX((s.envelope_end-s.envelope_start)-(s.alignment_end-s.alignment_start)) maximum_added_residues
    FROM segments s JOIN ('''+candidates+''') q USING(model_key,hit_id) GROUP BY s.pfam_accession'''
    independent=pd.read_sql_query(query,c).set_index('pfam_accession').sort_index()
    pd.testing.assert_frame_equal(summary.sort_index(),independent,check_dtype=False)
    if int(summary.candidate_pairs.sum())!=br['candidate_model_hit_pairs'] or int(summary.added_ge10.sum())!=br['pairs_with_at_least_10_added_residues'] or int(summary.added_ge20pct.sum())!=br['pairs_with_at_least_20_percent_added_residues']:
        raise ValueError('Global summary differs')
    summary['fraction_added_ge10']=summary.added_ge10/summary.candidate_pairs
    summary['fraction_added_ge20pct']=summary.added_ge20pct/summary.candidate_pairs
    summary=summary.sort_values(['added_ge20pct','candidate_pairs'],ascending=False).reset_index()
    a.output.mkdir(parents=True,exist_ok=False);target=a.output/'pfam_boundary_sensitivity.tsv';summary.to_csv(target,sep='\t',index=False)
    if sha(db)!=database_hash or sha(table)!=br['artifacts'][table.name]:raise ValueError('Sources changed during analysis')
    r={'status':'complete_pfam_boundary_sensitivity_with_full_sql_aggregation_readback','candidate_pairs':len(joined),'pfam_accessions':len(summary),'pfams_with_at_least_one_20_percent_extension':int(summary.added_ge20pct.gt(0).sum()),'source_hashes':{str(db):database_hash,str(a.registry/'receipt.json'):sha(a.registry/'receipt.json'),str(a.registry_readback):sha(a.registry_readback),str(a.boundaries/'receipt.json'):sha(a.boundaries/'receipt.json'),str(a.boundaries/'readback.json'):sha(a.boundaries/'readback.json'),str(table):sha(table)},'script_sha256':sha(__file__),'artifacts':{target.name:sha(target)},'scope':'All candidate model/hit pairs joined exactly to audited Pfam accessions and original endpoints; all per-Pfam integer aggregates independently recomputed in SQL. Counts are unique model/hit observations, not species, independent transitions or confidence-qualified domains. Descriptive ranking by absolute count, without enrichment or evolutionary inference.'}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='source_hashes'},indent=2));print(summary.head(10).to_string(index=False))


if __name__=='__main__':main()
