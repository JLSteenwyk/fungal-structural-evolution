#!/usr/bin/env python3
"""Export all screened whole-protein model partitions and descriptive weighting summaries."""
import csv,itertools,json,subprocess
from pathlib import Path
import duckdb
from screen_duplication_alignment_reuse import sha

METRICS=['rmsd_difference','mean_endpoint_tm_divergence_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference','target_rmsd_recomputed','background_rmsd_recomputed','target_sequence_distance','background_sequence_distance','sequence_distance_difference','log_positive_sequence_distance_difference']
def main():
 bindings={}
 def js(path):path=Path(path);bindings[str(path)]=sha(path);return json.loads(path.read_text())
 def bindroot(root):
  r=js(root/'receipt.json')
  for n,h in r['artifacts'].items():assert sha(root/n)==h;bindings[str(root/n)]=h
  return r
 source=Path('results/structural_comparisons/matched-whole-protein-contrasts-20260928-v1');r=bindroot(source);proof=js('metadata/matched_whole_protein_contrasts_completed_readback_20260928.json');assert proof['status']=='passed_full_matched_whole_protein_contrast_readback' and proof['source_receipt_sha256']==sha(source/'receipt.json')
 unit='fungal-matched-whole-protein-contrast-readback-20260928.service';state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 selection=Path('results/orthology/background-control-selection-20260927-v1');sr=bindroot(selection);attr=Path('results/structural_comparisons/fixed-match-structural-attrition-20260928-v1');ar=bindroot(attr);ap=js('metadata/fixed_match_structural_attrition_completed_readback_20260928.json');assert ap['source_receipt_sha256']==sha(attr/'receipt.json')
 expected={tuple(x[k] for k in ['guide','policy','scenario_id','mask','screen']):x for x in csv.DictReader((attr/'attrition_summary.tsv').open(),delimiter='\t')}
 out=Path('results/structural_comparisons/matched-whole-protein-records-20260928-v1');out.mkdir(exist_ok=False);(out/'pair_values').mkdir()
 db=duckdb.connect(str(out/'summaries.duckdb'));db.execute('SET threads=1');db.execute("SET memory_limit='12GB'");db.execute("SET max_temp_directory_size='20GB'")
 db.execute("CREATE TABLE selections AS SELECT target_id,background_id,policy,scenario_id FROM read_csv(?,delim='\t',header=true,all_varchar=true)",[str(selection/'selections.tsv.gz')]);assert db.execute('SELECT count(*) FROM selections').fetchone()[0]==2786912
 fixed=['target_id','background_id','guide','family','focal_taxon','species_pattern_id','mask','target_order','background_order','contrast_status','sequence_distance_status']
 flags=[s['id']+'_'+suffix for s in r['screens'] for suffix in ['mask_pass','both_masks_pass']]
 db.execute('CREATE TABLE measurements AS SELECT '+','.join(fixed+flags+["CAST(NULLIF("+m+",'') AS DOUBLE) AS "+m for m in METRICS])+" FROM read_csv(?,delim='\t',header=true,all_varchar=true)",[str(source/'whole_protein_contrasts.tsv.gz')]);assert db.execute('SELECT count(*) FROM measurements').fetchone()[0]==421400
 groups=sorted({k[:3] for k in expected});assert len(groups)==432;db.execute('CREATE TABLE strata(guide VARCHAR,policy VARCHAR,scenario_id VARCHAR)');db.executemany('INSERT INTO strata VALUES (?,?,?)',groups)
 keys='guide,policy,scenario_id';settings=list(itertools.product(['full','plddt70'],['same_mask','both_masks'],[s['id'] for s in r['screens']],['0','1'],['0','1']));assert len(settings)==96
 fields=['guide','policy','scenario_id','mask','cohort','screen','target_order','background_order','metadata_matched','retained_records','families','taxa','backgrounds','positive_sequence_pairs','positive_sequence_families','positive_sequence_taxa']+[m+'_'+suffix for m in METRICS for suffix in ['record_mean','family_equal_mean','taxon_equal_mean']]
 parts=[];total=0
 with (out/'record_summary.tsv').open('w') as f:
  writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
  for ix,(mask,cohort,screen,to,bo) in enumerate(settings):
   flag=screen+('_mask_pass' if cohort=='same_mask' else '_both_masks_pass')
   db.execute('CREATE OR REPLACE TABLE pair_values AS SELECT '+','.join(fixed[:6]+['sequence_distance_status']+METRICS)+' FROM measurements WHERE mask=? AND target_order=? AND background_order=? AND '+flag+"='1'",[mask,to,bo]);n=db.execute('SELECT count(*) FROM pair_values').fetchone()[0]
   assert db.execute('SELECT count(*) FROM (SELECT target_id,background_id,count(*) n FROM pair_values GROUP BY target_id,background_id HAVING n<>1)').fetchone()[0]==0
   parquet=out/'pair_values'/f'{ix:03d}.parquet';db.execute('COPY pair_values TO ? (FORMAT PARQUET)',[str(parquet)])
   db.execute('CREATE OR REPLACE VIEW records AS SELECT s.*,p.* EXCLUDE(target_id,background_id) FROM selections s JOIN pair_values p USING(target_id,background_id)')
   db.execute('CREATE OR REPLACE TABLE record_stats AS SELECT '+keys+',count(*) retained_records,count(DISTINCT family) families,count(DISTINCT focal_taxon) taxa,count(DISTINCT background_id) backgrounds,count(log_positive_sequence_distance_difference) positive_sequence_pairs,count(DISTINCT family) FILTER (WHERE log_positive_sequence_distance_difference IS NOT NULL) positive_sequence_families,count(DISTINCT focal_taxon) FILTER (WHERE log_positive_sequence_distance_difference IS NOT NULL) positive_sequence_taxa,'+','.join('avg('+m+') '+m+'_record_mean' for m in METRICS)+' FROM records GROUP BY '+keys)
   for field,label in [('family','family'),('focal_taxon','taxon')]:
    db.execute('CREATE OR REPLACE TABLE '+label+'_stats AS SELECT '+keys+','+','.join('avg('+m+') '+m+'_'+label+'_equal_mean' for m in METRICS)+' FROM (SELECT '+keys+','+field+','+','.join('avg('+m+') '+m for m in METRICS)+' FROM records GROUP BY '+keys+','+field+') GROUP BY '+keys)
   cursor=db.execute('SELECT s.*,r.* EXCLUDE(guide,policy,scenario_id),f.* EXCLUDE(guide,policy,scenario_id),t.* EXCLUDE(guide,policy,scenario_id) FROM strata s LEFT JOIN record_stats r USING('+keys+') LEFT JOIN family_stats f USING('+keys+') LEFT JOIN taxon_stats t USING('+keys+') ORDER BY s.guide,s.policy,s.scenario_id');names=[x[0] for x in cursor.description];rows=cursor.fetchall();assert len(rows)==432
   for values in rows:
    row=dict(zip(names,values));e=expected[row['guide'],row['policy'],row['scenario_id'],mask if cohort=='same_mask' else 'both_masks',screen]
    for field in ['retained_records','families','taxa','backgrounds','positive_sequence_pairs','positive_sequence_families','positive_sequence_taxa']:row[field]=row[field] or 0
    assert row['retained_records']==int(e['both_pass']) and row['positive_sequence_pairs']<=row['retained_records']
    row.update(mask=mask,cohort=cohort,screen=screen,target_order=to,background_order=bo,metadata_matched=int(e['metadata_matched']));writer.writerow(row);total+=1
   parts.append(dict(index=ix,mask=mask,cohort=cohort,screen=screen,target_order=to,background_order=bo,pairs=n,path=str(parquet.relative_to(out)),sha256=sha(parquet)));print('Prepared whole protein setting',ix+1,'/96',flush=True)
 assert total==41472;db.close()
 for p,h in bindings.items():assert sha(p)==h,p
 (out/'partition_manifest.json').write_text(json.dumps(parts,indent=2)+'\n')
 receipt=dict(status='complete_matched_whole_protein_record_summaries_pending_independent_readback',script_sha256=sha(__file__),source_hashes=bindings,settings=len(parts),summary_rows=total,metrics=METRICS,duckdb_version=duckdb.__version__,artifacts={n:sha(out/n) for n in ['record_summary.tsv','partition_manifest.json']},scope='All96 mask/cohort/screen/order settings and432 matching strata incl empty. Record/family-equal/taxon-equal descriptive means, not effect inference; log-positive summaries use only both-positive records and report their count. Original selection reuse not pooled across scenarios. Pair parquet partitions preserve model inputs and species-pattern IDs; no phylogenetic fit, uncertainty interval or significance claim.')
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('Complete whole protein record summaries',flush=True)
if __name__=='__main__':main()
