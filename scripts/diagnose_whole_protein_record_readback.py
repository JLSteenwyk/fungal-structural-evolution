#!/usr/bin/env python3
"""Diagnostic full replay of a non-reproduced whole-protein summary audit mismatch."""
import csv,itertools,json,math,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from screen_duplication_alignment_reuse import sha
METRICS=['rmsd_difference','mean_endpoint_tm_divergence_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference','target_rmsd_recomputed','background_rmsd_recomputed','target_sequence_distance','background_sequence_distance','sequence_distance_difference','log_positive_sequence_distance_difference']
KEYS=['guide','policy','scenario_id'];SETTING=['mask','cohort','screen','target_order','background_order']
def main():
 launchpath=Path('metadata/matched_whole_protein_records_launch_20260928.json');launch=json.loads(launchpath.read_text());lh=sha(launchpath)
 while psutil.pid_exists(launch['pid']):
  try:
   p=psutil.Process(launch['pid'])
   if abs(p.create_time()-launch['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
   assert p.cmdline()==launch['cmdline']
  except psutil.NoSuchProcess:break
  time.sleep(30)
 state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 root=Path('results/structural_comparisons/matched-whole-protein-records-20260928-v1');r=json.loads((root/'receipt.json').read_text());rh=sha(root/'receipt.json');assert r['status']=='complete_matched_whole_protein_record_summaries_pending_independent_readback' and r['metrics']==METRICS
 bindings=dict(r['source_hashes'])
 for n,h in r['artifacts'].items():bindings[str(root/n)]=h
 for path,h in bindings.items():assert sha(path)==h,path
 parts=json.loads((root/'partition_manifest.json').read_text());screens=['n30_c50','n30_c70','n30_c90','n50_c50','n50_c70','n50_c90'];settings=set(itertools.product(['full','plddt70'],['same_mask','both_masks'],screens,['0','1'],['0','1']));assert len(parts)==96 and {tuple(p[k] for k in SETTING) for p in parts}==settings
 sr=json.loads(Path('results/orthology/background-control-selection-20260927-v1/receipt.json').read_text());groups={tuple(k.split('|')[:3]) for k in sr['counts'] if k.endswith('|targets')};assert len(groups)==432
 table={}
 for row in csv.DictReader((root/'record_summary.tsv').open(),delimiter='\t'):
  key=tuple(row[k] for k in SETTING+KEYS);assert key not in table;table[key]=row
 assert set(table)=={(*s,*g) for s in settings for g in groups} and len(table)==r['summary_rows']==41472
 columns=['target_id','background_id','guide','family','focal_taxon','species_pattern_id','sequence_distance_status']+METRICS
 flags=[s+'_'+cohort+'_pass' for s in screens for cohort in ['mask','both_masks']]
 d=pd.read_csv('results/structural_comparisons/matched-whole-protein-contrasts-20260928-v1/whole_protein_contrasts.tsv.gz',sep='\t',usecols=columns+['mask','target_order','background_order']+flags,dtype={'target_order':str,'background_order':str});assert len(d)==421400
 selected=pd.read_csv('results/orthology/background-control-selection-20260927-v1/selections.tsv.gz',sep='\t',usecols=['target_id','background_id','policy','scenario_id']);assert len(selected)==2786912
 values_checked=0;partrows=0
 for ix,part in enumerate(parts):
  setting=tuple(part[k] for k in SETTING);mask,cohort,screen,to,bo=setting;flag=screen+('_mask_pass' if cohort=='same_mask' else '_both_masks_pass')
  expected=d.loc[d['mask'].eq(mask)&d.target_order.eq(to)&d.background_order.eq(bo)&d[flag].eq(1),columns].copy();assert not expected.duplicated(['target_id','background_id']).any()
  path=root/part['path'];assert sha(path)==part['sha256'];actual=pd.read_parquet(path);assert len(actual)==len(expected)==part['pairs'] and set(actual)==set(columns)
  sort=['target_id','background_id'];pd.testing.assert_frame_equal(actual[columns].sort_values(sort).reset_index(drop=True),expected.sort_values(sort).reset_index(drop=True),check_dtype=False,check_exact=False,rtol=1e-12,atol=1e-12);partrows+=len(actual)
  records=selected.merge(expected,on=sort,validate='many_to_one');grouped=records.groupby(KEYS);counts=grouped.size();families=grouped.family.nunique();taxa=grouped.focal_taxon.nunique();backgrounds=grouped.background_id.nunique()
  positive=records.loc[records.log_positive_sequence_distance_difference.notna()].groupby(KEYS);poscounts=positive.size();posfamilies=positive.family.nunique();postaxa=positive.focal_taxon.nunique()
  means={'record_mean':grouped[METRICS].mean(),'family_equal_mean':records.groupby(KEYS+['family'])[METRICS].mean().groupby(level=KEYS).mean(),'taxon_equal_mean':records.groupby(KEYS+['focal_taxon'])[METRICS].mean().groupby(level=KEYS).mean()}
  for group in groups:
   row=table[(*setting,*group)];expectedcounts={'metadata_matched':sr['counts'].get('|'.join(group)+'|matched',0),'retained_records':counts.get(group,0),'families':families.get(group,0),'taxa':taxa.get(group,0),'backgrounds':backgrounds.get(group,0),'positive_sequence_pairs':poscounts.get(group,0),'positive_sequence_families':posfamilies.get(group,0),'positive_sequence_taxa':postaxa.get(group,0)}
   for k,v in expectedcounts.items():assert int(row[k])==int(v),(setting,group,k,row[k],v)
   for suffix,frame in means.items():
    for metric in METRICS:
     value=frame.loc[group,metric] if group in frame.index else np.nan;observed=row[metric+'_'+suffix]
     if pd.isna(value):assert observed==''
     else:
      if not (math.isfinite(value) and math.isclose(float(observed),value,rel_tol=1e-9,abs_tol=1e-10)):
       local=records
       for field,label in zip(KEYS,group):local=local.loc[local[field].eq(label)]
       if suffix=='record_mean':items=local[metric].dropna().tolist()
       else:
        label='family' if suffix=='family_equal_mean' else 'focal_taxon';items=[]
        for _,sub in local.groupby(label):
         finite=sub[metric].dropna().tolist()
         if finite:items.append(math.fsum(finite)/len(finite))
       independent=math.fsum(items)/len(items) if items else None
       detail=dict(setting=setting,group=group,metric=metric,weighting=suffix,stored=observed,pandas_group_mean=float(value),independent_fsum_mean=independent,local_records=len(local),pandas_version=pd.__version__,numpy_version=np.__version__)
       dp=Path('metadata/whole_protein_summary_mismatch_diagnostic_20260928.json');dp.write_text(json.dumps(detail,indent=2)+'\n');local.to_csv('results/structural_comparisons/whole-protein-summary-mismatch-records-20260928.tsv.gz',sep='\t',index=False);print(json.dumps(detail),flush=True)
       raise AssertionError('Full replay mismatch; diagnostic records preserved')
     values_checked+=1
  print('Verified whole protein summary setting',ix+1,'/96',flush=True)
 for path,h in bindings.items():assert sha(path)==h,path
 for part in parts:assert sha(root/part['path'])==part['sha256']
 assert sha(root/'receipt.json')==rh and sha(launchpath)==lh
 result=dict(status='passed_full_matched_whole_protein_record_diagnostic_replay',source_receipt_sha256=rh,checker_sha256=sha(__file__),producer_terminal_state=state,settings=len(parts),summary_rows=len(table),partition_rows_checked=partrows,weighted_mean_cells_checked=values_checked,scope='All96 pair partitions reconstructed from source contrasts and flags; all41472 strata rejoined to original selections; all counts, positive-log record/family/taxon denominators, three weighting schemes, means and empty cells independently recomputed with pandas rather than producer SQL. Descriptive diagnostics and model inputs only; no effect inference.')
 with Path('metadata/matched_whole_protein_records_diagnostic_replay_20260928.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result),flush=True)
if __name__=='__main__':main()
