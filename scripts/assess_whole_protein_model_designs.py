#!/usr/bin/env python3
"""Assess all whole-protein sequence/structure working designs before mixed-model fitting."""
import json,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from scipy.linalg import qr
from screen_duplication_alignment_reuse import sha
KEYS=['guide','policy','scenario_id'];SETTING=['mask','cohort','screen','target_order','background_order']
NUISANCE=['original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']
VARIANTS=['gene_distance_linear','gene_distance_quadratic','gene_distance_cubic','positive_log_gene_distance_linear','alignment_identity_linear']
def assess(x,names):
 assert x.ndim==2 and x.shape[1]==len(names) and np.isfinite(x).all()
 if not len(x):return dict(records=0,status='no_records',columns=names,active_columns=[],constant_columns=[],nonzero_constant_columns=[],centered_rank=0,intercept_design_rank=0,residual_df=0,scaled_singular_values=[],scaled_condition=None,zero_in_each_marginal_range=None,minimum=[],maximum=[],mean=[])
 lo=x.min(axis=0);hi=x.max(axis=0);active=hi-lo>1e-12;center=x[:,active]-x[:,active].mean(axis=0);scale=np.std(x[:,active],axis=0)
 z=center/scale/np.sqrt(len(x));s=np.linalg.svd(z,compute_uv=False);rank=int((s>1e-7).sum())
 if active.any():
  _,tri,_=qr(z,mode='economic',pivoting=True);s2=np.linalg.svd(tri,compute_uv=False);np.testing.assert_allclose(s,s2,rtol=1e-9,atol=1e-12);assert rank==int((s2>1e-7).sum())
 nonzero=[names[i] for i in range(len(names)) if not active[i] and abs(x[:,i].mean())>1e-12];df=len(x)-rank-1
 status='rank_deficient' if rank<int(active.sum()) else 'no_residual_df' if df<=0 else 'nonzero_constant_changes_intercept' if nonzero else 'full_rank_with_positive_residual_df'
 return dict(records=len(x),status=status,columns=names,active_columns=[n for n,b in zip(names,active) if b],constant_columns=[n for n,b in zip(names,active) if not b],nonzero_constant_columns=nonzero,centered_rank=rank,intercept_design_rank=rank+1,residual_df=df,scaled_singular_values=s.tolist(),scaled_condition=float(s[0]/s[-1]) if len(s) and rank==int(active.sum()) else None,zero_in_each_marginal_range=bool(np.all((lo<=0)&(hi>=0))),minimum=lo.tolist(),maximum=hi.tolist(),mean=x.mean(axis=0).tolist())
def matrix(frame,variant):
 if variant=='positive_log_gene_distance_linear':frame=frame.loc[(frame.target_sequence_distance>0)&(frame.background_sequence_distance>0)]
 t=frame.target_sequence_distance.to_numpy();b=frame.background_sequence_distance.to_numpy()
 if variant.startswith('gene_distance_'):
  degree={'gene_distance_linear':1,'gene_distance_quadratic':2,'gene_distance_cubic':3}[variant];terms=[t**k-b**k for k in range(1,degree+1)];names=['gene_distance_power_'+str(k)+'_difference' for k in range(1,degree+1)]
 elif variant=='positive_log_gene_distance_linear':terms=[np.log(t)-np.log(b)];names=['log_gene_distance_difference']
 else:assert variant=='alignment_identity_linear';terms=[frame.identity_difference.to_numpy()];names=['alignment_identity_difference']
 return frame,np.column_stack(terms+[frame[k].to_numpy() for k in NUISANCE]),names+NUISANCE

def main():
 lp=Path('metadata/matched_whole_protein_record_readback_launch_20260928.json');launch=json.loads(lp.read_text());lh=sha(lp)
 while psutil.pid_exists(launch['pid']):
  try:
   p=psutil.Process(launch['pid'])
   if abs(p.create_time()-launch['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
   assert p.cmdline()==launch['cmdline']
  except psutil.NoSuchProcess:break
  time.sleep(30)
 state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 root=Path('results/structural_comparisons/matched-whole-protein-records-20260928-v1');rp=root/'receipt.json';r=json.loads(rp.read_text());ap=Path('metadata/matched_whole_protein_records_completed_readback_20260928.json');a=json.loads(ap.read_text());assert a['status']=='passed_full_matched_whole_protein_record_readback' and a['source_receipt_sha256']==sha(rp)
 bindings={str(rp):sha(rp),str(ap):sha(ap),str(lp):lh}
 for name,h in r['artifacts'].items():assert sha(root/name)==h;bindings[str(root/name)]=h
 sp=Path('results/orthology/background-control-selection-20260927-v1/selections.tsv.gz');assert sha(sp)==r['source_hashes'][str(sp)];bindings[str(sp)]=sha(sp)
 selected=pd.read_csv(sp,sep='\t',usecols=['target_id','background_id','policy','scenario_id']);summary=pd.read_csv(root/'record_summary.tsv',sep='\t',dtype={'target_order':str,'background_order':str});parts=json.loads((root/'partition_manifest.json').read_text())
 out=Path('results/structural_comparisons/whole-protein-model-designs-20260928-v1');out.mkdir(exist_ok=False);counts={};total=0
 with (out/'designs.jsonl').open('w') as f:
  for part in parts:
   path=root/part['path'];assert sha(path)==part['sha256'];bindings[str(path)]=sha(path);values=pd.read_parquet(path);records=selected.merge(values,on=['target_id','background_id'],validate='many_to_one');groups=records.groupby(KEYS).indices
   subset=summary
   for k in SETTING:subset=subset.loc[subset[k].eq(part[k])]
   assert len(subset)==432
   for row in subset.to_dict('records'):
    key=tuple(row[k] for k in KEYS);frame=records.iloc[groups.get(key,[])];assert len(frame)==row['retained_records']
    for variant in VARIANTS:
     chosen,x,names=matrix(frame,variant);result=dict(**{k:row[k] for k in KEYS+SETTING},variant=variant,screen_retained_records=len(frame),sequence_excluded_records=len(frame)-len(chosen),**assess(x,names))
     if variant=='positive_log_gene_distance_linear':assert len(chosen)==row['positive_sequence_pairs']
     result.update(targets=int(chosen.target_id.nunique()),backgrounds=int(chosen.background_id.nunique()),families=int(chosen.family.nunique()),taxa=int(chosen.focal_taxon.nunique()),species_patterns=int(chosen.species_pattern_id.nunique()),maximum_background_reuse=int(chosen.background_id.value_counts().max()) if len(chosen) else 0)
     f.write(json.dumps(result,separators=(',',':'),allow_nan=False)+'\n');total+=1;counts[result['status']]=counts.get(result['status'],0)+1
   print('Assessed whole protein setting',part['index']+1,'/96',flush=True)
 assert total==96*432*5==207360
 for path,h in bindings.items():assert sha(path)==h,path
 receipt=dict(status='complete_whole_protein_model_design_assessment_pending_independent_readback',script_sha256=sha(__file__),source_hashes=bindings,designs=total,variants=VARIANTS,counts=counts,artifacts={'designs.jsonl':sha(out/'designs.jsonl')},scope='Every screened stratum and five prespecified working designs, including empty and zero-distance exclusions. SVD rank cross-checked by pivoted QR; constants, nonzero constants, residual df and marginal zero overlap explicit. Marginal overlap is not joint support; numerical rank is not model adequacy, calibrated uncertainty or phylogenetic effect inference. Sequence distances are relative gene-tree distances, not dated rates.')
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='source_hashes'}),flush=True)
if __name__=='__main__':main()
