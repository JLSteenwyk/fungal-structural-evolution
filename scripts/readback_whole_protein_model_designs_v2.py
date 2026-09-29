#!/usr/bin/env python3
"""Independently reconstruct every sequence-aware whole-protein model design."""
import itertools,json,math,subprocess,time
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from scipy.linalg import svdvals
from screen_duplication_alignment_reuse import sha
KEYS=['guide','policy','scenario_id'];SETTING=['mask','cohort','screen','target_order','background_order']
NUISANCE=['original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']
VARIANTS=['gene_distance_linear','gene_distance_quadratic','gene_distance_cubic','positive_log_gene_distance_linear','alignment_identity_linear']
def reconstruct(frame,variant):
 if variant=='positive_log_gene_distance_linear':frame=frame.loc[(frame.target_sequence_distance>0)&(frame.background_sequence_distance>0)]
 t=frame.target_sequence_distance.to_numpy();b=frame.background_sequence_distance.to_numpy();delta=t-b
 if variant.startswith('gene_distance_'):
  degree=VARIANTS.index(variant)+1;terms=[delta,delta*(t+b),delta*(t*t+t*b+b*b)][:degree];names=['gene_distance_power_'+str(k)+'_difference' for k in range(1,degree+1)]
 elif variant.startswith('positive_log'):terms=[np.log(t/b)];names=['log_gene_distance_difference']
 else:assert variant=='alignment_identity_linear';terms=[frame.identity_difference.to_numpy()];names=['alignment_identity_difference']
 names+=NUISANCE;x=pd.DataFrame(np.column_stack(terms+[frame[k].to_numpy() for k in NUISANCE]),columns=names);n=len(x)
 result=dict(records=n,columns=names,targets=int(frame.target_id.nunique()),backgrounds=int(frame.background_id.nunique()),families=int(frame.family.nunique()),taxa=int(frame.focal_taxon.nunique()),species_patterns=int(frame.species_pattern_id.nunique()),maximum_background_reuse=int(frame.background_id.value_counts().max()) if n else 0)
 if not n:return dict(result,status='no_records',active_columns=[],constant_columns=[],nonzero_constant_columns=[],centered_rank=0,intercept_design_rank=0,residual_df=0,scaled_singular_values=[],scaled_condition=None,zero_in_each_marginal_range=None,minimum=[],maximum=[],mean=[])
 assert np.isfinite(x.to_numpy()).all();lo=x.min();hi=x.max();avg=x.mean();active=[k for k in names if hi[k]-lo[k]>1e-12];constant=[k for k in names if k not in active];nonzero=[k for k in constant if abs(avg[k])>1e-12]
 z=(x[active]-avg[active])/x[active].std(ddof=0)/math.sqrt(n);singular=svdvals(z.to_numpy());rank=int(sum(singular>1e-7));df=n-rank-1
 status='rank_deficient' if rank<len(active) else 'no_residual_df' if df<=0 else 'nonzero_constant_changes_intercept' if nonzero else 'full_rank_with_positive_residual_df'
 return dict(result,status=status,active_columns=active,constant_columns=constant,nonzero_constant_columns=nonzero,centered_rank=rank,intercept_design_rank=rank+1,residual_df=df,scaled_singular_values=singular.tolist(),scaled_condition=float(singular[0]/singular[-1]) if len(singular) and rank==len(active) else None,zero_in_each_marginal_range=bool(((lo<=0)&(hi>=0)).all()),minimum=lo.tolist(),maximum=hi.tolist(),mean=avg.tolist())
def compare(actual,expected):
 assert set(actual)==set(expected)
 for k,v in expected.items():
  if isinstance(v,list) and v and isinstance(v[0],float):np.testing.assert_allclose(actual[k],v,rtol=1e-8,atol=1e-10,err_msg=k)
  elif isinstance(v,float):assert math.isclose(actual[k],v,rel_tol=1e-6 if k=='scaled_condition' else 1e-9,abs_tol=1e-10),(k,actual[k],v)
  else:assert actual[k]==v,(k,actual[k],v)
def main():
 lp=Path('metadata/whole_protein_model_design_v2_launch_20260928.json');launch=json.loads(lp.read_text());lh=sha(lp)
 while psutil.pid_exists(launch['pid']):
  try:
   p=psutil.Process(launch['pid'])
   if abs(p.create_time()-launch['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
   assert p.cmdline()==launch['cmdline']
  except psutil.NoSuchProcess:break
  time.sleep(30)
 state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 root=Path('results/structural_comparisons/whole-protein-model-designs-20260928-v2');r=json.loads((root/'receipt.json').read_text());rh=sha(root/'receipt.json');assert r['status']=='complete_whole_protein_model_design_assessment_pending_independent_readback';bindings=dict(r['source_hashes'])
 for name,h in r['artifacts'].items():bindings[str(root/name)]=h
 for path,h in bindings.items():assert sha(path)==h,path
 source=Path('results/structural_comparisons/matched-whole-protein-records-20260928-v1');parts=json.loads((source/'partition_manifest.json').read_text());summary=pd.read_csv(source/'record_summary.tsv',sep='\t',dtype={'target_order':str,'background_order':str});selected=pd.read_csv('results/orthology/background-control-selection-20260927-v1/selections.tsv.gz',sep='\t',usecols=['target_id','background_id','policy','scenario_id']);assert len(selected)==2786912
 total=0;counts=Counter()
 with (root/'designs.jsonl').open() as f:
  for part in parts:
   values=pd.read_parquet(source/part['path']);records=selected.merge(values,on=['target_id','background_id'],validate='many_to_one');groups=records.groupby(KEYS).indices;grid=summary
   for k in SETTING:grid=grid.loc[grid[k].eq(part[k])]
   assert len(grid)==432;batch={}
   for _ in range(432*5):
    row=json.loads(next(f));key=tuple(row[k] for k in KEYS)+ (row['variant'],);assert key not in batch and all(row[k]==part[k] for k in SETTING);batch[key]=row
   assert set(batch)=={(*tuple(row[k] for k in KEYS),v) for row in grid.to_dict('records') for v in VARIANTS}
   for s in grid.to_dict('records'):
    group=tuple(s[k] for k in KEYS);frame=records.iloc[groups.get(group,[])];assert len(frame)==s['retained_records']
    for variant in VARIANTS:
     expected=reconstruct(frame,variant);expected.update(**{k:s[k] for k in KEYS+SETTING},variant=variant,screen_retained_records=len(frame),sequence_excluded_records=len(frame)-expected['records'])
     if variant.startswith('positive_log'):assert expected['records']==s['positive_sequence_pairs']
     compare(batch[(*group,variant)],expected);counts[expected['status']]+=1;total+=1
   print('Verified whole protein design setting',part['index']+1,'/96',flush=True)
  assert f.read()==''
 assert total==r['designs']==207360 and dict(counts)==r['counts']
 for path,h in bindings.items():assert sha(path)==h,path
 assert sha(root/'receipt.json')==rh and sha(lp)==lh
 result=dict(status='passed_full_whole_protein_model_design_readback',source_receipt_sha256=rh,checker_sha256=sha(__file__),producer_terminal_state=state,designs=total,counts=dict(counts),scope='All stratum memberships, five predictor designs, factored polynomial differences, positive-log exclusions, constant/rank/conditioning/marginal-overlap/reuse summaries independently reconstructed. Condition numbers use explicit numerical tolerance; no joint-support, fitted-effect or calibration claim.')
 with Path('metadata/whole_protein_model_designs_v2_completed_readback_20260928.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result),flush=True)
if __name__=='__main__':main()
