#!/usr/bin/env python3
"""Independently reconstruct all screened balance cells from audited fixed matches."""
import argparse,csv,json,subprocess,time
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from readback_background_control_balance import FEATURES,values,reconstruct,compare
from run_ortholog_pair_guide_comparison import sha

KEY=['guide','policy','scenario_id','mask','screen']
EXTRA=dict(metadata_matched_baseline_targets='baseline_targets',metadata_matched_target_mean='baseline_target_mean',retention_mean_shift='selection_mean_shift',retention_shift_in_metadata_matched_sd='selection_shift_in_baseline_sd',retention_shift_status='selection_shift_status')
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
 def verify():
  assert sha(args.plan)==ph
  for path,h in plan['pins'].items():assert sha(path)==h,path
 verify();dep=plan['producer']
 while psutil.pid_exists(dep['pid']):
  try:
   proc=psutil.Process(dep['pid'])
   if abs(proc.create_time()-dep['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
   assert proc.cmdline()==dep['cmdline']
  except psutil.NoSuchProcess:break
  time.sleep(30)
 state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 verify();root=Path(plan['source']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
 assert r['status']=='complete_screened_match_balance_pending_independent_readback'
 assert r['script_sha256']==plan['pins']['scripts/assess_screened_match_balance.py'] and r['balance_helper_sha256']==plan['pins']['scripts/background_control_balance.py']
 def sources():
  for path,h in r['source_hashes'].items():assert sha(path)==h,path
  for name,h in r['artifacts'].items():assert sha(root/name)==h,name
 sources();aroot=Path(plan['attrition']);ar=json.loads((aroot/'receipt.json').read_text());proof=json.loads(Path(plan['attrition_readback']).read_text())
 assert proof['status']=='passed_full_fixed_match_structural_attrition_readback' and proof['source_receipt_sha256']==sha(aroot/'receipt.json')
 for name,h in ar['artifacts'].items():assert sha(aroot/name)==h
 nodes={};vectors={}
 for kind in ['target','background']:
  path=Path(plan['graph'])/(kind+'_nodes.jsonl');assert sha(path)==ar['source_hashes'][str(path)]
  index={};vs={}
  for line in path.open():
   n=json.loads(line);k=n['node_id'];assert k not in index;vs[k]=values(n);index[k]={f:n[f] for f in ['guide','family','same_model']+(['taxon_id'] if kind=='target' else [])}
  nodes[kind]=index;vectors[kind]=pd.DataFrame.from_dict(vs,orient='index')
 baseline={g:vectors['target'].loc[[k for k,n in nodes['target'].items() if n['guide']==g]].to_numpy() for g in ['profile','mafft']}
 masks=['full','plddt70','both_masks'];screens=[s['id'] for s in ar['screens']]
 selected=pd.read_csv(aroot/'selection_eligibility.tsv.gz',sep='\t',usecols=['target_id','background_id',*KEY[:3]]+[side+'_'+mask+'_pass_bits' for side in ['target','background'] for mask in masks])
 assert len(selected)==r['selected_records']==ar['selected_records']
 assert not selected.duplicated(['target_id','policy','scenario_id']).any()
 for side,kind in [('target','target'),('background','background')]:
  assert selected[side+'_id'].map({k:n['guide'] for k,n in nodes[kind].items()}).equals(selected.guide)
 grouped=selected.groupby(KEY[:3]).indices
 def table(path,keys):
  result={}
  for row in csv.DictReader(path.open(),delimiter='\t'):
   k=tuple(row[f] for f in keys);assert k not in result;result[k]=row
  return result
 attr=table(aroot/'attrition_summary.tsv',KEY);cov=table(root/'coverage.tsv',KEY);bal=table(root/'balance.tsv',KEY+['feature'])
 groups={k[:3] for k in attr};expected={(*g,m,s) for g in groups for m in masks for s in screens}
 assert len(groups)==432 and set(grouped)<=groups and set(attr)==set(cov)==expected
 assert set(bal)=={(*k,f) for k in expected for f in FEATURES}
 for count,group in enumerate(sorted(groups),1):
  chunk=selected.iloc[grouped.get(group,[])];x=vectors['target'].loc[chunk.target_id].to_numpy();y=vectors['background'].loc[chunk.background_id].to_numpy();total=len(baseline[group[0]])
  for mask in masks:
   for bit,screen in enumerate(screens):
    # Decode each endpoint independently, then intersect boolean membership.
    keep=(chunk['target_'+mask+'_pass_bits'].to_numpy()//(2**bit)%2==1)&(chunk['background_'+mask+'_pass_bits'].to_numpy()//(2**bit)%2==1)
    retained=chunk.loc[keep];xx=x[keep];yy=y[keep];n=len(retained);key=(*group,mask,screen);e=attr[key]
    assert n==int(e['both_pass']) and total==int(e['original_targets']) and len(chunk)==int(e['metadata_matched']) and total-len(chunk)==int(e['metadata_unmatched'])
    tn=[nodes['target'][k] for k in retained.target_id];bn=[nodes['background'][k] for k in retained.background_id];assert not any(v['same_model'] for v in tn+bn)
    reuse=Counter(retained.background_id);expected_cov=dict(zip(KEY,key));expected_cov.update(original_targets=total,metadata_matched=len(chunk),metadata_unmatched=total-len(chunk),retained_matches=n,screen_excluded_matches=len(chunk)-n,retained_target_taxa=len({v['taxon_id'] for v in tn}),retained_target_families=len({v['family'] for v in tn}),distinct_controls=len(reuse),maximum_control_reuse=max(reuse.values(),default=0),top_five_control_fraction=sum(sorted(reuse.values(),reverse=True)[:5])/n if n else '',zero_sequence_distance_targets=int(sum(xx[:,0]==0)),zero_sequence_distance_controls=int(sum(yy[:,0]==0)))
    compare(cov[key],expected_cov)
    for j,f in enumerate(FEATURES):
     stats=reconstruct(xx[:,j],yy[:,j],baseline[group[0]][:,j]);retention=reconstruct(xx[:,j],yy[:,j],x[:,j]);expected_bal=dict(zip(KEY,key));expected_bal.update(feature=f,**stats,**{new:retention[old] for new,old in EXTRA.items()});compare(bal[(*key,f)],expected_bal)
  if count%12==0:print('Verified screened balance strata',count,'/432',flush=True)
 assert len(cov)==r['coverage_rows']==7776 and len(bal)==r['balance_rows']==62208 and r['features']==FEATURES
 verify();sources();assert sha(rp)==rh
 result=dict(status='passed_full_screened_match_balance_readback',plan_sha256=ph,source_receipt_sha256=rh,producer_terminal_state=state,selected_records_checked=len(selected),coverage_rows_checked=len(cov),balance_rows_checked=len(bal),scope='All fixed-match identities joined to source nodes; endpoint bits independently decoded; all feature moments, quantiles, nonestimable statuses, shifts against both baselines, taxon/family representation and reuse reconstructed without producer feature/statistics helpers. Descriptive balance only; no evolutionary effect or calibration claim.')
 with Path(plan['output']).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result),flush=True)
if __name__=='__main__':main()
