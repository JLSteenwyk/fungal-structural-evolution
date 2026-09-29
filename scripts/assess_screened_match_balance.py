#!/usr/bin/env python3
"""Assess covariate balance and representation after fixed-match structural screening."""
import csv,gzip,json,subprocess
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from background_control_balance import FEATURES,feature_vector,balance
from screen_duplication_alignment_reuse import sha

def main():
 bindings={}
 def js(path):path=Path(path);bindings[str(path)]=sha(path);return json.loads(path.read_text())
 cp=Path('metadata/fixed_match_structural_attrition_completed_20260928.json');c=js(cp)
 root=Path('results/structural_comparisons/fixed-match-structural-attrition-20260928-v1');r=js(root/'receipt.json');assert c['receipt_sha256']==sha(root/'receipt.json')
 proof=Path('metadata/fixed_match_structural_attrition_completed_readback_20260928.json');assert sha(proof)==c['readback_sha256'];a=js(proof);assert a['source_receipt_sha256']==sha(root/'receipt.json')
 for name,h in r['artifacts'].items():assert sha(root/name)==h;bindings[str(root/name)]=h
 graph=Path('results/orthology/background-match-graph-20260927-v1');nodes={};features={}
 for role in ['target','background']:
  path=graph/(role+'_nodes.jsonl');assert sha(path)==r['source_hashes'][str(path)];bindings[str(path)]=sha(path)
  index={};vectors={}
  for line in path.open():
   n=json.loads(line);assert n['node_id'] not in index;index[n['node_id']]={k:n[k] for k in ['guide','family','same_model']+(['taxon_id'] if role=='target' else ['taxon_a','taxon_b'])};vectors[n['node_id']]=feature_vector(n)
  nodes[role]=index;features[role]=vectors
 baseline={g:np.array([features['target'][k] for k,v in nodes['target'].items() if v['guide']==g]) for g in ['profile','mafft']}
 groups=defaultdict(list);n=0
 with gzip.open(root/'selection_eligibility.tsv.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   key=row['guide'],row['policy'],row['scenario_id'];tid,bid=row['target_id'],row['background_id'];assert nodes['target'][tid]['guide']==nodes['background'][bid]['guide']==key[0]
   bits=tuple(int(row['target_'+mask+'_pass_bits'])&int(row['background_'+mask+'_pass_bits']) for mask in ['full','plddt70','both_masks'])
   groups[key].append((tid,bid,*bits));n+=1
 assert n==r['selected_records']
 expected={}
 for row in csv.DictReader((root/'attrition_summary.tsv').open(),delimiter='\t'):
  key=tuple(row[k] for k in ['guide','policy','scenario_id','mask','screen']);assert key not in expected;expected[key]=row
 strata=sorted({k[:3] for k in expected});assert len(strata)==432
 out=Path('results/structural_comparisons/screened-match-balance-20260928-v1');out.mkdir(exist_ok=False);coverage_rows=0;balance_rows=0
 with (out/'coverage.tsv').open('w') as cf,(out/'balance.tsv').open('w') as bf:
  cw=bw=None
  for group in strata:
   records=groups.get(group,[]);x=np.array([features['target'][v[0]] for v in records]).reshape(-1,len(FEATURES));y=np.array([features['background'][v[1]] for v in records]).reshape(-1,len(FEATURES));flags=np.array([v[2:] for v in records],dtype=np.int64).reshape(-1,3)
   for j,mask in enumerate(['full','plddt70','both_masks']):
    for i,spec in enumerate(r['screens']):
     key=(*group,mask,spec['id']);e=expected[key];idx=np.flatnonzero(flags[:,j]&(1<<i));xx=x[idx];yy=y[idx];assert len(idx)==int(e['both_pass'])
     tids=[records[k][0] for k in idx];bids=[records[k][1] for k in idx];reuse=Counter(bids)
     coverage=dict(guide=group[0],policy=group[1],scenario_id=group[2],mask=mask,screen=spec['id'],original_targets=int(e['original_targets']),metadata_matched=int(e['metadata_matched']),metadata_unmatched=int(e['metadata_unmatched']),retained_matches=len(idx),screen_excluded_matches=len(records)-len(idx),retained_target_taxa=len({nodes['target'][k]['taxon_id'] for k in tids}),retained_target_families=len({nodes['target'][k]['family'] for k in tids}),distinct_controls=len(reuse),maximum_control_reuse=max(reuse.values(),default=0),top_five_control_fraction=sum(sorted(reuse.values(),reverse=True)[:5])/len(idx) if len(idx) else '',zero_sequence_distance_targets=int(sum(xx[:,0]==0)),zero_sequence_distance_controls=int(sum(yy[:,0]==0)))
     assert not any(nodes['target'][k]['same_model'] for k in tids) and not any(nodes['background'][k]['same_model'] for k in bids)
     if cw is None:cw=csv.DictWriter(cf,fieldnames=list(coverage),delimiter='\t',lineterminator='\n');cw.writeheader()
     cw.writerow(coverage);coverage_rows+=1
     for col,feature in enumerate(FEATURES):
      stats=balance(xx[:,col],yy[:,col],baseline[group[0]][:,col]);retention=balance(xx[:,col],yy[:,col],x[:,col])
      extra={new:retention[old] for new,old in [('metadata_matched_baseline_targets','baseline_targets'),('metadata_matched_target_mean','baseline_target_mean'),('retention_mean_shift','selection_mean_shift'),('retention_shift_in_metadata_matched_sd','selection_shift_in_baseline_sd'),('retention_shift_status','selection_shift_status')]}
      row=dict(guide=group[0],policy=group[1],scenario_id=group[2],mask=mask,screen=spec['id'],feature=feature,**stats,**extra)
      if bw is None:bw=csv.DictWriter(bf,fieldnames=list(row),delimiter='\t',lineterminator='\n');bw.writeheader()
      bw.writerow(row);balance_rows+=1
   print('Screened balance strata',coverage_rows//18,'/432',flush=True)
 assert coverage_rows==7776 and balance_rows==62208
 for path,h in bindings.items():assert sha(path)==h
 result=dict(status='complete_screened_match_balance_pending_independent_readback',source_hashes=bindings,script_sha256=sha(__file__),balance_helper_sha256=sha('scripts/background_control_balance.py'),selected_records=n,coverage_rows=coverage_rows,balance_rows=balance_rows,features=FEATURES,artifacts={name:sha(out/name) for name in ['coverage.tsv','balance.tsv']},scope='All432 fixed matching strata x3 masks x6 screens. Eight covariate summaries, shifts relative to original modeled targets and metadata-matched targets, retained taxa/families and control reuse. Zero-distance log exclusions and nonestimable variance explicit; no automatic balance threshold, independence, phylogenetic correction or effect inference.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print('Complete screened balance; independent readback required',flush=True)
if __name__=='__main__':main()
