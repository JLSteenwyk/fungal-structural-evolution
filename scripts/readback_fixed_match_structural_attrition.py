#!/usr/bin/env python3
"""Replay all fixed-match flags and attrition cells from the original coverage sources."""
import csv,gzip,json,subprocess
from collections import defaultdict,Counter
from itertools import zip_longest
from pathlib import Path
import numpy as np
from screen_duplication_alignment_reuse import sha

def main():
 unit='fungal-fixed-match-structural-attrition-20260928.service';state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 root=Path('results/structural_comparisons/fixed-match-structural-attrition-20260928-v1');rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
 for path,h in r['source_hashes'].items():assert sha(path)==h
 for name,h in r['artifacts'].items():assert sha(root/name)==h
 screens=[x['id'] for x in r['screens']];masks=['full','plddt70','both_masks'];tables={}
 for kind,path in [('target','results/structural_comparisons/matched-target-coverage-20260928-v1/target_mask_coverage.tsv'),('background','results/structural_comparisons/background-whole-protein-coverage-20260928-v1/candidate_mask_coverage.tsv')]:
  index={}
  for row in csv.DictReader(open(path),delimiter='\t'):
   key=(row['target_id'],row['mask']) if kind=='target' else (row['gene_a'],row['gene_b'],row['mask'])
   assert key not in index;values=tuple(int(row[s+'_pass']) for s in screens);assert set(values)<={0,1};index[key]=(values,row['pair_key'])
  tables[kind]=index
 nodes={}
 for kind in ['target','background']:
  index={}
  for line in open('results/orthology/background-match-graph-20260927-v1/'+kind+'_nodes.jsonl'):
   row=json.loads(line);values=[]
   for mask in masks[:2]:
    key=(row['node_id'],mask) if kind=='target' else (row['gene_a'],row['gene_b'],mask)
    flags,pair=tables[kind][key];assert pair==row['pair_key'];values.append(flags)
   values.append(tuple(int(a and b) for a,b in zip(*values)))
   bits=tuple(sum(v*2**i for i,v in enumerate(flags)) for flags in values)
   assert row['node_id'] not in index;index[row['node_id']]=(row['guide'],tuple(v for flags in values for v in flags),bits)
  nodes[kind]=index
 totals=defaultdict(lambda:np.zeros((18,4),dtype=np.int64));selected=Counter();n=0;positions=np.arange(18);class_cache={}
 selection=Path('results/orthology/background-control-selection-20260927-v1')
 with gzip.open(selection/'selections.tsv.gz','rt') as f,gzip.open(root/'selection_eligibility.tsv.gz','rt') as g:
  for original,row in zip_longest(csv.DictReader(f,delimiter='\t'),csv.DictReader(g,delimiter='\t')):
   assert original is not None and row is not None and all(row[k]==v for k,v in original.items())
   t=nodes['target'][original['target_id']];b=nodes['background'][original['background_id']];assert row['guide']==t[0]==b[0]
   for side,node in [('target',t),('background',b)]:
    for mask,value in zip(masks,node[2]):assert int(row[side+'_'+mask+'_pass_bits'])==value
   signature=t[1],b[1]
   if signature not in class_cache:class_cache[signature]=np.array([a+2*c for a,c in zip(*signature)],dtype=np.int64)
   group=row['guide'],row['policy'],row['scenario_id'];totals[group][positions,class_cache[signature]]+=1;selected[group]+=1;n+=1
   if n%250000==0:print('Verified fixed selections',n,'/',r['selected_records'],flush=True)
 assert n==r['selected_records']==2786912
 sr=json.load(open(selection/'receipt.json'));groups={tuple(k.split('|')[:3]) for k in sr['counts'] if k.endswith('|targets')};assert len(groups)==432
 seen=set()
 for row in csv.DictReader((root/'attrition_summary.tsv').open(),delimiter='\t'):
  group=row['guide'],row['policy'],row['scenario_id'];assert group in groups;key=(*group,row['mask'],row['screen']);assert key not in seen;seen.add(key)
  i=masks.index(row['mask'])*6+screens.index(row['screen']);v=totals[group][i];prefix='|'.join(group)
  assert int(row['original_targets'])==sr['counts'][prefix+'|targets']
  assert int(row['metadata_matched'])==selected[group]==sr['counts'].get(prefix+'|matched',0)
  assert int(row['metadata_unmatched'])==sr['counts'].get(prefix+'|unmatched',0)
  for field,value in zip(['neither_pass','target_only_pass','background_only_pass','both_pass'],v):assert int(row[field])==value
  assert sum(v)==selected[group] and int(row['metadata_matched'])+int(row['metadata_unmatched'])==int(row['original_targets'])
 assert seen=={(*group,mask,screen) for group in groups for mask in masks for screen in screens} and len(seen)==r['summary_rows']==7776
 assert sha(rp)==rh
 for name,h in r['artifacts'].items():assert sha(root/name)==h
 result=dict(status='passed_full_fixed_match_structural_attrition_readback',source_receipt_sha256=rh,checker_sha256=sha(__file__),producer_terminal_state=state,selected_records_checked=n,eligibility_bitsets_checked=n*6,summary_cells_checked=len(seen),scope='All original selection fields preserved, every target/control flag independently reconstructed from original screen tables and graph identities, all7776 cells including empty strata and unmatched denominators independently aggregated. No rematching, independence, calibration or evolutionary effect claim.')
 with Path('metadata/fixed_match_structural_attrition_completed_readback_20260928.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
