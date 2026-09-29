#!/usr/bin/env python3
"""Check every target coverage row independently against original nodes and summaries."""
import csv,json
from collections import Counter
from decimal import Decimal,ROUND_CEILING
from pathlib import Path
from screen_duplication_alignment_reuse import sha

def main():
 root=Path('results/structural_comparisons/matched-target-coverage-20260928-v1');rp=root/'receipt.json';r=json.loads(rp.read_text())
 for path,h in r['source_hashes'].items():assert sha(path)==h
 for name,h in r['artifacts'].items():assert sha(root/name)==h
 nodes={}
 for line in Path('results/orthology/background-match-graph-20260927-v1/target_nodes.jsonl').open():
  row=json.loads(line);assert row['node_id'] not in nodes;nodes[row['node_id']]=row
 summary={(x['pair_key'],x['mask']):x for x in csv.DictReader(Path('results/structural_comparisons/primary-usable-orders-20260927-v1/pair_mask_order_summary.tsv').open(),delimiter='\t')}
 seen=set();counts=Counter();both=Counter();same=Counter();flags={};decisions=0
 for row in csv.DictReader((root/'target_mask_coverage.tsv').open(),delimiter='\t'):
  key=row['target_id'],row['mask'];assert key not in seen;seen.add(key);node=nodes[key[0]]
  for field in ['guide','family','gene_a','gene_b','taxon_id','pair_key','same_model']:assert row[field]==str(node[field])
  flags[key]={}
  if node['same_model']:same[row['guide']+':'+row['mask']]+=1
  for spec in r['screens']:
   why=[]
   if node['same_model']:why=['identical_model_no_alignment']
   else:
    n=summary[node['pair_key'],row['mask']];lengths=[]
    for order in [0,1]:
     if n[f'order{order}_status']!='aligned':why.append(f'order{order}_not_numerically_usable')
     else:lengths.append(int(Decimal(n[f'order{order}_aligned_length'])))
    if any(v<spec['minimum_aligned_residues'] for v in lengths):why.append('short_alignment')
    threshold=max(int((Decimal(node['length_'+s])*Decimal(str(spec['minimum_original_coverage']))).to_integral_value(rounding=ROUND_CEILING)) for s in ['a','b'])
    if any(v<threshold for v in lengths):why.append('low_original_protein_coverage')
   name=spec['id'];assert row[name+'_exclusions']==';'.join(why) and int(row[name+'_pass'])==int(not why)
   flags[key][name]=not why;counts[row['guide']+':'+row['mask']+':'+name]+=not why;decisions+=1
 assert seen=={(key,mask) for key in nodes for mask in ['full','plddt70']} and len(seen)==r['target_mask_rows']
 for key,node in nodes.items():
  for spec in r['screens']:
   name=spec['id'];both[node['guide']+':'+name]+=flags[key,'full'][name] and flags[key,'plddt70'][name]
 assert dict(counts)==r['pass_counts'] and dict(both)==r['both_masks_pass_counts'] and dict(same)==r['same_model_counts']
 result=dict(status='passed_full_matching_target_coverage_readback',receipt_sha256=sha(rp),checker_sha256=sha(__file__),targets=len(nodes),target_mask_rows=len(seen),screen_decisions_checked=decisions,pass_counts=dict(counts),both_masks_pass_counts=dict(both),scope='Every target identity and six coverage decisions independently replayed with decimal ceiling cutoffs against original nodes and audited summaries. Same-model exclusions and exact both-mask intersections checked; selected-match attrition and evolutionary inference remain pending.')
 with Path('metadata/matched_target_coverage_completed_20260928.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
