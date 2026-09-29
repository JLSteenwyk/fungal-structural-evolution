#!/usr/bin/env python3
"""Independently check every background coverage decision using integer ceiling cutoffs."""
import csv,json
from collections import Counter
from decimal import Decimal,ROUND_CEILING
from itertools import zip_longest
from pathlib import Path
from screen_duplication_alignment_reuse import sha

def main():
 root=Path('results/structural_comparisons/background-whole-protein-coverage-20260928-v1');rp=root/'receipt.json';r=json.loads(rp.read_text())
 for path,h in r['source_hashes'].items():assert sha(path)==h
 for name,h in r['artifacts'].items():assert sha(root/name)==h
 source=Path('results/structural_comparisons/background-candidate-measurements-20260928-v1/candidate_mask_measurements.tsv')
 seen=set();counts=Counter();exclusions=Counter();flags={};checked=0
 for original,row in zip_longest(csv.DictReader(source.open(),delimiter='\t'),csv.DictReader((root/'candidate_mask_coverage.tsv').open(),delimiter='\t')):
  assert original is not None and row is not None
  key=row['native_pair_key'],row['mask'];assert key not in seen;seen.add(key)
  for name in ['native_pair_key','gene_a','gene_b','taxon_a','taxon_b','pair_key','mask','measurement_disposition','candidate_and_native_ortholog_both','both_guides_and_parents_unreported']:assert row[name]==original[name]
  flags[key]={}
  for spec in r['screens']:
   why=[];name=spec['id'];disposition=original['measurement_disposition']
   if disposition not in ['new_model_pair','already_in_reference_queue']:why=[disposition]
   else:
    lengths=[]
    for order in [0,1]:
     if original[f'fit_order{order}_status']!='aligned':why.append(f'order{order}_not_numerically_usable')
     else:lengths.append(int(Decimal(original[f'fit_order{order}_aligned_length'])))
    if any(n<spec['minimum_aligned_residues'] for n in lengths):why.append('short_alignment')
    required=max(int((Decimal(original['length_'+side])*Decimal(str(spec['minimum_original_coverage']))).to_integral_value(rounding=ROUND_CEILING)) for side in ['a','b'])
    if any(n<required for n in lengths):why.append('low_original_protein_coverage')
   assert row[name+'_exclusions']==';'.join(why) and int(row[name+'_pass'])==int(not why)
   counts[key[1]+':'+name]+=not why
   for reason in why:exclusions[key[1]+':'+name+':'+reason]+=1
   flags[key][name]=not why;checked+=1
 keys={k[0] for k in seen};assert len(seen)==r['candidate_mask_rows']==156744 and seen=={(k,m) for k in keys for m in ['full','plddt70']}
 assert dict(counts)==r['pass_counts'] and dict(exclusions)==r['exclusion_counts']
 both={s['id']:sum(flags[k,'full'][s['id']] and flags[k,'plddt70'][s['id']] for k in keys) for s in r['screens']};assert both==r['both_masks_pass_counts']
 result=dict(status='passed_full_background_original_protein_coverage_readback',receipt_sha256=sha(rp),checker_sha256=sha(__file__),candidate_mask_rows=len(seen),screen_decisions_checked=checked,pass_counts=dict(counts),both_masks_pass_counts=both,scope='Every row, original identity/flags, six screen decisions and overlapping exclusion lists independently replayed with decimal ceiling coverage cutoffs. Exact both-mask intersections checked; no scientific acceptance.')
 with Path('metadata/background_whole_protein_coverage_completed_20260928.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
