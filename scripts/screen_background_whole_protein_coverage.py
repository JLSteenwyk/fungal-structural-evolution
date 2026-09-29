#!/usr/bin/env python3
"""Apply existing six original-protein coverage screens to every background candidate."""
import csv,json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from screen_duplication_alignment_reuse import sha

def reasons(row,spec):
 if row['measurement_disposition'] not in ['new_model_pair','already_in_reference_queue']:return [row['measurement_disposition']]
 result=[];lengths=[]
 for order in [0,1]:
  if row[f'fit_order{order}_status']!='aligned':result.append(f'order{order}_not_numerically_usable')
  else:lengths.append(int(float(row[f'fit_order{order}_aligned_length'])))
 cutoff=Fraction(str(spec['minimum_original_coverage']));original=[int(row['length_a']),int(row['length_b'])]
 if any(x<=0 for x in original):raise ValueError('Invalid original protein length')
 if any(n>min(original) for n in lengths):raise ValueError('Aligned length exceeds protein')
 if lengths and min(lengths)<spec['minimum_aligned_residues']:result.append('short_alignment')
 if any(n*cutoff.denominator<d*cutoff.numerator for n in lengths for d in original):result.append('low_original_protein_coverage')
 return result

def main():
 source=Path('results/structural_comparisons/background-candidate-measurements-20260928-v1');rp=source/'receipt.json';r=json.loads(rp.read_text());cp=Path('metadata/background_candidate_measurements_completed_20260928.json');c=json.loads(cp.read_text());assert c['receipt_sha256']==sha(rp)
 table=source/'candidate_mask_measurements.tsv';assert sha(table)==r['artifacts'][table.name]
 pp=Path('metadata/whole_protein_common_fits_plan_20260927.json');screens=json.loads(pp.read_text())['screens']
 bindings={str(p):sha(p) for p in [rp,cp,table,pp]}
 out=Path('results/structural_comparisons/background-whole-protein-coverage-20260928-v1');out.mkdir(exist_ok=False)
 fields=['native_pair_key','gene_a','gene_b','taxon_a','taxon_b','pair_key','mask','measurement_disposition','candidate_and_native_ortholog_both','both_guides_and_parents_unreported']
 fields+=[s['id']+suffix for s in screens for suffix in ['_pass','_exclusions']]
 counts=Counter();exclusions=Counter();flags={};seen=set()
 with (out/'candidate_mask_coverage.tsv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
  for row in csv.DictReader(table.open(),delimiter='\t'):
   key=row['native_pair_key'],row['mask'];assert key not in seen;seen.add(key)
   result={k:row[k] for k in fields if k in row};flags[key]={}
   for spec in screens:
    name=spec['id'];why=reasons(row,spec);passed=not why;result[name+'_pass']=int(passed);result[name+'_exclusions']=';'.join(why);flags[key][name]=passed
    counts[row['mask']+':'+name]+=passed
    for reason in why:exclusions[row['mask']+':'+name+':'+reason]+=1
   writer.writerow(result)
 assert len(seen)==r['candidate_mask_rows']==156744
 candidates={k[0] for k in seen};assert len(candidates)==78372 and seen=={(key,mask) for key in candidates for mask in ['full','plddt70']}
 joint={s['id']:sum(flags[key,'full'][s['id']] and flags[key,'plddt70'][s['id']] for key in candidates) for s in screens}
 for p,h in bindings.items():assert sha(p)==h
 result=dict(status='complete_background_original_protein_coverage_screen_pending_readback',source_hashes=bindings,script_sha256=sha(__file__),screens=screens,candidates=len(candidates),candidate_mask_rows=len(seen),pass_counts=dict(counts),both_masks_pass_counts=joint,exclusion_counts=dict(exclusions),artifacts={'candidate_mask_coverage.tsv':sha(out/'candidate_mask_coverage.tsv')},scope='Existing30/50-residue and50/70/90-percent original-protein coverage screens; both native orders must be numerically usable and pass each threshold. Exact integer fraction comparisons; confidence-masked coverage uses original full protein lengths, never retained lengths. All unmeasured candidates and overlapping exclusions retained. Both-mask summaries are sensitivity counts, not independent replicates, domain-orientation validation or evolutionary significance. Confidence calibration and phylogenetic matching remain downstream.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['exclusion_counts','source_hashes']},indent=2))
if __name__=='__main__':main()
