#!/usr/bin/env python3
"""Measure structural eligibility loss in every frozen metadata-matching scenario."""
import csv,gzip,json
from collections import Counter,defaultdict
from pathlib import Path
from screen_duplication_alignment_reuse import sha

def main():
 bindings={}
 def bind(p):p=Path(p);bindings[str(p)]=sha(p);return p
 def js(p):return json.loads(bind(p).read_text())
 selection=Path('results/orthology/background-control-selection-20260927-v1');sc=js('metadata/background_control_selection_completed_20260927.json');sr=js(selection/'receipt.json');assert sc['source_receipt_sha256']==sha(selection/'receipt.json')
 assert sha(bind(sc['readback']))==sc['readback_sha256']
 for name,h in sr['artifacts'].items():assert sha(bind(selection/name))==h
 roots={'target':Path('results/structural_comparisons/matched-target-coverage-20260928-v1'),'background':Path('results/structural_comparisons/background-whole-protein-coverage-20260928-v1')}
 evidence={'target':'metadata/matched_target_coverage_completed_20260928.json','background':'metadata/background_whole_protein_coverage_completed_20260928.json'}
 indexes={};screens=None
 for kind,root in roots.items():
  r=js(root/'receipt.json');proof=js(evidence[kind]);assert proof['receipt_sha256']==sha(root/'receipt.json')
  if screens is None:screens=r['screens']
  else:assert screens==r['screens']
  file='target_mask_coverage.tsv' if kind=='target' else 'candidate_mask_coverage.tsv';path=bind(root/file);assert sha(path)==r['artifacts'][file];index={}
  for row in csv.DictReader(path.open(),delimiter='\t'):
   key=(row['target_id'],row['mask']) if kind=='target' else (row['gene_a'],row['gene_b'],row['mask'])
   assert key not in index;bits=0
   for i,s in enumerate(screens):
    value=int(row[s['id']+'_pass']);assert value in [0,1];bits|=value<<i
   index[key]=(bits,row['pair_key'])
  indexes[kind]=index
 spath=Path('metadata/background_control_selection_plan_20260927.json');sp=js(spath);assert sr['plan_sha256']==sha(spath)
 graph=Path(sp['graph']);gr=js(graph/'receipt.json');assert sha(graph/'receipt.json')==sp['pins'][str(graph/'receipt.json')];nodes={}
 for kind in ['target','background']:
  path=bind(graph/(kind+'_nodes.jsonl'));assert sha(path)==gr['artifacts'][path.name]==sp['pins'][str(path)]
  index={}
  for line in path.open():
   row=json.loads(line);key=row['node_id'];assert key not in index
   masks=[]
   for mask in ['full','plddt70']:
    k=(key,mask) if kind=='target' else (row['gene_a'],row['gene_b'],mask)
    value,pair=indexes[kind][k];assert pair==row['pair_key'];masks.append(value)
   index[key]=(row['guide'],masks[0],masks[1],masks[0]&masks[1])
  nodes[kind]=index
 out=Path('results/structural_comparisons/fixed-match-structural-attrition-20260928-v1');out.mkdir(exist_ok=False)
 hist=defaultdict(Counter);selected=Counter();n=0
 extras=['guide']+[kind+'_'+mask+'_pass_bits' for kind in ['target','background'] for mask in ['full','plddt70','both_masks']]
 with gzip.open(selection/'selections.tsv.gz','rt') as f,gzip.open(out/'selection_eligibility.tsv.gz','wt',compresslevel=1) as dest:
  reader=csv.DictReader(f,delimiter='\t');writer=csv.DictWriter(dest,fieldnames=reader.fieldnames+extras,delimiter='\t',lineterminator='\n');writer.writeheader()
  for row in reader:
   t=nodes['target'][row['target_id']];b=nodes['background'][row['background_id']];assert t[0]==b[0];group=t[0],row['policy'],row['scenario_id'];selected[group]+=1
   result=dict(row,guide=t[0])
   for i,mask in enumerate(['full','plddt70','both_masks'],1):
    result['target_'+mask+'_pass_bits']=t[i];result['background_'+mask+'_pass_bits']=b[i];hist[(*group,mask)][(t[i],b[i])]+=1
   writer.writerow(result);n+=1
   if n%250000==0:print('Linked frozen selections',n,'/',sr['selected_records'],flush=True)
 assert n==sr['selected_records']==2786912
 groups=[tuple(k.split('|')[:3]) for k in sr['counts'] if k.endswith('|targets')];assert len(groups)==432 and len(set(groups))==432
 summary=[]
 for group in sorted(groups):
  prefix='|'.join(group);total=sr['counts'][prefix+'|targets'];matched=sr['counts'].get(prefix+'|matched',0);unmatched=sr['counts'].get(prefix+'|unmatched',0);assert total==matched+unmatched and selected[group]==matched
  for mask in ['full','plddt70','both_masks']:
   h=hist[(*group,mask)];assert sum(h.values())==matched
   for i,s in enumerate(screens):
    cells=Counter()
    for (t,b),count in h.items():cells[bool(t&(1<<i)),bool(b&(1<<i))]+=count
    summary.append(dict(guide=group[0],policy=group[1],scenario_id=group[2],mask=mask,screen=s['id'],original_targets=total,metadata_matched=matched,metadata_unmatched=unmatched,both_pass=cells[True,True],target_only_pass=cells[True,False],background_only_pass=cells[False,True],neither_pass=cells[False,False]))
 with (out/'attrition_summary.tsv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(summary)
 for path,h in bindings.items():assert sha(path)==h
 result=dict(status='complete_fixed_match_structural_attrition_pending_independent_readback',script_sha256=sha(__file__),source_hashes=bindings,selected_records=n,matching_strata=len(groups),summary_rows=len(summary),screens=screens,bit_encoding='Bit i corresponds to screens[i]; both_masks bitsets are intersections of full and plddt70.',artifacts={name:sha(out/name) for name in ['selection_eligibility.tsv.gz','attrition_summary.tsv']},scope='All original metadata-selected matches retained unchanged, with target/control eligibility under every screen and mask; all432 strata including empty cells retain original target and metadata-unmatched denominators. No rematching based on structural outcomes. Attrition is descriptive and can change the estimand; retained matches are not independent or automatically biologically qualified.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print('Completed',n,'fixed selections;',len(summary),'summary cells',flush=True)
if __name__=='__main__':main()
