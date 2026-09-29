#!/usr/bin/env python3
"""Independently verify matched outcome identity, source fields and species contrasts."""
import csv,gzip,hashlib,json,subprocess
from collections import Counter
from pathlib import Path
from screen_duplication_alignment_reuse import sha

def main():
 unit='fungal-matched-structural-measurements-20260928.service';state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 root=Path('results/structural_comparisons/matched-structural-measurements-20260928-v1');r=json.loads((root/'receipt.json').read_text());rh=sha(root/'receipt.json');assert r['status']=='complete_matched_structural_measurement_assembly_pending_independent_readback';bindings=dict(r['source_hashes'])
 for n,h in r['artifacts'].items():bindings[str(root/n)]=h
 for p,h in bindings.items():assert sha(p)==h,p
 def records(path):
  with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:yield from csv.DictReader(f,delimiter='\t')
 species=Path('results/phylogeny/matched-species-contrasts-20260927-v1');taxa=json.loads((species/'taxa.json').read_text());index={t:i for i,t in enumerate(taxa)};pairs={}
 for row in records(species/'selected_pair_patterns.tsv.gz'):
  key=row['target_id'],row['background_id'];assert key not in pairs;pairs[key]=row
 nodes={}
 for kind,pos in [('target',0),('background',1)]:
  wanted={k[pos] for k in pairs};nodes[kind]={}
  for line in Path('results/orthology/background-match-graph-20260927-v1',kind+'_nodes.jsonl').open():
   n=json.loads(line)
   if n['node_id'] in wanted:assert n['node_id'] not in nodes[kind];nodes[kind][n['node_id']]=n
  assert set(nodes[kind])==wanted
 target={}
 for row in records('results/structural_comparisons/matching-target-measurements-20260928-v1/target_mask_measurements.tsv.gz'):
  if row['node_id'] in nodes['target']:key=row['node_id'],row['mask'];assert key not in target;target[key]=row
 background={};genes={(n['gene_a'],n['gene_b']) for n in nodes['background'].values()}
 for row in records('results/structural_comparisons/background-candidate-measurements-20260928-v1/candidate_mask_measurements.tsv'):
  if (row['gene_a'],row['gene_b']) in genes:key=row['gene_a'],row['gene_b'],row['mask'];assert key not in background;background[key]=row
 selected=Counter();flags={};flagfields=[side+'_'+mask+'_pass_bits' for side in ['target','background'] for mask in ['full','plddt70','both_masks']]
 for row in records('results/structural_comparisons/fixed-match-structural-attrition-20260928-v1/selection_eligibility.tsv.gz'):
  key=row['target_id'],row['background_id'];v=tuple(row[k] for k in flagfields);assert key not in flags or flags[key]==v;flags[key]=v;selected[key]+=1
 assert set(selected)==set(pairs) and sum(selected.values())==r['original_selection_uses']==2786912
 seen=set();counts=Counter();checked=0
 for row in records(root/'unique_pair_mask_measurements.tsv.gz'):
  key=row['target_id'],row['background_id'];mask=row['mask'];assert (*key,mask) not in seen and mask in ['full','plddt70'];seen.add((*key,mask));original=pairs[key]
  tn=nodes['target'][key[0]];bn=nodes['background'][key[1]];assert tn['guide']==bn['guide'] and tn['family']==bn['family']
  assert all(row[k]==v for k,v in original.items()) and int(row['selected_record_uses'])==selected[key]
  ti=index[tn['taxon_id']];ai=index[bn['taxon_a']];bi=index[bn['taxon_b']]
  assert [int(row[k]) for k in ['focal_index','background_index_a','background_index_b']]==[ti,ai,bi]
  weights=Counter();weights[ti]+=2;weights[ai]-=1;weights[bi]-=1;pattern=sorted((i,w) for i,w in weights.items() if w);assert sum(w for _,w in pattern)==0
  assert row['species_pattern_id']==hashlib.sha256(json.dumps(pattern,separators=(',',':')).encode()).hexdigest()
  expected=dict(guide=tn['guide'],family=tn['family'],focal_taxon=tn['taxon_id'],background_taxon_a=bn['taxon_a'],background_taxon_b=bn['taxon_b'],target_gene_node=tn['gene_node'],mask=mask,target_pair_key=tn['pair_key'],background_pair_key=bn['pair_key'],target_sequence_distance=str(tn['sequence_distance']),background_sequence_distance=str(bn['sequence_distance']),**dict(zip(flagfields,flags[key])))
  t=target[key[0],mask];b=background[bn['gene_a'],bn['gene_b'],mask]
  expected.update(target_endpoint_order=t['target_endpoint_order'],background_endpoint_order=b['candidate_endpoint_order'])
  fields=['measurement_disposition','canonical_model_a','canonical_version_a','canonical_model_b','canonical_version_b']+[k for k in t if k.startswith('fit_')]
  for role,source in [('target',t),('background',b)]:
   for field in fields:expected[role+'_'+field]=source[field];checked+=1
  assert set(row)==set(original)|set(expected)
  for k,v in expected.items():assert row[k]==v,(key,mask,k)
  counts[mask+':'+t['measurement_disposition']+':'+b['measurement_disposition']]+=1
  if len(seen)%25000==0:print('Verified matched measurement rows',len(seen),'/105350',flush=True)
 assert seen=={(*k,m) for k in pairs for m in ['full','plddt70']} and len(seen)==r['pair_mask_rows']==105350 and len(pairs)==r['unique_pairs']==52675 and dict(counts)==r['counts']
 for p,h in bindings.items():assert sha(p)==h,p
 assert sha(root/'receipt.json')==rh
 result=dict(status='passed_full_matched_structural_measurement_readback',source_receipt_sha256=rh,checker_sha256=sha(__file__),producer_terminal_state=state,unique_pairs=len(pairs),pair_mask_rows=len(seen),source_measurement_fields_checked=checked,original_selection_uses=sum(selected.values()),scope='All source structural fields, node identities, eligibility flags and original selection multiplicities checked. Every species contrast index/hash independently reconstructed from focal/background taxa. Exact full pair/mask membership; no independent-replicate, weighting, effect, calibrated covariance or biological inference claim.')
 with Path('metadata/matched_structural_measurements_completed_readback_20260928.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result),flush=True)
if __name__=='__main__':main()
