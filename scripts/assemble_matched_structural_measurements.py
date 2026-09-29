#!/usr/bin/env python3
"""Join fixed matched identities, structural measurements and species contrast IDs."""
import csv,gzip,json,subprocess
from collections import Counter
from pathlib import Path
from screen_duplication_alignment_reuse import sha

def main():
 bindings={}
 def load(path):path=Path(path);bindings[str(path)]=sha(path);return json.loads(path.read_text())
 def bindroot(root):
  r=load(root/'receipt.json')
  for n,h in r['artifacts'].items():assert sha(root/n)==h;bindings[str(root/n)]=h
  return r
 unit='fungal-matching-target-measurement-readback-20260928.service';state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 tr=Path('results/structural_comparisons/matching-target-measurements-20260928-v1');tproof=load('metadata/matching_target_measurements_completed_readback_20260928.json');bindroot(tr);assert tproof['status']=='passed_full_matching_target_measurement_readback' and tproof['source_receipt_sha256']==sha(tr/'receipt.json')
 br=Path('results/structural_comparisons/background-candidate-measurements-20260928-v1');bproof=load('metadata/background_candidate_measurements_completed_20260928.json');bindroot(br);assert bproof['receipt_sha256']==sha(br/'receipt.json')
 species=Path('results/phylogeny/matched-species-contrasts-20260927-v1');sr=bindroot(species);assert sr['status']=='complete_matched_species_contrast_design_with_full_pair_readback'
 attr=Path('results/structural_comparisons/fixed-match-structural-attrition-20260928-v1');ar=bindroot(attr);ap=load('metadata/fixed_match_structural_attrition_completed_readback_20260928.json');assert ap['source_receipt_sha256']==sha(attr/'receipt.json')
 pairs={}
 with gzip.open(species/'selected_pair_patterns.tsv.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   key=row['target_id'],row['background_id'];assert key not in pairs;pairs[key]=row
 assert len(pairs)==sr['unique_selected_pairs']==52675
 selected=Counter();flags={};flagfields=[side+'_'+mask+'_pass_bits' for side in ['target','background'] for mask in ['full','plddt70','both_masks']]
 with gzip.open(attr/'selection_eligibility.tsv.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   key=row['target_id'],row['background_id'];v={k:row[k] for k in flagfields};assert key not in flags or flags[key]==v;flags[key]=v;selected[key]+=1
 assert set(selected)==set(pairs) and sum(selected.values())==sr['selected_records']==2786912
 for key,row in pairs.items():assert selected[key]==int(row['selected_record_uses'])
 graph=Path('results/orthology/background-match-graph-20260927-v1');gr=bindroot(graph);nodes={}
 for kind,index in [('target',0),('background',1)]:
  wanted={key[index] for key in pairs};found={}
  for line in (graph/(kind+'_nodes.jsonl')).open():
   n=json.loads(line)
   if n['node_id'] in wanted:assert n['node_id'] not in found;found[n['node_id']]=n
  assert set(found)==wanted;nodes[kind]=found
 target={}
 with gzip.open(tr/'target_mask_measurements.tsv.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   if row['node_id'] in nodes['target']:key=row['node_id'],row['mask'];assert key not in target;target[key]=row
 backgrounds={};wantedgenes={(n['gene_a'],n['gene_b']) for n in nodes['background'].values()}
 for row in csv.DictReader((br/'candidate_mask_measurements.tsv').open(),delimiter='\t'):
  if (row['gene_a'],row['gene_b']) in wantedgenes:key=row['gene_a'],row['gene_b'],row['mask'];assert key not in backgrounds;backgrounds[key]=row
 assert len(target)==len(nodes['target'])*2 and len(backgrounds)==len(wantedgenes)*2
 out=Path('results/structural_comparisons/matched-structural-measurements-20260928-v1');out.mkdir(exist_ok=False);counts=Counter();nrows=0;output=out/'unique_pair_mask_measurements.tsv.gz'
 fields=['measurement_disposition','canonical_model_a','canonical_version_a','canonical_model_b','canonical_version_b']
 fitfields=[k for k in next(iter(target.values())) if k.startswith('fit_')]
 with gzip.open(output,'wt',compresslevel=3) as f:
  writer=None
  for key,original in pairs.items():
   tn=nodes['target'][key[0]];bn=nodes['background'][key[1]];assert tn['guide']==bn['guide'] and tn['family']==bn['family']
   for mask in ['full','plddt70']:
    t=target[key[0],mask];b=backgrounds[bn['gene_a'],bn['gene_b'],mask];assert t['pair_key']==tn['pair_key'] and b['pair_key']==bn['pair_key']
    for side in ['a','b']:assert b['model_id_'+side]==bn['model_id_'+side] and int(b['version_'+side])==bn['version_'+side]
    record=dict(original,guide=tn['guide'],family=tn['family'],focal_taxon=tn['taxon_id'],background_taxon_a=bn['taxon_a'],background_taxon_b=bn['taxon_b'],target_gene_node=tn['gene_node'],mask=mask,target_pair_key=tn['pair_key'],background_pair_key=bn['pair_key'],target_sequence_distance=tn['sequence_distance'],background_sequence_distance=bn['sequence_distance'],**flags[key])
    record['target_endpoint_order']=t['target_endpoint_order'];record['background_endpoint_order']=b['candidate_endpoint_order']
    for role,source in [('target',t),('background',b)]:
     for field in fields+fitfields:record[role+'_'+field]=source[field]
    if writer is None:writer=csv.DictWriter(f,fieldnames=list(record),delimiter='\t',lineterminator='\n');writer.writeheader()
    writer.writerow(record);counts[mask+':'+t['measurement_disposition']+':'+b['measurement_disposition']]+=1;nrows+=1
 assert nrows==105350
 for path,h in bindings.items():assert sha(path)==h,path
 result=dict(status='complete_matched_structural_measurement_assembly_pending_independent_readback',script_sha256=sha(__file__),source_hashes=bindings,unique_pairs=len(pairs),pair_mask_rows=nrows,original_selection_uses=sum(selected.values()),screens=ar['screens'],counts=dict(counts),artifacts={output.name:sha(output)},scope='All unique selected pairs with both masks, unchanged original fit fields, per-endpoint screening bitsets, family/gene-node/taxon identities and existing species-contrast pattern IDs. Scenario multiplicities retained as provenance, not statistical weights or independent observations. No outcome averaging, effect estimation, rematching, eligibility threshold change or evolutionary inference.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)
if __name__=='__main__':main()
