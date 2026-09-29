#!/usr/bin/env python3
"""Bind all matched whole-protein pairs to verified family and species covariance indices."""
import csv,gzip,hashlib,json
from pathlib import Path
import numpy as np
from screen_duplication_alignment_reuse import sha

def main():
 bindings={}
 def load(path):path=Path(path);bindings[str(path)]=sha(path);return json.loads(path.read_text())
 def root(path):
  path=Path(path);r=load(path/'receipt.json')
  for name,h in r['artifacts'].items():assert sha(path/name)==h;bindings[str(path/name)]=h
  return path,r
 nodes,nr=root('results/orthology/selected-taxon-inputs-20260927-v1');a=load('metadata/selected_taxon_inputs_readback_20260927.json');assert a['status']=='passed_full_selected_taxon_input_readback' and a['producer_receipt_sha256']==sha(nodes/'receipt.json')
 factors,fr=root('results/phylogeny/matched-species-covariance-factors-20260927-v1');fa=load('metadata/matched_species_factor_readback_20260927.json');assert fa['status']=='passed_full_matched_species_factor_tree_readback' and fa['source_receipt_sha256']==sha(factors/'receipt.json')
 species,sr=root('results/phylogeny/matched-species-contrasts-20260927-v1');assert fr['pins'][str(species/'receipt.json')]==sha(species/'receipt.json')
 measured,mr=root('results/structural_comparisons/matched-structural-measurements-20260928-v1');ma=load('metadata/matched_structural_measurements_completed_readback_20260928.json');assert ma['source_receipt_sha256']==sha(measured/'receipt.json')
 node={(n['role'],n['node_id']):n for n in csv.DictReader((nodes/'selected_node_taxa.tsv').open(),delimiter='\t')};assert len(node)==a['selected_nodes']
 patterns={n['species_pattern_id']:int(n['row_index']) for n in csv.DictReader((factors/'patterns.tsv').open(),delimiter='\t')};assert len(patterns)==fr['patterns']==4568
 quadratic={(n['tree'],n['species_pattern_id']):float(n['kernel_quadratic_form']) for n in csv.DictReader((species/'kernel_quadratic_forms.tsv').open(),delimiter='\t')};diagonal={}
 for tree in sorted(fr['trees']):
  with np.load(factors/(tree+'.npz')) as z:
   factor=z['factor'];assert factor.shape==(4568,242);diagonal[tree]=np.einsum('ij,ij->i',factor,factor)
  for pid,i in patterns.items():assert np.isclose(diagonal[tree][i],quadratic[tree,pid],rtol=1e-10,atol=1e-12)
 out=Path('results/structural_comparisons/whole-protein-covariance-index-20260928-v1');out.mkdir(exist_ok=False);table=out/'pair_covariance_index.tsv';seen={};family_for_background={}
 with gzip.open(measured/'unique_pair_mask_measurements.tsv.gz','rt') as source,table.open('w') as dest:
  writer=None
  for row in csv.DictReader(source,delimiter='\t'):
   key=row['target_id'],row['background_id'];t=node['target',key[0]];b=node['background',key[1]]
   assert t['guide']==b['guide']==row['guide'] and t['family']==b['family']==row['family'] and t['family_component']==b['family_component']
   pid=row['species_pattern_id'];i=patterns[pid];assert key[1] not in family_for_background or family_for_background[key[1]]==t['family_component'];family_for_background[key[1]]=t['family_component']
   identifier=hashlib.sha256(json.dumps([*key,t['family_component'],pid],separators=(',',':')).encode()).hexdigest()
   result=dict(target_id=key[0],background_id=key[1],guide=t['guide'],family=t['family'],family_component=t['family_component'],species_pattern_id=pid,species_pattern_row=i,row_identity=identifier)
   result.update({tree+'_kernel_quadratic_form':float(diagonal[tree][i]) for tree in diagonal})
   if key in seen:assert seen[key]==result;continue
   seen[key]=result
   if writer is None:writer=csv.DictWriter(dest,list(result),delimiter='\t',lineterminator='\n');writer.writeheader()
   writer.writerow(result)
 assert len(seen)==mr['unique_pairs']==52675
 checked=0;serialized=set()
 for row in csv.DictReader(table.open(),delimiter='\t'):
  key=row['target_id'],row['background_id'];assert key not in serialized;serialized.add(key);assert row=={k:str(v) for k,v in seen[key].items()};checked+=len(row)
 assert serialized==set(seen)
 for path,h in bindings.items():assert sha(path)==h,path
 receipt=dict(status='complete_whole_protein_covariance_index_with_full_serialized_readback',script_sha256=sha(__file__),source_hashes=bindings,pairs=len(seen),serialized_fields_checked=checked,trees=list(diagonal),patterns=len(patterns),factor_rank=242,artifacts={table.name:sha(table)},scope='All unique original matched pairs bound to audited shared-entity family components and exact species-factor rows; all five kernel diagonals checked, full exported fields read back. Background nesting verified. Covariance remains a working nuisance model, not an estimated process or complete dependence correction; future per-setting filtering and fit calibration required.')
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='source_hashes'}),flush=True)
if __name__=='__main__':main()
