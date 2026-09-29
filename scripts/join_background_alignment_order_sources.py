#!/usr/bin/env python3
"""Combine audited background/reference measurements for the full control inventory."""
import csv,json
from collections import Counter
from pathlib import Path
from screen_duplication_alignment_reuse import sha

def main():
 bindings={}
 def bind(p):
  p=Path(p);bindings[str(p)]=sha(p);return p
 def js(p):return json.loads(bind(p).read_text())
 def table(p):return list(csv.DictReader(bind(p).open(),delimiter='\t'))
 inventory=Path('results/structural_comparisons/background-measurement-inventory-20260927-v1')
 ir=js(inventory/'receipt.json');proof=js('metadata/background_measurement_inventory_completed_readback_20260927.json')
 assert proof['producer_receipt_sha256']==sha(inventory/'receipt.json')
 for name,h in ir['artifacts'].items():assert sha(bind(inventory/name))==h
 pairs=table(inventory/'model_pairs.tsv');assert len(pairs)==ir['distinct_eligible_model_pairs']==71461
 classes=Counter(x['work_disposition'] for x in pairs);assert classes==dict(new_model_pair=71450,already_in_reference_queue=11)
 plans={};summaries={};receipts={}
 for kind,pp,ap in [('new','metadata/background_alignment_usable_orders_plan_20260928.json','metadata/background_alignment_usable_orders_readback_20260928.json'),('reference','metadata/duplication_reference_usable_orders_plan_20260927.json','metadata/duplication_reference_usable_orders_readback_20260927.json')]:
  plan=js(pp);audit=js(ap)
  for path,h in plan['pins'].items():assert sha(bind(path))==h
  root=Path(plan['output']);r=js(root/'receipt.json');assert r['plan_sha256']==sha(pp) and audit['source_receipt_sha256']==sha(root/'receipt.json')
  assert audit['status']==('passed_full_background_alignment_usable_order_summary_readback' if kind=='new' else 'passed_full_reference_usable_order_summary_readback')
  assert audit['geometry_audit_sha256']==sha(bind(plan['geometry_audit']))
  path=root/'pair_mask_order_summary.tsv';assert sha(bind(path))==r['artifacts'][path.name]
  rows=table(path);index={(x['pair_key'],x['mask']):x for x in rows};assert len(index)==len(rows)==audit['pair_mask_rows']
  summaries[kind]=index;plans[kind]=plan;receipts[kind]=sha(root/'receipt.json')
 bp=js('metadata/background_alignment_plan_20260927.json');rp=js('metadata/duplication_reference_alignment_plan_20260926.json')
 assert bp['usalign']==rp['usalign'] and bp['options']==rp['options']
 # Both recorded plans pin the same production executable, and its bytes still match.
 binary=bp['usalign'];assert bp['pins'][binary]==rp['pins'][binary]==sha(bind(binary))
 reused={x['pair_key']:x for x in pairs if x['work_disposition']=='already_in_reference_queue'}
 inputs={}
 for spec in bp['input_sources'].values():
  folder=Path(spec['inputs']);receipt=js(folder/'receipt.json');manifest=bind(folder/'inputs.jsonl');assert sha(manifest)==receipt['artifacts']['inputs.jsonl']
  for line in manifest.open():
   row=json.loads(line);key=(row['model_id'],row['version'],row['mask']);assert key not in inputs;inputs[key]=row
 native=Path(plans['reference']['native']);nr=js(native/'receipt.json');cm=native/'checkpoint_manifest.tsv';assert sha(bind(cm))==nr['artifacts'][cm.name]
 hashes={x['path']:x['sha256'] for x in table(cm)};reuse_checks=[]
 for key,row in reused.items():
  ends=[(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))]
  for mask in ['full','plddt70']:
   for order in [0,1]:
    rel=f'pairs/{key[:2]}/{key}-{mask}-{order}.json';path=bind(native/rel);assert sha(path)==hashes[rel];record=js(path)
    directed=ends if order==0 else ends[::-1];current=[inputs[(*e,mask)] for e in directed]
    identities=[dict(model_id=e[0],version=e[1],status=v['status'],path=v.get('path'),sha256=v.get('sha256')) for e,v in zip(directed,current)]
    assert record['inputs']==identities and (record['pair_key'],record['mask'],record['order'])==(key,mask,order)
    ready=all(v['status']=='ready' for v in current)
    assert record['command']==([binary,*[v['path'] for v in current],*bp['options']] if ready else None)
    for v in current:
     if v['status']=='ready':assert sha(bind(v['path']))==v['sha256']
    reuse_checks.append(dict(pair_key=key,mask=mask,order=order,checkpoint=str(path),sha256=sha(path)))
 out=Path('results/structural_comparisons/background-combined-orders-20260928-v1');out.mkdir(exist_ok=False)
 combined=[];seen=set();counts=Counter()
 for pair in pairs:
  kind='new' if pair['work_disposition']=='new_model_pair' else 'reference'
  for mask in ['full','plddt70']:
   key=pair['pair_key'],mask;assert key not in seen;seen.add(key)
   row=dict(summaries[kind][key],measurement_source=kind,source_summary_receipt_sha256=receipts[kind])
   combined.append(row);counts[kind+':'+mask+':'+row['order_summary_status']]+=1
 assert len(combined)==142922
 with (out/'pair_mask_order_summary.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(combined[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(combined)
 # Exact readback: every source field, including blanks and reasons, must survive.
 read=table(out/'pair_mask_order_summary.tsv');assert len(read)==len(combined)
 for row in read:
  original=summaries[row['measurement_source']][row['pair_key'],row['mask']]
  assert all(row[k]==v for k,v in original.items())
 for path,h in bindings.items():assert sha(path)==h
 result=dict(status='complete_combined_background_order_sources',script_sha256=sha(__file__),source_hashes=bindings,distinct_pairs=len(pairs),pair_mask_rows=len(read),source_pair_counts=dict(classes),counts=dict(counts),reused_checkpoint_checks=reuse_checks,artifacts={'pair_mask_order_summary.tsv':sha(out/'pair_mask_order_summary.tsv')},scope='All71461 distinct eligible background model pairs, including11 previously measured reference pairs, both masks and orders. Exact source fields retained and serialized readback checked. Reused44 checkpoints match frozen endpoint identities, input PDB bytes, masks, executable and options. Candidate/event linkage, coverage/confidence qualification and evolutionary tests remain downstream; no discrepancy cleared.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','reused_checkpoint_checks']},indent=2))
if __name__=='__main__':main()
