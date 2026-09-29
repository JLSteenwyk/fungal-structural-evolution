#!/usr/bin/env python3
"""Join all matching targets to audited pair metrics without dropping excluded targets."""
import csv,gzip,hashlib,json
from collections import Counter
from pathlib import Path
from screen_duplication_alignment_reuse import sha

def main():
 bindings={}
 def load(path):path=Path(path);bindings[str(path)]=sha(path);return json.loads(path.read_text())
 proof=load('metadata/primary_order_summary_completed_20260927.json');assert proof['status']=='complete_verified_primary_order_summary'
 for path,h in proof['source_hashes'].items():assert sha(path)==h;bindings[path]=h
 root=Path('results/structural_comparisons/primary-usable-orders-20260927-v1');r=load(root/'receipt.json');table=root/'pair_mask_order_summary.tsv';assert sha(table)==r['artifacts'][table.name]
 graph=Path('results/orthology/background-match-graph-20260927-v1');gr=load(graph/'receipt.json');nodes=graph/'target_nodes.jsonl';assert sha(nodes)==gr['artifacts'][nodes.name];bindings[str(nodes)]=sha(nodes)
 selection=load('metadata/background_control_selection_plan_20260927.json');assert selection['pins'][str(nodes)]==sha(nodes)
 measurements={}
 for row in csv.DictReader(table.open(),delimiter='\t'):
  key=row['pair_key'],row['mask'];assert key not in measurements;measurements[key]=row
 assert len(measurements)==proof['pair_mask_rows']==206400
 fields=list(next(iter(measurements.values())));out=Path('results/structural_comparisons/matching-target-measurements-20260928-v1');out.mkdir(exist_ok=False)
 originals={};counts=Counter();used=set();rows=0;path=out/'target_mask_measurements.tsv.gz'
 with gzip.open(path,'wt',compresslevel=3) as f:
  writer=None
  for line in nodes.open():
   n=json.loads(line);nid=n['node_id'];assert nid not in originals;originals[nid]=n
   endpoints=[(n['model_id_'+x],int(n['version_'+x])) for x in ['a','b']];canonical=sorted(endpoints);assert hashlib.sha256(json.dumps(canonical,separators=(',',':')).encode()).hexdigest()==n['pair_key'];assert bool(n['same_model'])==(endpoints[0]==endpoints[1])
   orientation='same_as_canonical' if endpoints==canonical else 'reversed_from_canonical'
   disposition='identical_model_no_alignment' if n['same_model'] else 'measured_pair'
   for mask in ['full','plddt70']:
    summary={k:'' for k in fields} if n['same_model'] else measurements[n['pair_key'],mask]
    if not n['same_model']:used.add(n['pair_key'])
    row=dict(n,mask=mask,measurement_disposition=disposition,target_endpoint_order=orientation,canonical_model_a=canonical[0][0],canonical_version_a=canonical[0][1],canonical_model_b=canonical[1][0],canonical_version_b=canonical[1][1],**{'fit_'+k:v for k,v in summary.items()})
    if writer is None:writer=csv.DictWriter(f,fieldnames=list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
    writer.writerow(row);rows+=1;counts[n['guide']+':'+mask+':'+(summary['order_summary_status'] if not n['same_model'] else disposition)]+=1
   if len(originals)%25000==0:print('Linked matching targets',len(originals),'/218473',flush=True)
 assert len(originals)==218473 and rows==436946 and used=={key[0] for key in measurements}
 seen=set()
 with gzip.open(path,'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   key=row['node_id'],row['mask'];assert key not in seen;seen.add(key);n=originals[key[0]]
   assert all(row[k]==str(v) for k,v in n.items())
   summary={k:'' for k in fields} if n['same_model'] else measurements[n['pair_key'],key[1]]
   assert all(row['fit_'+k]==v for k,v in summary.items())
 assert seen=={(nid,mask) for nid in originals for mask in ['full','plddt70']}
 for source,h in bindings.items():assert sha(source)==h,source
 result=dict(status='complete_matching_target_measurement_links_pending_independent_readback',script_sha256=sha(__file__),source_hashes=bindings,targets=len(originals),target_mask_rows=rows,distinct_measured_pairs=len(used),counts=dict(counts),artifacts={path.name:sha(path)},scope='All matching target node fields and both masks preserved. Fits retain canonical model endpoints and both native input orders, with gene-to-canonical orientation explicit. Identical model pairs have blank fits rather than imputed zero; numerical exclusions and missing inputs retained. No rematching, structural-effect estimate or phylogenetic correction. Frozen original matching atlas, not expanded September28 targets.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)
if __name__=='__main__':main()
