#!/usr/bin/env python3
"""Screen every original matching target using audited structural measurements."""
import csv,json
from collections import Counter
from pathlib import Path
from screen_duplication_alignment_reuse import sha
from screen_background_whole_protein_coverage import reasons

def main():
 bindings={}
 def js(path):
  path=Path(path);bindings[str(path)]=sha(path);return json.loads(path.read_text())
 complete=js('metadata/primary_order_summary_completed_20260927.json')
 for path,h in complete['source_hashes'].items():assert sha(path)==h;bindings[path]=h
 graph=Path('results/orthology/background-match-graph-20260927-v1');gr=js(graph/'receipt.json');nodes=graph/'target_nodes.jsonl';assert sha(nodes)==gr['artifacts'][nodes.name];bindings[str(nodes)]=sha(nodes)
 sp=js('metadata/background_control_selection_plan_20260927.json');assert sp['pins'][str(nodes)]==sha(nodes)
 pp=Path('metadata/whole_protein_common_fits_plan_20260927.json');screens=js(pp)['screens']
 root=Path('results/structural_comparisons/primary-usable-orders-20260927-v1');r=js(root/'receipt.json');table=root/'pair_mask_order_summary.tsv';assert sha(table)==r['artifacts'][table.name]
 numeric={}
 for row in csv.DictReader(table.open(),delimiter='\t'):
  key=row['pair_key'],row['mask'];assert key not in numeric;numeric[key]=row
 assert len(numeric)==complete['pair_mask_rows']
 out=Path('results/structural_comparisons/matched-target-coverage-20260928-v1');out.mkdir(exist_ok=False);counts=Counter();same=Counter();seen=set();joint=Counter();n=0
 fields=['target_id','guide','family','gene_a','gene_b','taxon_id','pair_key','same_model','mask']+[s['id']+suffix for s in screens for suffix in ['_pass','_exclusions']]
 with (out/'target_mask_coverage.tsv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
  for line in nodes.open():
   node=json.loads(line);assert node['node_id'] not in seen;seen.add(node['node_id']);passed={}
   for mask in ['full','plddt70']:
    row={k:node[k] for k in fields if k in node};row.update(target_id=node['node_id'],mask=mask)
    measured=not node['same_model'];source=dict(length_a=node['length_a'],length_b=node['length_b'],measurement_disposition='new_model_pair' if measured else 'identical_model_no_alignment')
    if measured:source.update({'fit_'+k:v for k,v in numeric[node['pair_key'],mask].items()})
    else:same[node['guide']+':'+mask]+=1
    passed[mask]={}
    for spec in screens:
     name=spec['id'];why=reasons(source,spec);row[name+'_pass']=int(not why);row[name+'_exclusions']=';'.join(why);passed[mask][name]=not why;counts[node['guide']+':'+mask+':'+name]+=not why
    writer.writerow(row);n+=1
   for spec in screens:joint[node['guide']+':'+spec['id']]+=passed['full'][spec['id']] and passed['plddt70'][spec['id']]
 assert len(seen)==218473 and n==436946
 for path,h in bindings.items():assert sha(path)==h
 result=dict(status='complete_matching_target_coverage_pending_independent_readback',source_hashes=bindings,script_sha256=sha(__file__),screen_helper_sha256=sha('scripts/screen_background_whole_protein_coverage.py'),targets=len(seen),target_mask_rows=n,screens=screens,pass_counts=dict(counts),both_masks_pass_counts=dict(joint),same_model_counts=dict(same),artifacts={'target_mask_coverage.tsv':sha(out/'target_mask_coverage.tsv')},scope='All original218473 modeled matching targets retained under both masks and six identical coverage criteria used for backgrounds. Both orders required, original protein denominator, same-model targets explicitly unmeasured. No rematching, selected-control attrition analysis, calibration or evolutionary inference yet.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','screens']},indent=2))
if __name__=='__main__':main()
