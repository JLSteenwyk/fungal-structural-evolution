#!/usr/bin/env python3
"""Link complete audited background measurements without dropping unmeasured candidates."""
import csv,json,hashlib
from collections import Counter
from pathlib import Path
from screen_duplication_alignment_reuse import sha

def main():
 sources={}
 def bind(path):
  path=Path(path);sources[str(path)]=sha(path);return path
 def js(path):return json.loads(bind(path).read_text())
 def rows(path):return csv.DictReader(bind(path).open(),delimiter='\t')
 union=Path('results/structural_comparisons/background-combined-orders-20260928-v1')
 completion=js('metadata/background_combined_orders_completed_20260928.json');ur=js(union/'receipt.json');assert completion['receipt_sha256']==sha(union/'receipt.json')
 assert sha(bind(union/'pair_mask_order_summary.tsv'))==ur['artifacts']['pair_mask_order_summary.tsv']
 inventory=Path('results/structural_comparisons/background-measurement-inventory-20260927-v1');ir=js(inventory/'receipt.json')
 assert ur['source_hashes'][str(inventory/'receipt.json')]==sha(inventory/'receipt.json')
 for name in ['model_pairs.tsv','candidate_measurement_links.tsv']:assert sha(bind(inventory/name))==ir['artifacts'][name]
 pairs={r['pair_key']:r for r in rows(inventory/'model_pairs.tsv')};assert len(pairs)==71461
 measurements={}
 for row in rows(union/'pair_mask_order_summary.tsv'):
  key=row['pair_key'],row['mask'];assert key not in measurements;measurements[key]=row
 assert set(measurements)=={(key,mask) for key in pairs for mask in ['full','plddt70']}
 candidates={}
 for row in rows(inventory/'candidate_measurement_links.tsv'):
  key=row['native_pair_key'];assert key not in candidates;candidates[key]=row
 assert len(candidates)==ir['candidates']==78372
 fields=list(next(iter(measurements.values())));out=Path('results/structural_comparisons/background-candidate-measurements-20260928-v1');out.mkdir(exist_ok=False)
 counts=Counter();orientations=Counter();used=set();n=0
 with (out/'candidate_mask_measurements.tsv').open('w') as f:
  writer=None
  for key,row in candidates.items():
   disposition=row['measurement_disposition'];measured=disposition in ['new_model_pair','already_in_reference_queue'];pair=row['pair_key']
   orientation='not_measured';canonical={k:'' for k in ['canonical_model_a','canonical_version_a','canonical_model_b','canonical_version_b']}
   if measured:
    p=pairs[pair];assert p['work_disposition']==disposition;used.add(pair)
    ends=[(p['model_a'],int(p['version_a'])),(p['model_b'],int(p['version_b']))];cand=[(row['model_id_a'],int(row['version_a'])),(row['model_id_b'],int(row['version_b']))]
    assert cand==ends or cand==ends[::-1];orientation='same_as_canonical' if cand==ends else 'reversed_from_canonical'
    canonical=dict(canonical_model_a=ends[0][0],canonical_version_a=str(ends[0][1]),canonical_model_b=ends[1][0],canonical_version_b=str(ends[1][1]))
   else:assert disposition in ['identical_model_no_alignment','neither_guide_eligible']
   for mask in ['full','plddt70']:
    summary=measurements[pair,mask] if measured else {k:'' for k in fields}
    record=dict(row,mask=mask,candidate_endpoint_order=orientation,**canonical,**{'fit_'+k:v for k,v in summary.items()})
    if writer is None:writer=csv.DictWriter(f,fieldnames=list(record),delimiter='\t',lineterminator='\n');writer.writeheader()
    writer.writerow(record);counts[mask+':'+disposition+':'+(summary['order_summary_status'] if measured else 'not_measured')]+=1;orientations[orientation]+=1;n+=1
 assert used==set(pairs) and n==156744
 seen=set()
 for row in rows(out/'candidate_mask_measurements.tsv'):
  key=row['native_pair_key'],row['mask'];assert key not in seen;seen.add(key);original=candidates[key[0]]
  assert all(row[k]==v for k,v in original.items())
  measured=original['measurement_disposition'] in ['new_model_pair','already_in_reference_queue']
  expected=measurements[original['pair_key'],key[1]] if measured else {k:'' for k in fields}
  assert all(row['fit_'+k]==v for k,v in expected.items())
 assert seen=={(key,mask) for key in candidates for mask in ['full','plddt70']}
 for path,h in sources.items():assert sha(path)==h
 result=dict(status='complete_all_background_candidate_measurement_links',script_sha256=sha(__file__),source_hashes=sources,candidates=len(candidates),candidate_mask_rows=n,distinct_measured_pairs=len(used),counts=dict(counts),endpoint_orientation_counts=dict(orientations),artifacts={'candidate_mask_measurements.tsv':sha(out/'candidate_mask_measurements.tsv')},scope='Every78372 inventory candidate retained under both masks with all original guide/native-orthology flags, sequence distances, length/confidence covariates and measurement dispositions. Measurements use explicitly named canonical endpoints; candidate orientation recorded. Identical-model and neither-guide cases retain blank fits, never imputed zero. All source fields and fit fields checked after serialization. Numerical usability is not biological qualification; coverage/confidence, phylogenetic matching and evolutionary tests remain required.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))
if __name__=='__main__':main()
