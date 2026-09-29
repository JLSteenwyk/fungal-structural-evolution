#!/usr/bin/env python3
"""Prepare whole-protein matched contrasts without selecting an alignment order or screen."""
import csv,gzip,json,math,subprocess
from collections import Counter
from pathlib import Path
from screen_duplication_alignment_reuse import sha

METRICS=['aligned_length','rmsd_recomputed','sequence_identity_exact','tm_a','tm_b','joint_plddt70_fraction']
def contrast(row,lengths,to,bo,screens):
 out={};usable=[];measurements={}
 for role,order in [('target',to),('background',bo)]:
  prefix=f'{role}_fit_order{order}_';status=row[prefix+'status'] or row[role+'_measurement_disposition'];out[role+'_status']=status;ok=status=='aligned';usable.append(ok)
  for metric in METRICS:out[role+'_'+metric]=float(row[prefix+metric]) if ok else ''
  n=out[role+'_aligned_length'];a,b=lengths[role]
  for side,length in [('a',a),('b',b)]:out[role+'_original_coverage_'+side]=n/length if ok else ''
  out[role+'_minimum_original_coverage']=n/max(a,b) if ok else ''
  out[role+'_mean_endpoint_tm']=(out[role+'_tm_a']+out[role+'_tm_b'])/2 if ok else ''
 ready=all(usable);out['contrast_status']='numerically_computable' if ready else 'excluded_or_unmeasured'
 fields=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference','mean_endpoint_tm_divergence_difference']
 out.update({k:'' for k in fields})
 if ready:
  out.update(rmsd_difference=out['target_rmsd_recomputed']-out['background_rmsd_recomputed'],identity_difference=out['target_sequence_identity_exact']-out['background_sequence_identity_exact'],original_coverage_difference=out['target_minimum_original_coverage']-out['background_minimum_original_coverage'],log_aligned_length_ratio=math.log(out['target_aligned_length']/out['background_aligned_length']),confidence_fraction_difference=out['target_joint_plddt70_fraction']-out['background_joint_plddt70_fraction'],mean_endpoint_tm_divergence_difference=out['background_mean_endpoint_tm']-out['target_mean_endpoint_tm'])
 t=float(row['target_sequence_distance']);b=float(row['background_sequence_distance']);assert t>=0 and b>=0
 out.update(sequence_distance_difference=t-b,log_positive_sequence_distance_difference=math.log(t/b) if t>0 and b>0 else '',sequence_distance_status='both_positive' if t>0 and b>0 else 'zero_distance_present')
 for i,spec in enumerate(screens):
  for cohort,mask in [('mask',row['mask']),('both_masks','both_masks')]:
   passed=all(int(row[role+'_'+mask+'_pass_bits'])&(1<<i) for role in ['target','background']);assert not passed or ready;out[spec['id']+'_'+cohort+'_pass']=int(passed)
 return out

def main():
 root=Path('results/structural_comparisons/matched-structural-measurements-20260928-v1');r=json.loads((root/'receipt.json').read_text());proof=Path('metadata/matched_structural_measurements_completed_readback_20260928.json');a=json.loads(proof.read_text());assert a['status']=='passed_full_matched_structural_measurement_readback' and a['source_receipt_sha256']==sha(root/'receipt.json')
 unit='fungal-matched-structural-measurement-readback-20260928.service';state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 bindings={str(root/'receipt.json'):sha(root/'receipt.json'),str(proof):sha(proof)}
 for name,h in r['artifacts'].items():assert sha(root/name)==h;bindings[str(root/name)]=h
 lengths={}
 for role in ['target','background']:
  p=Path('results/orthology/background-match-graph-20260927-v1')/(role+'_nodes.jsonl');assert sha(p)==r['source_hashes'][str(p)];bindings[str(p)]=sha(p);lengths[role]={}
  for line in p.open():
   n=json.loads(line);assert n['node_id'] not in lengths[role];lengths[role][n['node_id']]=(n['length_a'],n['length_b'])
 out=Path('results/structural_comparisons/matched-whole-protein-contrasts-20260928-v1');out.mkdir(exist_ok=False);path=out/'whole_protein_contrasts.tsv.gz';counts=Counter();total=0
 identity=['target_id','background_id','selected_record_uses','guide','family','focal_taxon','background_taxon_a','background_taxon_b','target_gene_node','species_pattern_id','target_pair_key','background_pair_key','target_sequence_distance','background_sequence_distance','mask']
 with gzip.open(root/'unique_pair_mask_measurements.tsv.gz','rt') as source,gzip.open(path,'wt',compresslevel=3) as dest:
  writer=None
  for row in csv.DictReader(source,delimiter='\t'):
   sizes={role:lengths[role][row[role+'_id']] for role in ['target','background']}
   for role in ['target','background']:
    if row[role+'_endpoint_order']=='reversed_from_canonical':sizes[role]=sizes[role][::-1]
   for to in [0,1]:
    for bo in [0,1]:
     record={k:row[k] for k in identity};record.update(target_order=to,background_order=bo,**contrast(row,sizes,to,bo,r['screens']))
     if writer is None:writer=csv.DictWriter(dest,fieldnames=list(record),delimiter='\t',lineterminator='\n');writer.writeheader()
     writer.writerow(record);total+=1;counts[row['mask']+':'+record['contrast_status']]+=1
   if total%100000==0:print('Prepared whole protein contrast rows',total,'/421400',flush=True)
 assert total==r['pair_mask_rows']*4==421400
 for p,h in bindings.items():assert sha(p)==h,p
 receipt=dict(status='complete_matched_whole_protein_contrasts_pending_independent_readback',script_sha256=sha(__file__),source_hashes=bindings,rows=total,unique_pairs=r['unique_pairs'],counts=dict(counts),screens=r['screens'],artifacts={path.name:sha(path)},scope='All unique matched pairs x2 masks x4 input-order combinations. Original full-protein coverage denominators; raw and log-positive gene-tree sequence-distance differences with explicit zero cases; unchanged six screens. RMSD contrast compares separately fitted cores; mean-endpoint TM divergence is a descriptive symmetric score, not physical displacement. No preferred order, rematching, independent replication, phylogenetic fitting, calibration or effect inference.')
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
