#!/usr/bin/env python3
"""Reconstruct every whole-protein contrast from untransformed matched measurements."""
import csv,gzip,json,math,subprocess
from collections import Counter
from pathlib import Path
from screen_duplication_alignment_reuse import sha

def main():
 unit='fungal-matched-whole-protein-contrasts-20260928.service';state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 root=Path('results/structural_comparisons/matched-whole-protein-contrasts-20260928-v1');r=json.loads((root/'receipt.json').read_text());rh=sha(root/'receipt.json');assert r['status']=='complete_matched_whole_protein_contrasts_pending_independent_readback';bindings=dict(r['source_hashes'])
 for n,h in r['artifacts'].items():bindings[str(root/n)]=h
 for path,h in bindings.items():assert sha(path)==h,path
 lengths={}
 for kind in ['target','background']:
  for line in Path('results/orthology/background-match-graph-20260927-v1',kind+'_nodes.jsonl').open():
   node=json.loads(line)
   for side in ['a','b']:
    key=node['model_id_'+side],node['version_'+side];v=node['length_'+side];assert key not in lengths or lengths[key]==v;lengths[key]=v
 source={}
 with gzip.open('results/structural_comparisons/matched-structural-measurements-20260928-v1/unique_pair_mask_measurements.tsv.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   key=row['target_id'],row['background_id'],row['mask'];assert key not in source;source[key]=row
 identity=['target_id','background_id','selected_record_uses','guide','family','focal_taxon','background_taxon_a','background_taxon_b','target_gene_node','species_pattern_id','target_pair_key','background_pair_key','target_sequence_distance','background_sequence_distance','mask']
 metrics=['aligned_length','rmsd_recomputed','sequence_identity_exact','tm_a','tm_b','joint_plddt70_fraction'];counts=Counter();seen=set();numeric=0
 with gzip.open(root/'whole_protein_contrasts.tsv.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   key=tuple(row[k] for k in ['target_id','background_id','mask']);s=source[key];to,bo=int(row['target_order']),int(row['background_order']);assert to in [0,1] and bo in [0,1] and (*key,to,bo) not in seen;seen.add((*key,to,bo))
   expected={k:s[k] for k in identity};expected.update(target_order=str(to),background_order=str(bo));ok=[]
   for role,order in [('target',to),('background',bo)]:
    prefix=f'{role}_fit_order{order}_';status=s[prefix+'status'] if s[prefix+'status'] else s[role+'_measurement_disposition'];expected[role+'_status']=status;usable=status=='aligned';ok.append(usable)
    for metric in metrics:expected[role+'_'+metric]=float(s[prefix+metric]) if usable else ''
    cov=[]
    for side in ['a','b']:
     model=s[role+'_canonical_model_'+side];version=s[role+'_canonical_version_'+side]
     value=float(s[prefix+'aligned_length'])/lengths[model,int(version)] if usable else '';expected[role+'_original_coverage_'+side]=value;cov.append(value)
    expected[role+'_minimum_original_coverage']=min(cov) if usable else ''
    expected[role+'_mean_endpoint_tm']=sum(float(s[prefix+'tm_'+side]) for side in ['a','b'])/2 if usable else ''
   ready=all(ok);expected['contrast_status']='numerically_computable' if ready else 'excluded_or_unmeasured'
   for field,metric in [('rmsd_difference','rmsd_recomputed'),('identity_difference','sequence_identity_exact'),('original_coverage_difference','minimum_original_coverage'),('confidence_fraction_difference','joint_plddt70_fraction')]:expected[field]=expected['target_'+metric]-expected['background_'+metric] if ready else ''
   expected['log_aligned_length_ratio']=math.log(expected['target_aligned_length'])-math.log(expected['background_aligned_length']) if ready else ''
   expected['mean_endpoint_tm_divergence_difference']=(1-expected['target_mean_endpoint_tm'])-(1-expected['background_mean_endpoint_tm']) if ready else ''
   t=float(s['target_sequence_distance']);b=float(s['background_sequence_distance']);expected['sequence_distance_difference']=t-b;expected['sequence_distance_status']='both_positive' if t>0 and b>0 else 'zero_distance_present';expected['log_positive_sequence_distance_difference']=math.log(t)-math.log(b) if t>0 and b>0 else ''
   for i,spec in enumerate(r['screens']):
    for cohort,mask in [('mask',s['mask']),('both_masks','both_masks')]:
     passed=all(int(s[role+'_'+mask+'_pass_bits'])//2**i%2==1 for role in ['target','background']);assert not passed or ready;expected[spec['id']+'_'+cohort+'_pass']=str(int(passed))
   assert set(row)==set(expected)
   for field,value in expected.items():
    if isinstance(value,float):assert math.isfinite(value) and math.isclose(float(row[field]),value,abs_tol=1e-12,rel_tol=1e-10),(key,to,bo,field,row[field],value);numeric+=1
    else:assert row[field]==value,(key,field)
   counts[s['mask']+':'+expected['contrast_status']]+=1
   if len(seen)%100000==0:print('Verified whole protein contrast rows',len(seen),'/421400',flush=True)
 assert seen=={(*k,to,bo) for k in source for to in [0,1] for bo in [0,1]} and len(seen)==r['rows']==421400 and dict(counts)==r['counts']
 for path,h in bindings.items():assert sha(path)==h,path
 assert sha(root/'receipt.json')==rh
 result=dict(status='passed_full_matched_whole_protein_contrast_readback',source_receipt_sha256=rh,checker_sha256=sha(__file__),producer_terminal_state=state,rows=len(seen),numeric_values_checked=numeric,counts=dict(counts),scope='Every identity, status, raw metric, endpoint coverage and derived contrast checked; original protein lengths independently joined by canonical model/version, not producer orientation logic. All masks, four order combinations, screen flags and zero-distance exclusions retained. No effect inference, calibration or independent-replicate claim.')
 with Path('metadata/matched_whole_protein_contrasts_completed_readback_20260928.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result),flush=True)
if __name__=='__main__':main()
