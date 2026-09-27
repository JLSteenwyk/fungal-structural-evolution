#!/usr/bin/env python3
"""Check every matched-domain metric and complete four-order grid from audited sources."""
import argparse,csv,gzip,json,hashlib,math,subprocess
from pathlib import Path
from collections import Counter


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    for p,h in plan['pins'].items():assert sha(p)==h,p
    root=Path(plan['output']);r=json.loads((root/'receipt.json').read_text());assert r['status']=='complete_matched_domain_measurements_pending_independent_readback' and r['plan_sha256']==ph
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    inventory=Path(plan['inventory']);qualified=Path(plan['qualification'])
    assert sha(inventory/'receipt.json')==r['inventory_receipt_sha256'] and sha(qualified/'receipt.json')==r['qualification_receipt_sha256']
    expected={};seen={}
    for line in (inventory/'domain_configurations.jsonl').open():
        c=json.loads(line)
        for m in c['matches']:
            key=c['domain_config_id'],m['pfam_accession'],m['boundary'];assert key not in expected
            expected[key]=(m['target_domain_pair'],m['background_domain_pair'],c['status']);seen[key]=0
    sources={}
    for side,spec in plan['cohorts'].items():
        tables=[]
        for path,ordered in [(Path(spec['diagnostic'])/'numeric_readback.tsv',True),(Path(spec['coverage'])/'pair_mask_coverage.tsv',False),(qualified/(side+'.tsv'),False)]:
            data={}
            for x in csv.DictReader(path.open(),delimiter='\t'):
                key=(x['pair_key'],x['mask'],int(x['order'])) if ordered else (x['pair_key'],x['mask']);assert key not in data;data[key]=x
            tables.append(data)
        sources[side]=tables
    count=0;values_checked=0;flags_checked=0;counts=Counter()
    metrics=['aligned_length','rmsd_recomputed','sequence_identity_exact','tm_a','tm_b','retained_coverage_a','retained_coverage_b','original_coverage_a','original_coverage_b','joint_plddt70_fraction']
    with gzip.open(root/'matched_domain_contrasts.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=row['domain_config_id'],row['pfam_accession'],row['boundary'];ta,ba,configuration=expected[key];mask=row['mask'];to=int(row['target_order']);bo=int(row['background_order'])
            assert mask in ['full','plddt70'] and to in [0,1] and bo in [0,1]
            bit=1<<((0 if mask=='full' else 4)+2*to+bo);assert not seen[key]&bit;seen[key]|=bit
            assert (row['target_domain_pair'],row['background_domain_pair'])==(ta,ba)
            usable={};numeric={}
            for side,pair,order in [('target',ta,to),('background',ba,bo)]:
                ns,cs,qs=sources[side];n=ns.get((pair,mask,order));c=cs[pair,mask];q=qs[pair,mask]
                assert row[side+'_geometry_status']==q[f'order{order}_geometry_status']
                assert row[side+'_rmsd_status']==(n['rmsd_status'] if n else 'input_unavailable')
                ok=n is not None and float(n['rmsd_rounding_error'])<=.00501 and int(n['aligned_length'])>=3 and row[side+'_geometry_status']=='unique_at_numeric_tolerance';usable[side]=ok
                vals={}
                if ok:
                    for col in ['aligned_length','rmsd_recomputed','sequence_identity_exact','joint_plddt70_fraction']:vals[col]=float(n[col])
                    endpoints=['a','b'] if order==0 else ['b','a']
                    for native,end in zip(['left','right'],endpoints):
                        vals['tm_'+end]=float(n['tm_'+native+'_native']);vals['retained_coverage_'+end]=float(n['coverage_'+native])
                    vals['original_coverage_a']=int(n['aligned_length'])/int(c['length_a']);vals['original_coverage_b']=int(n['aligned_length'])/int(c['length_b'])
                    numeric[side]=vals
                for col in metrics:
                    value=row[side+'_'+col]
                    if ok:assert math.isfinite(float(value)) and float(value)==vals[col];values_checked+=1
                    else:assert value==''
            ready=all(usable.values()) and configuration=='shared_domain_comparisons'
            assert row['contrast_status']==('numerically_computable' if ready else 'excluded_numerically_or_identical')
            for metric,name in [('rmsd_recomputed','rmsd'),('sequence_identity_exact','sequence_identity'),('aligned_length','aligned_length')]:
                value=row[name+'_target_minus_background']
                if ready:assert float(value)==numeric['target'][metric]-numeric['background'][metric];values_checked+=1
                else:assert value==''
            for screen in plan['screens']:
                statuses={}
                for mm in ['full','plddt70']:
                    statuses[mm]=sources['target'][2][ta,mm][screen+'_pass']=='1' and sources['background'][2][ba,mm][screen+'_pass']=='1' and configuration=='shared_domain_comparisons'
                assert row[screen+'_mask_pass']==str(int(statuses[mask])) and row[screen+'_both_masks_pass']==str(int(statuses['full'] and statuses['plddt70']));flags_checked+=2
                assert not statuses[mask] or ready
            counts[mask+':'+row['contrast_status']]+=1;count+=1
            if count%200000==0:print('Verified matched measurements',count,'/',r['rows'],flush=True)
    assert all(v==255 for v in seen.values()) and count==len(expected)*8==r['rows'] and len(expected)==r['shared_boundary_matches']
    assert dict(counts)==r['counts']
    assert sha(a.plan)==ph
    for p,h in plan['pins'].items():assert sha(p)==h,p
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    result=dict(status='passed_full_matched_domain_measurement_readback',rows=count,shared_boundary_matches=len(expected),numeric_values_checked=values_checked,qualification_flags_checked=flags_checked,counts=dict(counts),source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),scope='All eight mask/order combinations per shared boundary match verified. Every measurement, endpoint orientation, original-length coverage, blank excluded field, numerical status, difference and qualification flag reconstructed without producer helpers. No new coordinate fits, common-residue four-protein fits or biological inference.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
