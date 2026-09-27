#!/usr/bin/env python3
"""Preserve matched-domain measurements for all masks and input-order combinations."""
import argparse,csv,gzip,json
from pathlib import Path
from collections import Counter
from screen_duplication_domain_alignment_coverage import sha

METRICS=['aligned_length','rmsd_recomputed','sequence_identity_exact','tm_a','tm_b','retained_coverage_a','retained_coverage_b','original_coverage_a','original_coverage_b','joint_plddt70_fraction']

def measurement(n,c,order):
    if n is None:return None
    row={k:float(n[k]) for k in ['aligned_length','rmsd_recomputed','sequence_identity_exact','joint_plddt70_fraction']}
    for end,other in [('a','b'),('b','a')]:
        side=('left' if end=='a' else 'right') if order==0 else ('right' if end=='a' else 'left')
        row['tm_'+end]=float(n['tm_'+side+'_native']);row['retained_coverage_'+end]=float(n['coverage_'+side]);row['original_coverage_'+end]=float(c[f'order{order}_original_coverage_{end}'])
    return row

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();inventory=Path(plan['inventory']);ir=json.loads((inventory/'receipt.json').read_text());ia=json.loads(Path(plan['inventory_audit']).read_text());assert ia['status']=='passed_full_selected_control_domain_inventory_readback' and ia['source_receipt_sha256']==sha(inventory/'receipt.json')
    assert sha(inventory/'domain_configurations.jsonl')==ir['artifacts']['domain_configurations.jsonl']
    qualified=Path(plan['qualification']);qr=json.loads((qualified/'receipt.json').read_text());qa=json.loads(Path(plan['qualification_audit']).read_text());assert qa['status']=='passed_full_domain_coverage_geometry_readback' and qa['source_receipt_sha256']==sha(qualified/'receipt.json')
    sources={}
    for label,spec in plan['cohorts'].items():
        numeric=Path(spec['diagnostic']);coverage=Path(spec['coverage']);nr=json.loads((numeric/'receipt.json').read_text());cr=json.loads((coverage/'receipt.json').read_text())
        assert cr['diagnostic_receipt_sha256']==sha(numeric/'receipt.json')
        assert qr['source_bindings'][label][str(coverage/'receipt.json')]==sha(coverage/'receipt.json')
        for root,receipt in [(numeric,nr),(coverage,cr),(qualified,qr)]:
            for name,h in receipt['artifacts'].items():assert sha(root/name)==h
        n={};c={};q={}
        for index,path,order in [(n,numeric/'numeric_readback.tsv',True),(c,coverage/'pair_mask_coverage.tsv',False),(q,qualified/(label+'.tsv'),False)]:
            for row in csv.DictReader(path.open(),delimiter='\t'):
                key=(row['pair_key'],row['mask'],int(row['order'])) if order else (row['pair_key'],row['mask']);assert key not in index;index[key]=row
        assert set(c)==set(q);sources[label]=(n,c,q)
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);counts=Counter();rows=0;matches=0
    fields=['domain_config_id','pfam_accession','boundary','mask','target_domain_pair','background_domain_pair','target_order','background_order','contrast_status','target_rmsd_status','background_rmsd_status','target_geometry_status','background_geometry_status']
    fields += [side+'_'+metric for side in ['target','background'] for metric in METRICS]
    fields += ['rmsd_target_minus_background','sequence_identity_target_minus_background','aligned_length_target_minus_background']
    for screen in plan['screens']:fields += [screen+'_mask_pass',screen+'_both_masks_pass']
    with gzip.open(out/'matched_domain_contrasts.tsv.gz','wt') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for line in (inventory/'domain_configurations.jsonl').open():
            cfg=json.loads(line)
            for m in cfg['matches']:
                matches+=1
                for mask in ['full','plddt70']:
                    for to in [0,1]:
                        for bo in [0,1]:
                            row=dict(domain_config_id=cfg['domain_config_id'],pfam_accession=m['pfam_accession'],boundary=m['boundary'],mask=mask,target_domain_pair=m['target_domain_pair'],background_domain_pair=m['background_domain_pair'],target_order=to,background_order=bo);usable=[];values={}
                            for label,order in [('target',to),('background',bo)]:
                                ns,cs,qs=sources[label];pair=m[label+'_domain_pair'];c=cs[pair,mask];q=qs[pair,mask];n=ns.get((pair,mask,order));assert (n is None)==(c[f'order{order}_status']=='input_unavailable')
                                row[label+'_rmsd_status']=n['rmsd_status'] if n else 'input_unavailable';row[label+'_geometry_status']=q[f'order{order}_geometry_status']
                                ok=n is not None and n['rmsd_status']=='within_printed_rounding' and int(n['aligned_length'])>=3 and q[f'order{order}_geometry_status']=='unique_at_numeric_tolerance';usable.append(ok)
                                if ok:
                                    values[label]=measurement(n,c,order)
                                    for metric,v in values[label].items():row[label+'_'+metric]=v
                            ready=all(usable) and cfg['status']=='shared_domain_comparisons';row['contrast_status']='numerically_computable' if ready else 'excluded_numerically_or_identical'
                            if ready:
                                for metric,name in [('rmsd_recomputed','rmsd'),('sequence_identity_exact','sequence_identity'),('aligned_length','aligned_length')]:row[name+'_target_minus_background']=values['target'][metric]-values['background'][metric]
                            for screen in plan['screens']:
                                local=all(sources[label][2][m[label+'_domain_pair'],mask][screen+'_pass']=='1' for label in ['target','background']) and cfg['status']=='shared_domain_comparisons'
                                common=all(sources[label][2][m[label+'_domain_pair'],mm][screen+'_pass']=='1' for label in ['target','background'] for mm in ['full','plddt70']) and cfg['status']=='shared_domain_comparisons'
                                assert not local or ready;row[screen+'_mask_pass']=int(local);row[screen+'_both_masks_pass']=int(common)
                            w.writerow(row);rows+=1;counts[mask+':'+row['contrast_status']]+=1
    assert rows==matches*8;verify()
    r=dict(status='complete_matched_domain_measurements_pending_independent_readback',plan_sha256=ph,shared_boundary_matches=matches,rows=rows,counts=dict(counts),inventory_receipt_sha256=sha(inventory/'receipt.json'),qualification_receipt_sha256=sha(qualified/'receipt.json'),artifacts={'matched_domain_contrasts.tsv.gz':sha(out/'matched_domain_contrasts.tsv.gz')},scope='All shared-domain configurations, two masks and four target/background input-order combinations. Unusable side metrics blank, numerical/geometry status retained, six same-mask and both-mask flags explicit. RMSD difference is a descriptive difference between separately aligned cores, not a matched-residue four-protein fit, error bound or biological effect. Endpoint TM scores and original/retained coverage remain separate; no preferred order or averaging. Event projection and phylogenetic/family-aware inference remain downstream.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2),flush=True)

if __name__=='__main__':main()
