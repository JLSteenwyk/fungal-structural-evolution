#!/usr/bin/env python3
"""Recompute all three proper-rigid distances on identical common residue triples."""
import argparse,csv,gzip,itertools,json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import numpy as np
from assess_domain_alignment_geometry import geometry
from duplication_alignment_numeric_readback import load_pdb
from readback_cross_clan_alignments import rmsd
from screen_duplication_domain_alignment_coverage import sha


def fit_triplet(coordinates,letters,confidence):
    n=len(coordinates[0])
    if n<3 or any(len(c)!=n for c in coordinates):raise ValueError('At least three common residues required')
    result={};statuses=[]
    for name,i,j in [('ab',0,1),('ar',0,2),('br',1,2)]:
        g=geometry(coordinates[i],coordinates[j]);result['rmsd_'+name]=rmsd(coordinates[i],coordinates[j])
        result['geometry_'+name]=g['geometry_status'];result['relative_curvature_'+name]=g['relative_rotation_curvature']
        result['sequence_identity_'+name]=sum(a==b for a,b in zip(letters[i],letters[j]))/n
        statuses.append(g['geometry_status'])
    result['rmsd_ar_minus_br']=result['rmsd_ar']-result['rmsd_br']
    for label,p in zip(['a','b','reference'],confidence):result['mean_plddt_'+label]=float(np.mean(p))
    result['joint_plddt70_fraction']=float(np.mean(np.logical_and.reduce([p>=70 for p in confidence])))
    result['fit_status']='computed_unique_at_numeric_tolerance' if all(s=='unique_at_numeric_tolerance' for s in statuses) else 'computed_degenerate_geometry'
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed source: '+p)
    verify();mp=json.loads(Path(plan['mapping_plan']).read_text());root=Path(mp['output']);mr=json.loads((root/'receipt.json').read_text());audit=json.loads(Path(plan['mapping_readback']).read_text())
    if audit['status']!='passed_full_domain_common_residue_readback' or audit['producer_receipt_sha256']!=sha(root/'receipt.json') or audit['plan_sha256']!=sha(plan['mapping_plan']):raise ValueError('Unbound common map audit')
    for name,h in mr['artifacts'].items():
        if sha(root/name)!=h:raise ValueError('Changed common maps')
    gp=json.loads(Path(mp['geometry_plan']).read_text());dp=json.loads(Path(gp['diagnostic_plan']).read_text());sp=json.loads(Path(dp['source_plan']).read_text());folder=Path(sp['inputs']);ir=json.loads((folder/'receipt.json').read_text())
    if sha(folder/'inputs.jsonl')!=ir['artifacts']['inputs.jsonl']:raise ValueError('Changed input manifest')
    inputs={}
    for line in (folder/'inputs.jsonl').open():
        r=json.loads(line);inputs[r['interval_id'],r['mask']]=r
    triads={}
    for line in (root/'triads.jsonl').open():
        r=json.loads(line);triads[r['triad_key']]=[r[f] for f in ['interval_a','interval_b','interval_reference']]
    @lru_cache(maxsize=256)
    def pdb(key):
        row=inputs[key];letters,xyz,p=load_pdb(row)
        return letters,xyz,p,{position:i for i,position in enumerate(row['original_positions'])}
    @lru_cache(maxsize=64)
    def calculate(tk,mask,triples):
        sources=[pdb((i,mask)) for i in triads[tk]];coords=[];letters=[];conf=[]
        for column,(seq,xyz,p,index) in enumerate(sources):
            ix=[index[t[column]] for t in triples]
            coords.append(xyz[ix]);letters.append(''.join(seq[i] for i in ix));conf.append(p[ix])
        return fit_triplet(coords,letters,conf)
    sample=fit_triplet([np.eye(3)]*3,['AAA']*3,[np.array([90.]*3)]*3)
    metrics=list(sample);fields=['triad_key','mask','order_ab','order_ar','order_br','mapping_definition','common_residues','length_a','length_b','length_reference','coverage_a','coverage_b','coverage_reference','source_exclusions','mapping_disagreement_count']+metrics
    for screen in plan['screens']:fields += [screen['id']+'_pass',screen['id']+'_exclusions']
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);seen=set();counts=Counter();passed=Counter();rows=0
    with gzip.open(root/'common_residue_maps.jsonl.gz','rt') as src,(out/'common_residue_fits.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for line in src:
            item=json.loads(line);tk=item['triad_key'];mask=item['mask'];orders=tuple(item['orders']);key=tk,mask,orders
            if key in seen:raise ValueError('Repeated mapping record')
            seen.add(key);lengths=item['original_interval_lengths'];source_problems=sorted(set(itertools.chain.from_iterable(item['edge_exclusions'])))
            if lengths!=[inputs[i,mask]['interval_length'] for i in triads[tk]]:raise ValueError('Original length mismatch')
            for definition,field in [('reference_common','reference_common_triples'),('cycle_consistent','cycle_consistent_triples')]:
                triples=tuple(tuple(x) for x in item[field]);n=len(triples)
                if any(len({t[c] for t in triples})!=n for c in [0,1,2]):raise ValueError('Repeated residue in common map')
                row=dict(triad_key=tk,mask=mask,order_ab=orders[0],order_ar=orders[1],order_br=orders[2],mapping_definition=definition,common_residues=n,
                    length_a=lengths[0],length_b=lengths[1],length_reference=lengths[2],coverage_a=n/lengths[0],coverage_b=n/lengths[1],coverage_reference=n/lengths[2],source_exclusions=';'.join(source_problems),mapping_disagreement_count=item['common_reference_count']-item['cycle_consistent_count'])
                if source_problems:row['fit_status']='source_excluded'
                elif n<3:row['fit_status']='fewer_than_three_common_residues'
                else:row.update(calculate(tk,mask,triples))
                for screen in plan['screens']:
                    why=list(source_problems)
                    if n<screen['minimum_aligned_residues']:why.append('short_common_core')
                    cutoff=Fraction(str(screen['minimum_original_coverage']))
                    if any(Fraction(n,l)<cutoff for l in lengths):why.append('low_original_interval_coverage')
                    if row['fit_status']=='computed_degenerate_geometry':why.append('degenerate_common_core_geometry')
                    if row['fit_status']=='fewer_than_three_common_residues':why.append('fewer_than_three_common_residues')
                    okay=not why
                    if okay and row['fit_status']!='computed_unique_at_numeric_tolerance':raise ValueError('Uncomputed fit passed screen')
                    row[screen['id']+'_pass']=int(okay);row[screen['id']+'_exclusions']=';'.join(why)
                    passed[mask+':'+definition+':'+screen['id']]+=okay
                writer.writerow(row);rows+=1;counts[mask+':'+definition+':'+row['fit_status']]+=1
            if len(seen)%10000==0:print('Common-core fit dispositions',len(seen),'/',mr['mask_order_dispositions'],flush=True)
    expected={(tk,mask,orders) for tk in triads for mask in ['full','plddt70'] for orders in itertools.product([0,1],repeat=3)}
    if seen!=expected or rows!=mr['mask_order_dispositions']*2:raise ValueError('Incomplete fit universe')
    verify()
    result=dict(status='complete_common_residue_domain_fits_pending_readback',plan_sha256=ph,mapping_receipt_sha256=sha(root/'receipt.json'),mapping_audit_sha256=sha(plan['mapping_readback']),mapping_dispositions=len(seen),fit_rows=rows,counts=dict(counts),screen_pass_counts=dict(passed),artifacts={'common_residue_fits.tsv':sha(out/'common_residue_fits.tsv')},scope='All common-reference and cycle-consistent mapping alternatives; RMSDs and sequence identities use identical residue triples for A-B, A-reference and B-reference. Signed RMSD difference retains A/B orientation. Common-core geometry and original-interval coverage re-evaluated; upstream exclusions preserved. No order selection, significance test, prediction-error calibration or claim of asymmetric evolution.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
