#!/usr/bin/env python3
"""Audit common-core RMSDs through quaternion rotations and complete source reconstruction."""
import argparse,csv,gzip,itertools,json,math
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.linalg import svd
from duplication_alignment_numeric_readback import load_pdb
from screen_duplication_domain_alignment_coverage import sha


def quaternion_fit(x,y):
    a=x-x.mean(0);b=y-y.mean(0);h=a.T@b;t=np.trace(h)
    z=np.array([h[1,2]-h[2,1],h[2,0]-h[0,2],h[0,1]-h[1,0]])
    k=np.empty((4,4));k[0,0]=t;k[0,1:]=z;k[1:,0]=z;k[1:,1:]=h+h.T-t*np.eye(3)
    values,vectors=np.linalg.eigh(k);q=vectors[:,-1];w=q[0];v=q[1:]
    cross=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    rotation=(w*w-v@v)*np.eye(3)+2*np.outer(v,v)+2*w*cross
    if not np.allclose(rotation.T@rotation,np.eye(3),atol=1e-12) or abs(np.linalg.det(rotation)-1)>1e-12:raise ValueError('Invalid proper quaternion rotation')
    distance=float(np.linalg.norm(a@rotation.T-b)/np.sqrt(len(a)))
    scale=float(svd(h,compute_uv=False,lapack_driver='gesvd')[0]);curvature=float((values[-1]-values[-2])/2)
    relative=curvature/scale if scale>0 else 0.
    unique=relative>np.finfo(float).eps*max(len(a),3)
    return distance,relative,'unique_at_numeric_tolerance' if unique else 'degenerate_at_numeric_tolerance'


def core_metrics(coords,letters,confidence):
    n=len(coords[0]);result={};statuses=[]
    for label,i,j in [('ab',0,1),('ar',0,2),('br',1,2)]:
        distance,curvature,status=quaternion_fit(coords[i],coords[j])
        result['rmsd_'+label]=distance;result['relative_curvature_'+label]=curvature;result['geometry_'+label]=status
        result['sequence_identity_'+label]=sum(a==b for a,b in zip(letters[i],letters[j]))/n;statuses.append(status)
    result['rmsd_ar_minus_br']=result['rmsd_ar']-result['rmsd_br']
    for name,p in zip(['a','b','reference'],confidence):result['mean_plddt_'+name]=sum(map(float,p))/n
    result['joint_plddt70_fraction']=sum(all(p[i]>=70 for p in confidence) for i in range(n))/n
    result['fit_status']='computed_unique_at_numeric_tolerance' if all(s=='unique_at_numeric_tolerance' for s in statuses) else 'computed_degenerate_geometry'
    if abs(result['rmsd_ar_minus_br'])>result['rmsd_ab']+1e-9:raise ValueError('Common-core RMSD metric bound violated')
    return result


def compare_row(actual,wanted):
    maximum=0.
    if set(actual)!=set(wanted):raise ValueError('Unexpected output columns')
    for field,value in wanted.items():
        if isinstance(value,float):
            observed=float(actual[field])
            if not math.isfinite(observed) or not math.isclose(observed,value,rel_tol=1e-9,abs_tol=1e-9):raise ValueError('Numeric mismatch: '+field)
            if field.startswith('rmsd_'):maximum=max(maximum,abs(observed-value))
        elif actual[field]!=str(value):raise ValueError('Field mismatch: '+field)
    return maximum


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed source: '+path)
    verify();out=Path(plan['output']);rp=out/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    if r['status']!='complete_common_residue_domain_fits_pending_readback' or r['plan_sha256']!=ph:raise ValueError('Wrong completed source')
    for name,h in r['artifacts'].items():
        if sha(out/name)!=h:raise ValueError('Changed fitted table')
    mp=json.loads(Path(plan['mapping_plan']).read_text());root=Path(mp['output']);mr=json.loads((root/'receipt.json').read_text())
    if sha(root/'receipt.json')!=r['mapping_receipt_sha256'] or sha(plan['mapping_readback'])!=r['mapping_audit_sha256']:raise ValueError('Unbound mapping sources')
    for name,h in mr['artifacts'].items():
        if sha(root/name)!=h:raise ValueError('Changed common mapping artifact')
    gp=json.loads(Path(mp['geometry_plan']).read_text());dp=json.loads(Path(gp['diagnostic_plan']).read_text());sp=json.loads(Path(dp['source_plan']).read_text());folder=Path(sp['inputs']);ir=json.loads((folder/'receipt.json').read_text())
    if sha(folder/'inputs.jsonl')!=ir['artifacts']['inputs.jsonl']:raise ValueError('Changed input manifest')
    inputs={}
    for line in (folder/'inputs.jsonl').open():
        x=json.loads(line);inputs[x['interval_id'],x['mask']]=x
    triads={}
    for line in (root/'triads.jsonl').open():
        x=json.loads(line);triads[x['triad_key']]=[x[f] for f in ['interval_a','interval_b','interval_reference']]
    @lru_cache(maxsize=256)
    def coordinates(key):
        row=inputs[key];seq,xyz,p=load_pdb(row)
        return {position:(xyz[i],seq[i],p[i]) for i,position in enumerate(row['original_positions'])}
    @lru_cache(maxsize=64)
    def metrics(tk,mask,triples):
        c=[coordinates((i,mask)) for i in triads[tk]];parts=[[c[j][row[j]] for row in triples] for j in range(3)]
        return core_metrics([np.array([v[0] for v in part]) for part in parts],[''.join(v[1] for v in part) for part in parts],[[v[2] for v in part] for part in parts])
    metric_fields=list(core_metrics([np.eye(3)]*3,['AAA']*3,[[90.]*3]*3))
    counts=Counter();passes=Counter();seen=set();rows=0;maximum=0.
    with (out/'common_residue_fits.tsv').open() as f,gzip.open(root/'common_residue_maps.jsonl.gz','rt') as src:
        reader=csv.DictReader(f,delimiter='\t')
        for line in src:
            item=json.loads(line);tk=item['triad_key'];mask=item['mask'];orders=tuple(item['orders']);key=tk,mask,orders
            if key in seen:raise ValueError('Repeated source mapping');
            seen.add(key);lengths=[inputs[i,mask]['interval_length'] for i in triads[tk]]
            source_problems=sorted(set(itertools.chain.from_iterable(item['edge_exclusions'])))
            for definition,field in [('reference_common','reference_common_triples'),('cycle_consistent','cycle_consistent_triples')]:
                triples=tuple(tuple(v) for v in item[field]);n=len(triples)
                wanted=dict(triad_key=tk,mask=mask,order_ab=orders[0],order_ar=orders[1],order_br=orders[2],mapping_definition=definition,common_residues=n,length_a=lengths[0],length_b=lengths[1],length_reference=lengths[2],coverage_a=n/lengths[0],coverage_b=n/lengths[1],coverage_reference=n/lengths[2],source_exclusions=';'.join(source_problems),mapping_disagreement_count=item['common_reference_count']-item['cycle_consistent_count'])
                wanted.update({k:'' for k in metric_fields})
                if source_problems:wanted['fit_status']='source_excluded'
                elif n<3:wanted['fit_status']='fewer_than_three_common_residues'
                else:wanted.update(metrics(tk,mask,triples))
                for screen in plan['screens']:
                    why=list(source_problems)
                    if n<screen['minimum_aligned_residues']:why.append('short_common_core')
                    cutoff=Fraction(str(screen['minimum_original_coverage']))
                    if any(n*cutoff.denominator<l*cutoff.numerator for l in lengths):why.append('low_original_interval_coverage')
                    if wanted['fit_status']=='computed_degenerate_geometry':why.append('degenerate_common_core_geometry')
                    if wanted['fit_status']=='fewer_than_three_common_residues':why.append('fewer_than_three_common_residues')
                    okay=not why;wanted[screen['id']+'_pass']=int(okay);wanted[screen['id']+'_exclusions']=';'.join(why)
                    passes[mask+':'+definition+':'+screen['id']]+=okay
                actual=next(reader,None)
                if actual is None:raise ValueError('Missing fit row')
                maximum=max(maximum,compare_row(actual,wanted));rows+=1;counts[mask+':'+definition+':'+wanted['fit_status']]+=1
            if len(seen)%10000==0:print('Independent common-core fit checks',len(seen),'/',r['mapping_dispositions'],flush=True)
        if next(reader,None) is not None:raise ValueError('Extra fit row')
    expected={(tk,m,o) for tk in triads for m in ['full','plddt70'] for o in itertools.product([0,1],repeat=3)}
    if seen!=expected or len(seen)!=r['mapping_dispositions'] or rows!=r['fit_rows'] or dict(counts)!=r['counts'] or dict(passes)!=r['screen_pass_counts']:raise ValueError('Incomplete fit readback')
    verify()
    if sha(rp)!=rh:raise ValueError('Producer receipt changed')
    for name,h in r['artifacts'].items():
        if sha(out/name)!=h:raise ValueError('Fit table changed during audit')
    result=dict(status='passed_full_common_residue_domain_fit_readback',plan_sha256=ph,producer_receipt_sha256=rh,checker_sha256=sha(__file__),fit_rows_checked=rows,mapping_dispositions=len(seen),maximum_absolute_rmsd_or_contrast_difference=maximum,counts=dict(counts),screen_pass_counts=dict(passes),scope='All records reconstructed from source residue maps and hashed PDBs; quaternion rotations independently reproduce SVD-fit RMSDs and signed differences. Every confidence/identity/geometry value, blank/exclusion and rational coverage decision checked. Original input exclusions preserved. No biological significance or asymmetric evolution inference.')
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
