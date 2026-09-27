#!/usr/bin/env python3
"""Measure coordinate rank and rigid-fit curvature for every recorded whole-protein reference alignment."""
import argparse,csv,json,time
import psutil
from collections import Counter
from functools import lru_cache
from pathlib import Path
import numpy as np
from duplication_alignment_numeric_readback import load_pdb
from duplication_alignment_numeric_diagnostic import check_alignment
from screen_duplication_domain_alignment_coverage import sha


def geometry(x,y):
    if x.shape!=y.shape or x.ndim!=2 or x.shape[1]!=3 or len(x)<1:raise ValueError('Invalid coordinate arrays')
    if not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('Nonfinite coordinates')
    a=x-x.mean(axis=0);b=y-y.mean(axis=0);n=len(a)
    sx=np.pad(np.linalg.svd(a,compute_uv=False),(0,max(0,3-min(a.shape))))
    sy=np.pad(np.linalg.svd(b,compute_uv=False),(0,max(0,3-min(b.shape))))
    u,s,vt=np.linalg.svd(a.T@b);sign=1 if np.linalg.det(u@vt)>=0 else -1
    tol=np.finfo(float).eps*max(n,3)
    rx=int(sum(sx>tol*sx[0]));ry=int(sum(sy>tol*sy[0]))
    # Rotation-Hessian eigenvalues are proportional to pairwise sums of the
    # determinant-corrected singular values; the smallest is s2 + sign*s3.
    curvature=float(s[1]+sign*s[2]);scale=float(s[0])
    unique=curvature>tol*scale
    status='unique_at_numeric_tolerance' if unique else 'degenerate_at_numeric_tolerance'
    return dict(aligned_length=n,rank_left=rx,rank_right=ry,
                rms_width1_left=float(sx[0]/np.sqrt(n)),rms_width2_left=float(sx[1]/np.sqrt(n)),rms_width3_left=float(sx[2]/np.sqrt(n)),
                rms_width1_right=float(sy[0]/np.sqrt(n)),rms_width2_right=float(sy[1]/np.sqrt(n)),rms_width3_right=float(sy[2]/np.sqrt(n)),
                width2_to_width1_left=float(sx[1]/sx[0]) if sx[0]>0 else 0.,width2_to_width1_right=float(sy[1]/sy[0]) if sy[0]>0 else 0.,
                cross_s1=float(s[0]),cross_s2=float(s[1]),cross_s3=float(s[2]),determinant_correction=sign,
                minimum_rotation_curvature=curvature,relative_rotation_curvature=curvature/scale if scale>0 else 0.,
                relative_numeric_tolerance=tol,geometry_status=status)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed pinned input: '+p)
    verify();identity=plan['producer']
    while True:
        try:
            process=psutil.Process(identity['pid'])
            if process.create_time()!=identity['created'] or process.status()==psutil.STATUS_ZOMBIE:break
            assert process.cmdline()==identity['cmdline']
        except psutil.NoSuchProcess:break
        print('waiting_for_complete_reference_diagnostic',identity['pid'],flush=True);time.sleep(30)
    verify();diagnostic=Path(plan['diagnostic']);dr=json.loads((diagnostic/'receipt.json').read_text())
    if dr['status']!='complete_reference_alignment_rmsd_diagnostic_not_scientific_acceptance':raise ValueError('Wrong diagnostic status')
    if sha(diagnostic/'numeric_readback.tsv')!=dr['artifacts']['numeric_readback.tsv']:raise ValueError('Changed diagnostic table')
    dp=json.loads(Path(plan['diagnostic_plan']).read_text());source=json.loads(Path(dp['source_plan']).read_text());root=Path(source['output'])
    if dr['plan_sha256']!=sha(plan['diagnostic_plan']) or dp['output']!=str(diagnostic):raise ValueError('Diagnostic provenance mismatch')
    if sha(root/'receipt.json')!=dr['producer_receipt_sha256']:raise ValueError('Changed alignment producer')
    producer=json.loads((root/'receipt.json').read_text());manifest=root/'checkpoint_manifest.tsv'
    if sha(manifest)!=producer['artifacts']['checkpoint_manifest.tsv']:raise ValueError('Changed checkpoint manifest')
    proofs={r['path']:r['sha256'] for r in csv.DictReader(manifest.open(),delimiter='\t')}
    inputs={}
    for spec in source['input_sources'].values():
        folder=Path(spec['inputs']);ir=json.loads((folder/'receipt.json').read_text())
        if sha(folder/'inputs.jsonl')!=ir['artifacts']['inputs.jsonl']:raise ValueError('Changed PDB manifest')
        for line in (folder/'inputs.jsonl').open():
            row=json.loads(line);key=row['model_id'],row['version'],row['mask']
            if key in inputs:raise ValueError('Duplicate input identity')
            inputs[key]=row
    @lru_cache(maxsize=256)
    def coords(key):return load_pdb(inputs[key])
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);counts=Counter();seen=set();maxerr=0.
    fields=['pair_key','mask','order','rmsd_status']+list(geometry(np.zeros((3,3)),np.zeros((3,3))))
    with (out/'alignment_geometry.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in csv.DictReader((diagnostic/'numeric_readback.tsv').open(),delimiter='\t'):
            pair,mask,order=row['pair_key'],row['mask'],int(row['order']);key=pair,mask,order
            if key in seen:raise ValueError('Repeated alignment');
            seen.add(key);rel=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json';path=root/rel
            if sha(path)!=proofs[rel]:raise ValueError('Changed checkpoint')
            record=json.loads(path.read_text())
            if record['status']!='aligned' or (record['pair_key'],record['mask'],record['order'])!=key:raise ValueError('Wrong checkpoint')
            left,right=[coords((r['model_id'],r['version'],mask)) for r in record['inputs']]
            numeric=check_alignment(record,left,right)
            for name,value in numeric.items():
                if name=='rmsd_status':
                    if row[name]!=value:raise ValueError('Changed diagnostic classification')
                elif float(row[name])!=float(value):raise ValueError('Changed diagnostic metric')
            x,y=record['metrics']['alignment_left'],record['metrics']['alignment_right'];i=j=0;ix=[];iy=[]
            for a,b in zip(x,y):
                if a!='-' and b!='-':ix.append(i);iy.append(j)
                i+=a!='-';j+=b!='-'
            g=geometry(left[1][ix],right[1][iy])
            if g['aligned_length']!=int(row['aligned_length']):raise ValueError('Mapping length mismatch')
            w.writerow(dict(pair_key=pair,mask=mask,order=order,rmsd_status=row['rmsd_status'],**g))
            counts[mask+':'+g['geometry_status']]+=1;maxerr=max(maxerr,numeric['rmsd_rounding_error'])
            if len(seen)%10000==0:print('Geometry checked',len(seen),'/',dr['numerically_checked_alignments'],flush=True)
    if len(seen)!=dr['numerically_checked_alignments']:raise ValueError('Incomplete geometry cohort')
    verify()
    result=dict(status='complete_reference_alignment_geometry_pending_independent_readback',plan_sha256=ph,diagnostic_receipt_sha256=sha(diagnostic/'receipt.json'),alignments=len(seen),counts=dict(counts),maximum_rmsd_discrepancy=maxerr,artifacts={'alignment_geometry.tsv':sha(out/'alignment_geometry.tsv')},scope='All successful directed whole-protein reference alignments; native mappings and numeric diagnostic rechecked from hashed PDBs, coordinate singular spectra and optimal proper-rotation curvature recorded. Numerical identifiability, not prediction accuracy or stability to coordinate/model uncertainty. Existing RMSD quarantine unchanged; no biological inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
