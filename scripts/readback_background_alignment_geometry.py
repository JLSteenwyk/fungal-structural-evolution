#!/usr/bin/env python3
"""Read back complete background geometry using an independent quaternion eigensystem."""
import argparse,csv,json,time,subprocess
from collections import Counter
from functools import lru_cache
from pathlib import Path
import numpy as np
import psutil
from scipy.linalg import svd
from duplication_alignment_numeric_readback import load_pdb
from screen_duplication_domain_alignment_coverage import sha


def quaternion_curvature(x,y):
    a=x-np.mean(x,axis=0);b=y-np.mean(y,axis=0);h=a.T@b;t=np.trace(h)
    z=np.array([h[1,2]-h[2,1],h[2,0]-h[0,2],h[0,1]-h[1,0]])
    k=np.empty((4,4));k[0,0]=t;k[0,1:]=z;k[1:,0]=z;k[1:,1:]=h+h.T-t*np.eye(3)
    e=np.linalg.eigvalsh(k)
    return float((e[-1]-e[-2])/2),float(e[-1])


def verify_row(row,x,y):
    n=len(x);tol=np.finfo(float).eps*max(n,3);a=x-x.mean(0);b=y-y.mean(0)
    if int(row['aligned_length'])!=n or n!=len(y):raise ValueError('Aligned count differs')
    spectra=[]
    for side,coords in [('left',a),('right',b)]:
        spectrum=np.pad(svd(coords,compute_uv=False,lapack_driver='gesvd'),(0,max(0,3-min(coords.shape))))
        spectra.append(spectrum)
        if int(row['rank_'+side])!=int(np.sum(spectrum>tol*spectrum[0])):raise ValueError('Coordinate rank differs')
        for i in range(3):
            observed=float(row[f'rms_width{i+1}_{side}']);expected=float(spectrum[i]/np.sqrt(n))
            if not np.isclose(observed,expected,rtol=1e-10,atol=1e-10):raise ValueError('Coordinate width differs')
        ratio=spectrum[1]/spectrum[0] if spectrum[0]>0 else 0.
        if not np.isclose(float(row['width2_to_width1_'+side]),ratio,rtol=1e-10,atol=1e-10):raise ValueError('Coordinate width ratio differs')
    cross=svd(a.T@b,compute_uv=False,lapack_driver='gesvd')
    if not np.allclose([float(row[f'cross_s{i+1}']) for i in range(3)],cross,rtol=1e-10,atol=1e-10):raise ValueError('Cross spectrum differs')
    qcurv,qopt=quaternion_curvature(x,y);reported=float(row['minimum_rotation_curvature']);scale=max(1.,float(cross[0]))
    err=abs(qcurv-reported)/scale
    if err>1e-10:raise ValueError('Quaternion curvature disagrees')
    delta=int(row['determinant_correction'])
    if delta not in [-1,1] or abs(qopt-(cross[0]+cross[1]+delta*cross[2]))/scale>1e-10:raise ValueError('Proper-rotation optimum differs')
    native_cross=[float(row[f'cross_s{i+1}']) for i in range(3)]
    if reported!=native_cross[1]+delta*native_cross[2]:raise ValueError('Curvature algebra differs')
    if float(row['relative_numeric_tolerance'])!=tol:raise ValueError('Numeric tolerance differs')
    relative=reported/native_cross[0] if native_cross[0]>0 else 0.
    if float(row['relative_rotation_curvature'])!=relative:raise ValueError('Relative curvature differs')
    status='unique_at_numeric_tolerance' if reported>tol*native_cross[0] else 'degenerate_at_numeric_tolerance'
    if row['geometry_status']!=status:raise ValueError('Geometry status differs')
    # Close-to-zero quaternion gaps are recorded; cross-method agreement at
    # 1e-10 scaled error is not a claim to resolve ranks at machine precision.
    return err, bool(qcurv/scale<=64*tol)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed readback plan')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed source: '+p)
    verify();dep=plan['producer']
    while True:
        try:
            process=psutil.Process(dep['pid'])
            if abs(process.create_time()-dep['created'])>.01 or process.status()==psutil.STATUS_ZOMBIE:break
            if process.cmdline()!=dep['cmdline']:raise ValueError('Changed producer identity')
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    if state!=dict(ActiveState='inactive',Result='success',ExecMainStatus='0'):raise ValueError('Geometry producer did not finish successfully: '+str(state))
    verify();gp=json.loads(Path(plan['source_plan']).read_text());root=Path(gp['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    if r['status']!='complete_background_alignment_geometry_pending_independent_readback' or r['plan_sha256']!=sha(plan['source_plan']):raise ValueError('Incomplete geometry source')
    table=root/'alignment_geometry.tsv'
    if sha(table)!=r['artifacts']['alignment_geometry.tsv']:raise ValueError('Changed geometry table')
    diagnostic=Path(gp['diagnostic']);dr=json.loads((diagnostic/'receipt.json').read_text())
    if sha(diagnostic/'receipt.json')!=r['diagnostic_receipt_sha256'] or sha(diagnostic/'numeric_readback.tsv')!=dr['artifacts']['numeric_readback.tsv']:raise ValueError('Changed diagnostic source')
    expected={}
    for row in csv.DictReader((diagnostic/'numeric_readback.tsv').open(),delimiter='\t'):
        key=row['pair_key'],row['mask'],int(row['order'])
        if key in expected:raise ValueError('Duplicate diagnostic identity')
        expected[key]=row
    if len(expected)!=dr['numerically_checked_alignments']:raise ValueError('Incomplete diagnostic table')
    if dr['status']!='complete_background_alignment_rmsd_diagnostic_not_scientific_acceptance' or dr['plan_sha256']!=sha(gp['diagnostic_plan']):raise ValueError('Diagnostic provenance mismatch')
    dp=json.loads(Path(gp['diagnostic_plan']).read_text());source=json.loads(Path(dp['source_plan']).read_text());alignments=Path(source['output'])
    if sha(alignments/'receipt.json')!=dr['producer_receipt_sha256']:raise ValueError('Changed alignment source')
    ar=json.loads((alignments/'receipt.json').read_text());manifest=alignments/'checkpoint_manifest.tsv'
    if sha(manifest)!=ar['artifacts']['checkpoint_manifest.tsv']:raise ValueError('Changed checkpoint manifest')
    proofs={row['path']:row['sha256'] for row in csv.DictReader(manifest.open(),delimiter='\t')}
    inputs={}
    for spec in source['input_sources'].values():
        folder=Path(spec['inputs']);ir=json.loads((folder/'receipt.json').read_text())
        if sha(folder/'inputs.jsonl')!=ir['artifacts']['inputs.jsonl']:raise ValueError('Changed PDB manifest')
        for line in (folder/'inputs.jsonl').open():
            row=json.loads(line);key=row['model_id'],row['version'],row['mask']
            if key in inputs:raise ValueError('Duplicate input identity')
            inputs[key]=row
    @lru_cache(maxsize=256)
    def coordinates(key):return load_pdb(inputs[key])
    seen=set();counts=Counter();near=0;maximum=0.;degenerate=[]
    for row in csv.DictReader(table.open(),delimiter='\t'):
        key=row['pair_key'],row['mask'],int(row['order'])
        if key in seen or key not in expected:raise ValueError('Unexpected/repeated record')
        seen.add(key);n=expected[key]
        if row['rmsd_status']!=n['rmsd_status'] or row['aligned_length']!=n['aligned_length']:raise ValueError('Diagnostic mapping differs')
        pair,mask,order=key;rel=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json';path=alignments/rel
        if sha(path)!=proofs[rel]:raise ValueError('Changed native checkpoint')
        record=json.loads(path.read_text());coords=[coordinates((v['model_id'],v['version'],mask)) for v in record['inputs']]
        strings=[record['metrics']['alignment_left'],record['metrics']['alignment_right']]
        nongap=[np.array(list(s))!='-' for s in strings];paired=nongap[0]&nongap[1]
        for text,c in zip(strings,coords):
            if text.replace('-','')!=c[0]:raise ValueError('Residue sequence mismatch')
        indexes=[(np.cumsum(p)-1)[paired] for p in nongap]
        err,near_boundary=verify_row(row,*[c[1][i] for c,i in zip(coords,indexes)])
        maximum=max(maximum,err);near+=near_boundary;counts[mask+':'+row['geometry_status']]+=1
        if row['geometry_status']=='degenerate_at_numeric_tolerance':degenerate.append(dict(pair_key=pair,mask=mask,order=order,aligned_length=int(row['aligned_length']),rmsd_status=row['rmsd_status']))
        if len(seen)%10000==0:print('Independent geometry checks',len(seen),'/',len(expected),flush=True)
    if seen!=set(expected) or len(seen)!=r['alignments'] or dict(counts)!=r['counts']:raise ValueError('Incomplete geometry readback')
    verify()
    if sha(rp)!=rh or sha(table)!=r['artifacts']['alignment_geometry.tsv']:raise ValueError('Geometry source changed during audit')
    result=dict(status='passed_full_background_geometry_readback',plan_sha256=ph,producer_receipt_sha256=rh,alignments_checked=len(seen),counts=dict(counts),maximum_scaled_quaternion_curvature_error=maximum,near_zero_quaternion_gaps=near,degenerate_alignments=degenerate,scope='Every geometry row checked from hashed native mapping/PDBs. Coordinate and cross spectra checked using LAPACK gesvd; curvature/optimum checked using independent 4x4 quaternion eigensystem. Status verified under original numeric tolerance; near-zero gaps recorded, not treated as cross-method proof of machine-precision rank. No prediction-uncertainty or biological inference claim.')
    out=Path(plan['output']);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='degenerate_alignments'},indent=2))


if __name__=='__main__':main()
