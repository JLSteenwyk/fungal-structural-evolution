#!/usr/bin/env python3
"""Read back four platform-check outputs and compare matched CA coordinates."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coords(path):
    lines = [s for s in path.read_text().splitlines() if s.startswith('ATOM') and s[12:16].strip() == 'CA']
    return [(s[17:20], s[21:27]) for s in lines], np.array([[float(s[a:b]) for a,b in [(30,38),(38,46),(46,54)]] for s in lines])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--predictions',type=Path,required=True)
    ap.add_argument('--inputs',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    inputs=json.loads(a.inputs.read_text());chunk=json.loads((a.predictions/'last_chunk.json').read_text())
    assert chunk['new_predictions']==4 and chunk['remaining_eligible']==0 and not chunk['interrupted'] and chunk['oom_deferred']==0
    rows=[]
    for item in inputs['sequences']:
        sid=item['sequence_id'];local=Path(item['local_receipt']);assert sha(local)==item['local_receipt_sha256']
        remote=a.predictions/(sid+'.json');pairs=[];receipts=[]
        for p in [local,remote]:
            r=json.loads(p.read_text());assert r['status']=='verified_prediction' and r['sequence_id']==sid
            for name,digest in r['artifacts'].items():assert sha(p.parent/name)==digest
            assert sha(p.parent/'config.json')==r['config_sha256']
            pairs.append(coords(p.parent/(sid+'.pdb')));receipts.append(r)
        (ids,x),(other,y)=pairs;assert ids==other and len(x)==item['length'];assert np.isfinite(x).all() and np.isfinite(y).all()
        x=x-x.mean(axis=0);y=y-y.mean(axis=0)
        u,_,vt=np.linalg.svd(x.T@y);d=np.eye(3);d[-1,-1]=np.sign(np.linalg.det(u@vt));rotation=u@d@vt
        rmsd=float(np.sqrt(np.mean(np.sum((x@rotation-y)**2,axis=1))))
        arrays=[]
        for p in [local,remote]:
            with np.load(p.parent/(sid+'.npz')) as z:arrays.append({k:z[k] for k in z.files})
        assert set(arrays[0])==set(arrays[1]);diff={}
        for k in arrays[0]:
            v,w=arrays[0][k],arrays[1][k];assert v.shape==w.shape
            if k=='sequence':
                assert str(v)==str(w) and hashlib.sha256(str(v).encode()).hexdigest()==sid[1:]
                continue
            assert np.isfinite(v).all() and np.isfinite(w).all()
            diff[k]={'mean_absolute_difference':float(np.mean(np.abs(v-w))),'maximum_absolute_difference':float(np.max(np.abs(v-w)))}
        rows.append({'sequence_id':sid,'length':len(x),'local_seconds':receipts[0]['inference_seconds'],'spark_seconds':receipts[1]['inference_seconds'],'ca_superposition_rmsd_angstrom':rmsd,'local_mean_plddt':receipts[0]['mean_ca_plddt'],'spark_mean_plddt':receipts[1]['mean_ca_plddt'],'array_differences':diff,'local_receipt_sha256':sha(local),'spark_receipt_sha256':sha(remote)})
    result={'status':'completed_four_sequence_platform_comparison','rows':rows,'source_input_sha256':sha(a.inputs),'script_sha256':sha(Path(__file__)),'spark_config_sha256':sha(a.predictions/'config.json'),'spark_chunk_sha256':sha(a.predictions/'last_chunk.json'),'interpretation':'Finite output and full-sequence checks on four length-selected proteins; independent artifact readback and rigid CA alignment. Not an accuracy assessment, complete platform equivalence proof or calibrated tolerance for evolutionary effects. Source/runtime differences remain explicit.'}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
