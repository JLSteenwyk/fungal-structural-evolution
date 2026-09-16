#!/usr/bin/env python3
"""Compare an audited full-sequence control across two attention chunk settings."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from audit_joint_path_uncertainty import checked, sha


def load(folder, audit, sid):
    ar = checked(audit)
    audited = [r for r in csv.DictReader((audit/'predictions.tsv').open(), delimiter='\t') if r['sequence_id'] == sid]
    if len(audited) != 1 or audited[0]['prediction_receipt_sha256'] != sha(folder/(sid+'.json')):
        raise ValueError('Prediction not bound to audited receipt table')
    config = json.loads((folder/'config.json').read_text())
    if ar['config_sha256'] != sha(folder/'config.json'):
        raise ValueError('Audit configuration mismatch')
    receipt = json.loads((folder/(sid+'.json')).read_text())
    if receipt['config_sha256'] != ar['config_sha256'] or receipt['sequence_sha256'] != sid[1:] or receipt['status'] != 'verified_prediction':
        raise ValueError('Prediction identity mismatch')
    for name, digest in receipt['artifacts'].items():
        if sha(folder/name) != digest:
            raise ValueError('Changed prediction artifact')
    model = PDBParser(QUIET=True).get_structure(sid,folder/(sid+'.pdb'))
    xyz = np.array([r['CA'].coord for r in model.get_residues()], dtype=np.float64)
    # Read serialized decimal coordinates directly, avoiding parser float32 rounding.
    raw = np.array([[float(line[x:y]) for x,y in [(30,38),(38,46),(46,54)]] for line in (folder/(sid+'.pdb')).read_text().splitlines() if line.startswith('ATOM') and line[12:16].strip()=='CA'])
    np.testing.assert_allclose(xyz,raw,rtol=1e-6,atol=1e-5)
    with np.load(folder/(sid+'.npz'),allow_pickle=False) as z:
        sequence=str(z['sequence']); confidence=z['ca_plddt'].astype(np.float64); pae=z['pae'].astype(np.float64)
    if len(sequence)!=len(raw) or confidence.shape!=(len(raw),) or pae.shape!=(len(raw),len(raw)):
        raise ValueError('Shape mismatch')
    return config,receipt,sequence,raw,confidence,pae


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ['baseline','alternative','baseline-audit','alternative-audit','output']:
        ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--sequence-id',required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    left=load(a.baseline,a.baseline_audit,a.sequence_id)
    right=load(a.alternative,a.alternative_audit,a.sequence_id)
    changed={k:[left[0].get(k),right[0].get(k)] for k in set(left[0])|set(right[0]) if left[0].get(k)!=right[0].get(k)}
    if set(changed)-{'script_sha256','input_receipt_sha256','attention_chunk_size','max_length'}:
        raise ValueError('Unexpected configuration differences')
    if left[2]!=right[2]:raise ValueError('Different sequences')
    rows=[]
    for cutoff in [0,70,90]:
        keep=(left[4]>=cutoff)&(right[4]>=cutoff)
        x,y=left[3][keep],right[3][keep]
        if len(x)<3:continue
        xc,yc=x-x.mean(0),y-y.mean(0)
        u,_,vt=np.linalg.svd(xc.T@yc)
        correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt)
        rotation=u@correction@vt
        residual=np.linalg.norm(xc@rotation-yc,axis=1)
        rmsd=float(np.sqrt(np.mean(residual**2)))
        independent=SVDSuperimposer();independent.set(y,x);independent.run()
        np.testing.assert_allclose(rmsd,independent.get_rms(),rtol=1e-10,atol=1e-10)
        rows.append({'joint_native_ca_plddt_cutoff':cutoff,'residues':int(keep.sum()),'ca_rmsd_angstrom':rmsd,'max_aligned_ca_displacement_angstrom':float(residual.max())})
    delta_conf=right[4]-left[4];delta_pae=right[5]-left[5]
    result={'status':'complete_chunk_setting_control_comparison','sequence_id':a.sequence_id,'length':len(left[2]),'configuration_differences':changed,'geometry':rows,'serialized_ca_coordinates_exactly_equal':bool(np.array_equal(left[3],right[3])),'confidence_mean_absolute_difference':float(abs(delta_conf).mean()),'confidence_max_absolute_difference':float(abs(delta_conf).max()),'pae_mean_absolute_difference_angstrom':float(abs(delta_pae).mean()),'pae_max_absolute_difference_angstrom':float(abs(delta_pae).max()),'confidence_threshold_crossings':{str(c):int(np.sum((left[4]>=c)!=(right[4]>=c))) for c in [50,70,90]},'baseline_inference_seconds':left[1]['inference_seconds'],'alternative_inference_seconds':right[1]['inference_seconds'],'baseline_peak_allocated_bytes':left[1]['peak_gpu_allocated_bytes'],'alternative_peak_allocated_bytes':right[1]['peak_gpu_allocated_bytes'],'pins':{str(p):sha(p) for p in [a.baseline/'config.json',a.alternative/'config.json',a.baseline/(a.sequence_id+'.json'),a.alternative/(a.sequence_id+'.json'),a.baseline_audit/'receipt.json',a.alternative_audit/'receipt.json',Path(__file__)]},'interpretation':'Single full-sequence numerical sensitivity control; coordinate RMSD cross-checked with independent Bio.SVDSuperimposer. Native confidence masks may differ from rounded PDB masks. PAE differences include all directed entries including diagonal. No experimental accuracy, equivalence across proteins, or general memory-scaling claim.'}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
