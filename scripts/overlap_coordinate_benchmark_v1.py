"""Exact-sequence predictor geometry on every original feature-confidence mask."""
import hashlib
import json
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from compare_marker_structures import geometry
from compare_predictor_alphabets import compare_states
from extract_domain_coordinates import load_atoms

METRICS = ['ca_superposition_rmsd_angstrom', 'local_distance_pairs',
           'local_distance_mean_absolute_change_angstrom', 'local_distance_rms_change_angstrom',
           'all_distance_pairs', 'all_distance_mean_absolute_change_angstrom',
           'all_distance_rms_change_angstrom', 'sequence_local_distance_pairs',
           'sequence_local_distance_mean_absolute_change_angstrom', 'sequence_local_distance_rms_change_angstrom']


def numeric_digest(value):
    return hashlib.sha256(value.dtype.str.encode()+str(value.shape).encode()+value.tobytes()).hexdigest()


def measure(x, y, positions):
    x,y,positions=np.asarray(x),np.asarray(y),np.asarray(positions)
    assert x.dtype==y.dtype==np.dtype('float64')
    assert x.shape==y.shape and x.shape==(len(positions),3) and len(positions)>=3
    assert np.isfinite(x).all() and np.isfinite(y).all()
    assert positions.ndim==1 and positions.dtype.kind in 'iu'
    assert np.all(positions>0) and np.all(np.diff(positions)>0)
    before=[numeric_digest(a) for a in [x,y,positions]]
    result=geometry(x,y,positions,positions)
    ii,jj=np.triu_indices(len(x),1)
    dx=np.linalg.norm(x[ii]-x[jj],axis=1);dy=np.linalg.norm(y[ii]-y[jj],axis=1)
    delta=dx-dy;gap=positions[jj]-positions[ii]
    for label,mask in [('all',np.ones(len(delta),dtype=bool)),('sequence_local',(gap>=3)&(gap<=10))]:
        selected=delta[mask]
        result[label+'_distance_pairs']=int(len(selected))
        result[label+'_distance_mean_absolute_change_angstrom']=float(np.mean(np.abs(selected))) if len(selected) else ''
        result[label+'_distance_rms_change_angstrom']=float(np.sqrt(np.mean(selected**2))) if len(selected) else ''
    assert before==[numeric_digest(a) for a in [x,y,positions]],'Input coordinate memory changed'
    assert set(result)==set(METRICS)
    assert all(np.isfinite(v) and v>=0 for v in result.values() if v!='')
    return result


def load_model(row, prefix, pins):
    ep=Path(row[prefix+'_encoding_path']);cp=Path(row[prefix+'_model_path'])
    assert sha(ep)==row[prefix+'_encoding_sha256'] and sha(cp)==row[prefix+'_model_sha256']
    pins[str(ep)]=row[prefix+'_encoding_sha256'];pins[str(cp)]=row[prefix+'_model_sha256']
    with np.load(ep,allow_pickle=False) as f:encoded={k:f[k].copy() for k in f.files}
    encoded['sequence']=str(encoded['sequence']);encoded['states_array']=np.asarray(list(str(encoded['states'])))
    assert hashlib.sha256(encoded['sequence'].encode()).hexdigest()==row[prefix+'_sequence_sha256']
    n=len(encoded['sequence']);assert n==int(row[prefix+'_length'])
    sequence,residues,plddt=load_atoms(dict(path=str(cp),sha256=row[prefix+'_model_sha256'],
        sequence_sha256=row[prefix+'_sequence_sha256'],length=n))
    assert sequence==encoded['sequence']
    xyz=np.asarray([next(atom[2:5] for atom in residues[i] if atom[0]=='CA') for i in range(1,n+1)],dtype=np.float64)
    assert np.array_equal(np.asarray([plddt[i] for i in range(1,n+1)]),encoded['ca_plddt'])
    assert xyz.shape==(n,3) and np.isfinite(xyz).all()
    xyz.flags.writeable=False
    return encoded,xyz


def identity(row):
    return tuple(row[k] for k in ['reference_model_id','reference_model_version','local_model_id','local_model_version'])


def pair_id(key):
    return hashlib.sha256(json.dumps(key,separators=(',',':')).encode()).hexdigest()
