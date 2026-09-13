#!/usr/bin/env python3
"""Bind locally predicted PAE to an exact completed structure-mapping snapshot."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT,sha
from retrieve_marker_pae import validate_pae


def encode_pae(sequence, pae, maximum, model):
    if len(sequence)!=model['length'] or hashlib.sha256(sequence.encode()).hexdigest()!=model['sequence_sha256']:
        raise ValueError('Local PAE sequence identity differs')
    if pae.shape!=(len(sequence),len(sequence)) or not np.isfinite(pae).all() or (pae<0).any() or not np.isfinite(maximum) or maximum<=0 or (pae>maximum+1e-4).any():
        raise ValueError('Invalid local PAE matrix or declared maximum')
    raw=json.dumps([{'predicted_aligned_error':pae.tolist(),'max_predicted_aligned_error':maximum}],separators=(',',':'),allow_nan=False).encode()
    compressed=gzip.compress(raw,mtime=0)
    restored=validate_pae(gzip.decompress(compressed),len(sequence))
    if not np.array_equal(restored,pae):raise ValueError('Local PAE changed during export')
    return raw,compressed


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable local PAE output')
    source=checked_receipt(a.snapshot)
    if source['source_policy']['provider']!='local' or source['source_policy']['tool']!='ESMFold v1':raise ValueError('Local ESMFold mapping required')
    models=json.loads((a.snapshot/'model_provenance.json').read_text())
    keys={(m['model_id'],str(m['version'])) for m in models}
    if len(keys)!=len(models) or not models:raise ValueError('Missing or repeated mapped models')
    a.output=a.output.resolve();a.output.mkdir(parents=True);folder=a.output/'matrices';folder.mkdir()
    records=[]
    for m in models:
        if m['provider']!='local' or m['tool']!='ESMFold v1':raise ValueError('Mixed model source')
        pr=ROOT/m['prediction_receipt_path'];npz=ROOT/m['local_pae_npz_path']
        if sha(pr)!=m['prediction_receipt_sha256'] or sha(npz)!=m['local_pae_npz_sha256']:raise ValueError('Changed local prediction source')
        pred=json.loads(pr.read_text())
        if pred['status']!='verified_prediction' or pred['config_sha256']!=m['prediction_config_sha256'] or pred['sequence_sha256']!=m['sequence_sha256'] or pred['artifacts'][npz.name]!=m['local_pae_npz_sha256']:raise ValueError('Prediction source linkage differs')
        with np.load(npz,allow_pickle=False) as data:raw,compressed=encode_pae(str(data['sequence']),data['pae'],float(data['max_pae']),m)
        target=folder/(m['model_id']+'-v'+str(m['version'])+'.json.gz');target.write_bytes(compressed)
        records.append({'model_id':m['model_id'],'version':m['version'],'sequence_sha256':m['sequence_sha256'],'length':m['length'],
            'status':'verified','source_type':'local_prediction_export','source_npz_path':m['local_pae_npz_path'],'source_npz_sha256':m['local_pae_npz_sha256'],
            'prediction_config_sha256':m['prediction_config_sha256'],'prediction_receipt_sha256':m['prediction_receipt_sha256'],
            'path':str(target.relative_to(ROOT)),'gzip_sha256':sha(target),'json_sha256':hashlib.sha256(raw).hexdigest(),'compressed_bytes':len(compressed),'json_bytes':len(raw)})
        if len(records)%500==0:print('Exported',len(records),'of',len(models),flush=True)
    manifest=a.output/'pae_manifest.json';manifest.write_text(json.dumps(records,indent=2)+'\n')
    result={'status':'complete_local_mapping_bound_pae_export','mapping_snapshot':str(a.snapshot),'mapping_receipt_sha256':sha(a.snapshot/'receipt.json'),
        'models_requested':len(models),'models_verified':len(records),'models_failed':0,'source_type':'local_prediction_export',
        'compressed_bytes':sum(r['compressed_bytes'] for r in records),'json_bytes':sum(r['json_bytes'] for r in records),'script_sha256':sha(Path(__file__)),
        'interpretation':'Exact lossless export of local ESMFold directional PAE into the downstream matrix schema. No download, calibration against AlphaFold, symmetrization or inference of experimental accuracy.',
        'artifacts':{'pae_manifest.json':sha(manifest)}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
