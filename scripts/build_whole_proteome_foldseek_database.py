#!/usr/bin/env python3
"""Build a full structural-search database after independent catalog readback."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import hashlib
import numpy as np
import psutil
from catalog_whole_proteome_structures import sha


def index(path):
    result={}
    with path.open() as handle:
        for line in handle:
            key,offset,length=map(int,line.split())
            if key in result or offset<0 or length<=0:raise ValueError('Invalid database index')
            result[key]=(offset,length)
    return result


def readback(prefix,models):
    expected={Path(m['path']).stem:m for m in models}
    if len(expected)!=len(models):raise ValueError('Ambiguous input filename')
    lookup={}
    with Path(str(prefix)+'.lookup').open() as handle:
        for line in handle:
            key,name,_=line.rstrip('\n').split('\t');key=int(key)
            if key in lookup or name not in expected:raise ValueError('Unexpected database model')
            lookup[key]=name
    if len(lookup)!=len(expected) or set(lookup.values())!=set(expected):raise ValueError('Missing or repeated database model')
    aa=index(Path(str(prefix)+'.index'));ss=index(Path(str(prefix)+'_ss.index'));ca=index(Path(str(prefix)+'_ca.index'))
    if any(set(i)!=set(lookup) for i in (aa,ss,ca)):raise ValueError('Database key grids differ')
    total=0
    with prefix.open('rb') as amino,Path(str(prefix)+'_ss').open('rb') as structural,Path(str(prefix)+'_ca').open('rb') as coordinates:
        for key,name in lookup.items():
            model=expected[name];length=model['length'];total+=length
            buffers=[]
            for handle,indices in ((amino,aa),(structural,ss),(coordinates,ca)):
                offset,size=indices[key];handle.seek(offset);value=handle.read(size)
                if len(value)!=size or not value.endswith(b'\0'):raise ValueError('Invalid database record boundary')
                buffers.append(value)
            sequence=buffers[0].removesuffix(b'\0').removesuffix(b'\n').upper()
            states=buffers[1].removesuffix(b'\0').removesuffix(b'\n').upper()
            if len(sequence)!=length or hashlib.sha256(sequence).hexdigest()!=model['sequence_sha256']:
                raise ValueError('Database amino-acid sequence differs')
            if len(states)!=length or not set(states)<=set(b'ACDEFGHIKLMNPQRSTVWYX'):
                raise ValueError('Structural-alphabet length or symbols differ')
            raw=buffers[2][:-1]
            if len(raw)!=12*length or not np.isfinite(np.frombuffer(raw,dtype=np.float32)).all():
                raise ValueError('Coordinate dimensions or finiteness differ')
    return dict(models=len(models),residues=total)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan_sha=sha(args.plan)
    def verify():
        if sha(args.plan)!=plan_sha:raise ValueError('Changed plan')
        for path,digest in plan['pins'].items():
            if sha(path)!=digest:raise ValueError('Changed pinned input')
    verify();dependency=plan['predecessor']
    while True:
        try:
            proc=psutil.Process(dependency['pid'])
            live=proc.create_time()==dependency['create_time'] and proc.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify();catalog=Path(plan['catalog']);audit_path=Path(plan['readback'])
    receipt=json.loads((catalog/'receipt.json').read_text());audit=json.loads(audit_path.read_text())
    if (audit['status']!='passed_full_proteome_sequence_and_model_selection_readback'
            or audit['producer_receipt_sha256']!=sha(catalog/'receipt.json')
            or audit['plan_sha256']!=sha(plan['readback_plan'])
            or receipt['status']!='complete_whole_representative_proteome_exact_sequence_catalog'
            or receipt['plan_sha256']!=sha(plan['catalog_plan']) or audit['models']!=receipt['unique_models']):
        raise ValueError('Incomplete or inconsistent full-catalog readback')
    for name,digest in receipt['artifacts'].items():
        if sha(catalog/name)!=digest:raise ValueError('Changed catalog artifact')
    with (catalog/'models.jsonl').open() as handle:models=[json.loads(line) for line in handle]
    if len(models)!=audit['models']:raise ValueError('Model grid differs')
    out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
    if psutil.virtual_memory().available<plan['resources']['minimum_available_memory_gib']*2**30:raise ValueError('Insufficient available memory')
    out.mkdir();paths=out/'model_paths.tsv'
    paths.write_text(''.join(str(Path(m['path']).resolve())+'\n' for m in models))
    prefix=out/'structures';command=[plan['foldseek'],'createdb',str(paths),str(prefix),'--threads','4','--gpu','0','--mask-bfactor-threshold','70','--coord-store-mode','1']
    dimensions=dict(models=len(models),residues=sum(m['length'] for m in models),command=command,
                    coordinate_store='float32',catalog_receipt_sha256=sha(catalog/'receipt.json'),readback_sha256=sha(audit_path))
    (out/'prelaunch.json').write_text(json.dumps(dimensions,indent=2)+'\n')
    with (out/'createdb.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,env=dict(os.environ,CUDA_VISIBLE_DEVICES=''))
    checked=readback(prefix,models);verify()
    result=dict(status='complete_whole_proteome_foldseek_database_and_sequence_readback',**checked,
                plan_sha256=plan_sha,catalog_receipt_sha256=dimensions['catalog_receipt_sha256'],
                catalog_readback_sha256=dimensions['readback_sha256'],artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},
                scope='All selected models converted; every lookup identity, AA sequence hash, structural-alphabet length/symbol set and float-coordinate shape/finiteness checked. Mask pLDDT70 affects seeding, not full confidence qualification. 3Di values and coordinates not independently reconstructed from CIF here. Search, clustering, domain and homology inference remain separate.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
