#!/usr/bin/env python3
"""Build and read back a domain search database after full original-atom audit."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time
import numpy as np
import psutil
from build_whole_proteome_foldseek_database import readback,index
from catalog_whole_proteome_structures import sha


def compare_all_coordinates(prefix,archives,expected):
    lookup={}
    with Path(str(prefix)+'.lookup').open() as f:
        for line in f:
            key,name,_=line.rstrip('\n').split('\t')
            if name in lookup:raise ValueError('Repeated database name')
            lookup[name]=int(key)
    indices=index(Path(str(prefix)+'_ca.index'));seen=set();residues=0
    with Path(str(prefix)+'_ca').open('rb') as handle:
        for path in archives:
            with tarfile.open(path,'r') as archive:
                for member in archive:
                    name=Path(member.name).stem
                    if name not in expected or name in seen or not member.isfile():raise ValueError('Unexpected archive member')
                    seen.add(name);data=archive.extractfile(member).read()
                    coordinates=np.array([[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]] for line in data.decode().splitlines() if line.startswith('ATOM  ') and line[12:16].strip()=='CA'],dtype=np.float32)
                    offset,length=indices[lookup[name]];handle.seek(offset);raw=handle.read(length)
                    if not raw.endswith(b'\0'):raise ValueError('Invalid coordinate terminator')
                    actual=np.frombuffer(raw[:-1],dtype=np.float32).reshape(3,-1).T
                    if coordinates.shape!=(expected[name],3) or not np.array_equal(actual,coordinates):raise ValueError('Native coordinates differ from audited domain PDB')
                    residues+=len(coordinates)
    if seen!=set(expected):raise ValueError('Incomplete coordinate readback')
    return residues


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Plan changed')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed bound source: '+p)
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    def state(status,**kw):(out/'state.json').write_text(json.dumps(dict(status=status,**kw),indent=2)+'\n')
    state('waiting_for_full_domain_atom_audit')
    while True:
        try:
            p=psutil.Process(plan['predecessor']['pid']);live=p.create_time()==plan['predecessor']['create_time'] and p.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify();audit_path=Path(plan['audit'])/'receipt.json';audit=json.loads(audit_path.read_text());extraction=Path(plan['extraction']);receipt_path=extraction/'receipt.json';receipt=json.loads(receipt_path.read_text())
    if audit['status']!='passed_all_exported_domain_atoms_and_full_disposition_scope' or audit['producer_receipt_sha256']!=sha(receipt_path) or audit['plan_sha256']!=sha(plan['audit_plan']) or receipt['plan_sha256']!=sha(plan['extraction_plan']):raise ValueError('Incomplete or inconsistent original-atom audit')
    if any(audit['counts'].get(k,0)!=v for k,v in receipt['counts'].items()):raise ValueError('Audit/extraction scope differs')
    if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30 or psutil.virtual_memory().available<plan['resources']['minimum_available_memory_gib']*2**30:raise ValueError('Insufficient resources')
    archives=[];models=[];expected={};shards=extraction/'shards';source_hashes={}
    for proof in sorted(receipt['proofs'],key=lambda p:p['job']):
        label=Path(proof['job']).stem
        for name,h in proof['receipt']['artifacts'].items():
            if sha(shards/name)!=h:raise ValueError('Changed validated archive/disposition')
            source_hashes[str(shards/name)]=h
        archives.append(shards/(label+'.tar'))
        with (shards/(label+'.jsonl')).open() as f:
            for line in f:
                row=json.loads(line)
                if row['status']!='exported':continue
                iid=row['interval_id']
                if iid in expected or row['member']!=iid+'.pdb':raise ValueError('Duplicate domain or member mismatch')
                expected[iid]=row['residues'];models.append({'path':row['member'],'length':row['residues'],'sequence_sha256':row['fragment_sequence_sha256']})
    if len(models)!=audit['counts']['exported']:raise ValueError('Exported domain count differs')
    paths=out/'archives.tsv';paths.write_text(''.join(str(p.resolve())+'\n' for p in archives));prefix=out/'domains'
    command=[plan['foldseek'],'createdb',str(paths),str(prefix),'--threads',str(plan['resources']['cpu']),'--gpu','0','--mask-bfactor-threshold','70','--coord-store-mode','1']
    state('building_domain_database',command=command,models=len(models))
    with (out/'createdb.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,env=dict(os.environ,CUDA_VISIBLE_DEVICES=''))
    state('checking_every_domain_sequence_and_coordinate',models=len(models));checked=readback(prefix,models)
    coordinate_residues=compare_all_coordinates(prefix,archives,expected)
    if coordinate_residues!=checked['residues']:raise ValueError('Sequence/coordinate residue count differs')
    for path,h in source_hashes.items():
        if sha(path)!=h:raise ValueError('Archive changed during conversion')
    verify();result={'status':'complete_domain_foldseek_database_with_full_sequence_and_coordinate_readback',**checked,'plan_sha256':ph,'extraction_receipt_sha256':sha(receipt_path),'audit_receipt_sha256':sha(audit_path),'excluded_rejected_intervals':receipt['counts']['rejected'],'exported_intervals_with_missing_backbone':receipt['counts']['missing_backbone_intervals'],'command':command,'artifacts':{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='state.json'},'scope':'All successfully exported intervals encoded afresh, retaining both boundary definitions. Every lookup/AA sequence, 3Di length/alphabet and C-alpha coordinate compared with audited PDB spans. pLDDT masking affects seeding, not full confidence qualification. Native 3Di states not independently reconstructed here; PAE, boundary sensitivity, clustering and evolutionary inference remain separate.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');state('complete',models=len(models))


if __name__=='__main__':main()
