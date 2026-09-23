#!/usr/bin/env python3
"""Verify all exported domain atoms against original CIF arrays and manifest identities."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
import csv
import hashlib
import io
import json
import math
import multiprocessing
from pathlib import Path
import tarfile
import time
import psutil
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1
from catalog_whole_proteome_structures import sha


def verify_pdb(data,cif,sequence,start,end):
    names=['label_seq_id','label_atom_id','label_comp_id','Cartn_x','Cartn_y','Cartn_z','occupancy','B_iso_or_equiv','type_symbol']
    expected={}
    for pos,atom,comp,x,y,z,occ,b,element in zip(*(cif['_atom_site.'+n] for n in names)):
        pos=int(pos)
        if start<=pos<=end:
            key=(pos-start+1,atom)
            if key in expected:raise ValueError('Duplicate source domain atom')
            expected[key]=(comp,float(x),float(y),float(z),float(occ),float(b),element)
    seen=set();ca={};backbone={}
    for line in data.decode().splitlines():
        if not line.startswith('ATOM  '):
            if line not in ['TER','END']:raise ValueError('Unexpected PDB record')
            continue
        pos=int(line[22:26]);atom=line[12:16].strip();key=(pos,atom)
        if key not in expected or key in seen or line[21]!='A':raise ValueError('Exported atom identity differs')
        seen.add(key);backbone.setdefault(pos,set()).add(atom);src=expected[key]
        if line[17:20]!=src[0] or line[76:78].strip()!=src[-1]:raise ValueError('Residue or element differs')
        values=[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54),(54,60),(60,66)]]
        if not all(math.isfinite(v) for v in values):raise ValueError('Nonfinite PDB value')
        for observed,original,tolerance in zip(values,src[1:6],[0.000501]*3+[0.005001]*2):
            if abs(observed-original)>tolerance:raise ValueError('Exported coordinate, occupancy or confidence differs')
        if atom=='CA':ca[pos]=(protein_letters_3to1[src[0]],src[5])
    if seen!=set(expected) or set(ca)!=set(range(1,end-start+2)):raise ValueError('Incomplete exported atom/residue grid')
    fragment=''.join(ca[p][0] for p in sorted(ca))
    if fragment!=sequence[start-1:end]:raise ValueError('Exported fragment sequence differs')
    confidence=[ca[p][1] for p in sorted(ca)]
    return {'fragment_sequence_sha256':hashlib.sha256(fragment.encode()).hexdigest(),'residues':len(confidence),'mean_ca_plddt':sum(confidence)/len(confidence),'fraction_ca_plddt_ge70':sum(x>=70 for x in confidence)/len(confidence),'missing_backbone_source_positions':[p+start-1 for p in sorted(backbone) if not {'N','CA','C','O'}<=backbone[p]]}


def check_shard(job_path,shard_dir,proof):
    job_path=Path(job_path);root=Path(shard_dir);label=job_path.stem
    if sha(job_path)!=proof['job_sha256']:raise ValueError('Changed source job')
    for name,h in proof['artifacts'].items():
        if sha(root/name)!=h:raise ValueError('Changed extraction artifact')
    with (root/(label+'.jsonl')).open() as f:
        rows=[json.loads(line) for line in f]
    output={r['interval_id']:r for r in rows}
    if len(output)!=len(rows):raise ValueError('Repeated output disposition')
    counts=Counter();members=set()
    with tarfile.open(root/(label+'.tar'),'r') as archive,job_path.open() as jobs:
        tar_members=archive.getmembers();tar_names=[m.name for m in tar_members]
        if len(set(tar_names))!=len(tar_names) or any(not m.isfile() for m in tar_members):raise ValueError('Invalid archive membership')
        for line in jobs:
            job=json.loads(line);counts['models']+=1;model=job['model'];exported=[output[iid] for iid,_,_ in job['intervals'] if output[iid]['status']=='exported']
            if exported:
                data=Path(model['path']).read_bytes()
                if hashlib.sha256(data).hexdigest()!=model['sha256']:raise ValueError('Changed original CIF')
                cif=MMCIF2Dict(io.StringIO(data.decode()));seqs=[''.join(s.split()) for s in cif['_entity_poly.pdbx_seq_one_letter_code_can']]
                if len(seqs)!=1 or len(seqs[0])!=model['length'] or hashlib.sha256(seqs[0].encode()).hexdigest()!=model['sequence_sha256']:raise ValueError('Source polymer differs')
                sequence=seqs[0]
            for iid,start,end in job['intervals']:
                row=output.pop(iid);counts['intervals']+=1
                if (row['model_key'],row['start'],row['end'],row['source_sha256'])!=(job['model_key'],start,end,model['sha256']):raise ValueError('Disposition source binding differs')
                if row['status']=='exported':
                    if row['member']!=iid+'.pdb':raise ValueError('Wrong member name')
                    members.add(row['member']);blob=archive.extractfile(row['member']).read()
                    if hashlib.sha256(blob).hexdigest()!=row['pdb_sha256']:raise ValueError('PDB hash differs')
                    checked=verify_pdb(blob,cif,sequence,start,end)
                    for key,value in checked.items():
                        if isinstance(value,float):
                            if abs(value-row[key])>1e-10:raise ValueError('Confidence summary differs')
                        elif value!=row[key]:raise ValueError('Domain summary differs')
                    counts['exported']+=1;counts['missing_backbone_intervals']+=bool(checked['missing_backbone_source_positions'])
                elif row['status'] in ['source_validation_rejected','serialization_rejected'] and row.get('reason'):counts['rejected']+=1
                else:raise ValueError('Unknown or unexplained disposition')
        if output or set(tar_names)!=members:raise ValueError('Extra dispositions or archive members')
    for key in proof['counts']:
        if counts[key]!=proof['counts'][key]:raise ValueError('Shard totals differ')
    return dict(counts)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Plan changed')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Pinned source changed: '+p)
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    def state(stage,**kw):(out/'state.json').write_text(json.dumps(dict(stage=stage,**kw),indent=2)+'\n')
    state('waiting_for_coordinate_extraction')
    while True:
        try:
            p=psutil.Process(plan['predecessor']['pid']);live=p.create_time()==plan['predecessor']['create_time'] and p.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify();producer=json.loads(Path(plan['producer_plan']).read_text());root=Path(producer['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    if receipt['status']!='complete_domain_coordinate_dispositions_pending_independent_readback' or receipt['plan_sha256']!=sha(plan['producer_plan']):raise ValueError('Incomplete producer')
    expected={}
    with Path(producer['intervals']).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=row['interval_id']
            if key in expected:raise ValueError('Duplicate input interval')
            expected[key]=(row['model_key'],int(row['start']),int(row['end']),row['full_sequence_sha256'],row['model_path'])
    models={}
    with Path(producer['catalog_models']).open() as f:
        for line in f:
            m=json.loads(line);models[Path(m['path']).stem]={k:m[k] for k in ['path','sha256','sequence_sha256','length']}
    jobs=set()
    for proof in receipt['proofs']:
        job=Path(proof['job'])
        if str(job) in jobs or sha(job)!=proof['receipt']['job_sha256']:raise ValueError('Duplicate or changed job')
        jobs.add(str(job))
        with job.open() as f:
            for line in f:
                j=json.loads(line)
                if j['model']!=models[j['model_key']]:raise ValueError('Job/catalog model differs')
                for iid,start,end in j['intervals']:
                    if expected.pop(iid)!=(j['model_key'],start,end,j['model']['sequence_sha256'],j['model']['path']):raise ValueError('Job/manifest interval differs')
    if expected:raise ValueError('Incomplete extraction jobs')
    del expected,models
    counts=Counter();completed=0
    with ProcessPoolExecutor(max_workers=plan['resources']['cpu'],mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(check_shard,p['job'],str(root/'shards'),p['receipt']) for p in receipt['proofs']]
        for future in as_completed(futures):
            try:counts.update(future.result())
            except BaseException:
                for pending in futures:pending.cancel()
                raise
            completed+=1;state('checking_exported_atoms',shards=completed,counts=dict(counts))
    if any(counts[k]!=v for k,v in receipt['counts'].items()):raise ValueError('Global counts differ')
    verify();result={'status':'passed_all_exported_domain_atoms_and_full_disposition_scope','counts':dict(counts),'shards':completed,'producer_receipt_sha256':sha(rp),'plan_sha256':ph,'scope':'All original manifest/job identities checked; every exported atom, residue, coordinate, occupancy and confidence read back against original CIF arrays and fragment sequence. Rejection identities and recorded reasons accounted for but rejection causes not independently adjudicated. Shared CIF lexical parser; no PAE qualification or biological domain-boundary validation.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
