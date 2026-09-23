#!/usr/bin/env python3
"""Extract all verified domain intervals into bounded, independently readable tar shards."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import csv
import hashlib
import io
import json
import math
import multiprocessing
from pathlib import Path
import shutil
import tarfile
import time
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1
from catalog_whole_proteome_structures import sha


def load_atoms(model):
    path=Path(model['path'])
    blob=path.read_bytes()
    if hashlib.sha256(blob).hexdigest()!=model['sha256']:raise ValueError('Source coordinate checksum mismatch')
    cif=MMCIF2Dict(io.StringIO(blob.decode()));polymers=[''.join(s.split()) for s in cif['_entity_poly.pdbx_seq_one_letter_code_can']]
    if len(polymers)!=1:raise ValueError('Expected one polymer entity')
    sequence=polymers[0]
    if len(sequence)!=model['length'] or hashlib.sha256(sequence.encode()).hexdigest()!=model['sequence_sha256']:raise ValueError('Source polymer sequence mismatch')
    keys=['group_PDB','label_atom_id','label_comp_id','label_seq_id','label_asym_id','label_alt_id','pdbx_PDB_model_num','Cartn_x','Cartn_y','Cartn_z','occupancy','B_iso_or_equiv','type_symbol']
    arrays=[cif['_atom_site.'+k] for k in keys]
    if len({len(a) for a in arrays})!=1:raise ValueError('Atom column lengths differ')
    residues={};seen=set();chains=set();ca={}
    for group,atom,comp,pos,chain,alt,number,x,y,z,occ,b,element in zip(*arrays):
        if group!='ATOM' or alt not in ('.','?') or number!='1':raise ValueError('Unsupported atom, alternate location or model')
        position=int(pos);chains.add(chain)
        if not 1<=position<=len(sequence) or protein_letters_3to1.get(comp)!=sequence[position-1]:raise ValueError('Atom residue identity mismatch or noncanonical residue')
        if (position,atom) in seen:raise ValueError('Repeated atom')
        seen.add((position,atom));values=list(map(float,(x,y,z,occ,b)))
        if not all(map(math.isfinite,values)) or not 0<=values[3]<=1 or not 0<=values[4]<=100:raise ValueError('Invalid coordinates, occupancy or confidence')
        if len(atom)>4 or len(element)>2:raise ValueError('Unsupported atom name')
        residues.setdefault(position,[]).append((atom,comp,*values,element))
        if atom=='CA':ca[position]=values[4]
    if len(chains)!=1 or set(ca)!=set(range(1,len(sequence)+1)):raise ValueError('Incomplete C-alpha coverage or multiple chains')
    return sequence,residues,ca


def pdb_bytes(sequence,residues,start,end):
    lines=[];serial=0;ca_sequence=[];missing_backbone=[]
    for original in range(start,end+1):
        position=original-start+1;atoms=residues.get(original,[]);names={a[0] for a in atoms}
        if not {'N','CA','C','O'}<=names:missing_backbone.append(original)
        for atom,comp,x,y,z,occupancy,b,element in atoms:
            serial+=1;name=atom if len(atom)==4 else ' '+atom.ljust(3)
            line=f'ATOM  {serial:5d} {name} {comp:3s} A{position:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{b:6.2f}          {element:>2s}  '
            if len(line)!=80:raise ValueError('PDB field overflow')
            # Independently decode the serialized fields to check rounding/identity.
            if int(line[22:26])!=position or line[12:16].strip()!=atom or line[17:20]!=comp:raise ValueError('PDB identity serialization mismatch')
            if any(abs(float(line[a:bnd])-v)>0.000501 for a,bnd,v in [(30,38,x),(38,46,y),(46,54,z)]):raise ValueError('PDB coordinate serialization mismatch')
            if abs(float(line[60:66])-b)>0.005001:raise ValueError('PDB confidence serialization mismatch')
            if atom=='CA':ca_sequence.append(protein_letters_3to1[comp])
            lines.append(line+'\n')
    fragment=sequence[start-1:end]
    if ''.join(ca_sequence)!=fragment:raise ValueError('Exported domain sequence mismatch')
    return (''.join(lines)+'TER\nEND\n').encode(),missing_backbone


def shard(job_path,out_dir,reserve_gib):
    job_path=Path(job_path);out=Path(out_dir);label=job_path.stem
    tarpath=out/(label+'.tar');manifest=out/(label+'.jsonl');receiptpath=out/(label+'.receipt.json')
    job_hash=sha(job_path)
    if receiptpath.exists():
        r=json.loads(receiptpath.read_text())
        if r['job_sha256']!=job_hash or any(sha(out/name)!=h for name,h in r['artifacts'].items()):raise ValueError('Changed completed shard')
        return r
    if tarpath.exists() or manifest.exists():raise FileExistsError('Unreceipted shard requires explicit recovery review')
    counts={'models':0,'intervals':0,'exported':0,'rejected':0,'missing_backbone_intervals':0}
    with tarfile.open(tarpath,'w',format=tarfile.USTAR_FORMAT) as archive,manifest.open('w') as log,job_path.open() as jobs:
        for line in jobs:
            if shutil.disk_usage(out).free<reserve_gib*2**30:raise RuntimeError('Disk reserve reached')
            job=json.loads(line);model=job['model'];counts['models']+=1
            try:sequence,residues,ca=load_atoms(model);problem=None
            except (ValueError,KeyError) as e:problem=f'{type(e).__name__}: {e}'
            for interval in job['intervals']:
                iid,start,end=interval;counts['intervals']+=1
                row={'interval_id':iid,'model_key':job['model_key'],'start':start,'end':end,'source_sha256':model['sha256']}
                if problem:
                    row.update(status='source_validation_rejected',reason=problem);counts['rejected']+=1
                else:
                    try:
                        data,missing=pdb_bytes(sequence,residues,start,end)
                        member=tarfile.TarInfo(iid+'.pdb');member.size=len(data);member.mtime=0;member.mode=0o644;archive.addfile(member,io.BytesIO(data))
                        confidence=[ca[p] for p in range(start,end+1)]
                        row.update(status='exported',member=member.name,pdb_sha256=hashlib.sha256(data).hexdigest(),fragment_sequence_sha256=hashlib.sha256(sequence[start-1:end].encode()).hexdigest(),residues=end-start+1,mean_ca_plddt=sum(confidence)/len(confidence),fraction_ca_plddt_ge70=sum(v>=70 for v in confidence)/len(confidence),missing_backbone_source_positions=missing)
                        counts['exported']+=1;counts['missing_backbone_intervals']+=bool(missing)
                    except ValueError as e:row.update(status='serialization_rejected',reason=str(e));counts['rejected']+=1
                log.write(json.dumps(row,sort_keys=True)+'\n')
            log.flush()
    result={'job_sha256':job_hash,'counts':counts,'artifacts':{p.name:sha(p) for p in [tarpath,manifest]}}
    receiptpath.write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Plan changed')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Pinned source changed: '+p)
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
    out.mkdir(parents=True);jobs=out/'jobs';jobs.mkdir();shards=out/'shards';shards.mkdir();groups={}
    with Path(plan['intervals']).open() as f:
        for r in csv.DictReader(f,delimiter='\t'):groups.setdefault(r['model_key'],[]).append((r['interval_id'],int(r['start']),int(r['end'])))
    models={}
    with Path(plan['catalog_models']).open() as f:
        for line in f:
            m=json.loads(line);key=Path(m['path']).stem
            if key in groups:models[key]={k:m[k] for k in ['path','sha256','sequence_sha256','length']}
    if set(models)!=set(groups):raise ValueError('Model coverage differs')
    if len(models)!=plan['expected_models'] or sum(map(len,groups.values()))!=plan['expected_intervals']:raise ValueError('Manifest dimensions differ from verified plan')
    keys=sorted(groups);job_paths=[]
    for index,start in enumerate(range(0,len(keys),plan['models_per_shard'])):
        path=jobs/f'shard_{index:05d}.jsonl'
        with path.open('w') as f:
            for key in keys[start:start+plan['models_per_shard']]:f.write(json.dumps({'model_key':key,'model':models[key],'intervals':groups[key]},separators=(',',':'))+'\n')
        job_paths.append(path)
    expected_models=len(models);expected_intervals=sum(map(len,groups.values()));del models,groups
    counts={k:0 for k in ['models','intervals','exported','rejected','missing_backbone_intervals']};started=time.monotonic();proofs=[]
    with ProcessPoolExecutor(max_workers=plan['resources']['cpu'],mp_context=multiprocessing.get_context('spawn')) as pool:
        futures={pool.submit(shard,str(path),str(shards),plan['resources']['emergency_free_disk_gib']):path for path in job_paths}
        for future in as_completed(futures):
            try:r=future.result()
            except BaseException:
                for pending in futures:pending.cancel()
                raise
            proofs.append({'job':str(futures[future]),'receipt':r})
            for key,value in r['counts'].items():counts[key]+=value
            (out/'state.json').write_text(json.dumps({'stage':'extracting','completed_shards':len(proofs),'total_shards':len(job_paths),'counts':counts,'elapsed_seconds':time.monotonic()-started},indent=2)+'\n')
    if counts['models']!=expected_models or counts['intervals']!=expected_intervals or counts['exported']+counts['rejected']!=expected_intervals:raise ValueError('Extraction scope differs')
    verify();result={'status':'complete_domain_coordinate_dispositions_pending_independent_readback','counts':counts,'shards':len(proofs),'plan_sha256':ph,'proofs':proofs,'scope':'Full manifest attempted; source hashes, full polymer and atom residue identity checked. All-atom PDB spans retain source coordinates with standard rounding, renumbered from one; serialized sequence/coordinate/confidence fields checked. Rejections and missing backbone retained explicitly. No PAE qualification, validated boundaries, homology or evolutionary inference. Independent archive readback required.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
