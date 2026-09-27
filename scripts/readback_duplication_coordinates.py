#!/usr/bin/env python3
"""Independently compare all exported C-alpha records to their frozen raw CIFs."""
import argparse,gzip,hashlib,io,json,math,time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import psutil
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1
from run_ortholog_pair_guide_comparison import sha


def check_record(record,model,blob):
    for key,value in [('model_id',model['model_id']),('version',model['version']),('source_path',model['path']),('source_sha256',model['sha256']),('sequence_sha256',model['sequence_sha256'])]:
        if record[key]!=value:raise ValueError('Record identity differs: '+key)
    if hashlib.sha256(blob).hexdigest()!=model['sha256']:raise ValueError('Raw source bytes changed')
    if record['status']=='rejected_content':
        if not record.get('reason'):raise ValueError('Missing rejection reason')
        return 0
    if record['status']!='validated':raise ValueError('Unexpected disposition')
    cif=MMCIF2Dict(io.StringIO(blob.decode()))
    columns=[cif['_atom_site.'+k] for k in ['label_atom_id','label_seq_id','label_comp_id','label_asym_id','Cartn_x','Cartn_y','Cartn_z','B_iso_or_equiv']]
    if len({len(c) for c in columns})!=1:raise ValueError('Unequal CIF columns')
    ca={};chains=set()
    for atom,index,res,chain,x,y,z,b in zip(*columns):
        if atom!='CA':continue
        i=int(index)
        if i in ca:raise ValueError('Repeated source C-alpha')
        xyz=[float(x),float(y),float(z)];v=float(b)
        if not all(math.isfinite(n) for n in xyz+[v]) or not 0<=v<=100:raise ValueError('Invalid source value')
        ca[i]=(protein_letters_3to1[res],xyz,v);chains.add(chain)
    n=model['length']
    if set(ca)!=set(range(1,n+1)) or len(chains)!=1:raise ValueError('Source residue grid differs')
    sequence=''.join(ca[i][0] for i in range(1,n+1));xyz=[ca[i][1] for i in range(1,n+1)];confidence=[ca[i][2] for i in range(1,n+1)]
    if sequence!=record['sequence'] or hashlib.sha256(sequence.encode()).hexdigest()!=model['sequence_sha256']:raise ValueError('Exported sequence differs')
    if record['ca_xyz']!=xyz or record['ca_plddt']!=confidence or record['length']!=n:raise ValueError('Exported numeric values differ')
    for threshold in [70,90]:
        if record['residues_ge'+str(threshold)]!=sum(v>=threshold for v in confidence):raise ValueError('Confidence threshold count differs')
    mean=sum(confidence)/n;fraction=sum(v<50 for v in confidence)/n
    for name,v in [('mean_ca_plddt',mean),('fraction_ca_plddt_below50',fraction)]:
        if not math.isclose(record[name],v,rel_tol=0,abs_tol=1e-12) or not math.isclose(model[name],v,rel_tol=0,abs_tol=1e-9):raise ValueError('Confidence summary differs')
    return n


def check_shard(job):
    folder,receipt,models,output,ph=job;folder=Path(folder);out=Path(output);dest=out/(receipt['output']+'.readback.json')
    source=folder/receipt['output'];before=sha(source)
    if before!=receipt['output_sha256']:raise ValueError('Changed shard output')
    counts=Counter();residues=0
    with gzip.open(source,'rt') as f:
        for index,model in enumerate(models):
            line=f.readline()
            if not line:raise ValueError('Missing exported model')
            row=json.loads(line);residues+=check_record(row,model,Path(model['path']).read_bytes());counts[row['status']]+=1
        if f.readline():raise ValueError('Extra exported model')
    if dict(counts)!=receipt['counts'] or len(models)!=receipt['models']:raise ValueError('Shard count mismatch')
    if sha(source)!=before:raise ValueError('Shard changed during readback')
    r=dict(plan_sha256=ph,source_sha256=before,models=len(models),counts=dict(counts),validated_residues=residues)
    if dest.exists():
        if json.loads(dest.read_text())!=r:raise ValueError('Prior proof differs')
    else:dest.write_text(json.dumps(r,indent=2)+'\n')
    return r


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed audit plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed audit pin '+path)
    verify();producer=plan['producer']
    while psutil.pid_exists(producer['pid']):
        try:
            proc=psutil.Process(producer['pid'])
            if abs(proc.create_time()-producer['created'])>0.01 or proc.status()==psutil.STATUS_ZOMBIE:break
            if proc.cmdline()!=producer['cmdline']:raise ValueError('Producer identity changed')
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();folder=Path(plan['source']);rp=folder/'receipt.json';r=json.loads(rp.read_text());sourceplan=json.loads(Path(plan['source_plan']).read_text());producer_hash=sha(rp)
    if r['status']!='complete_duplication_coordinate_validation_with_dispositions' or r['plan_sha256']!=sha(plan['source_plan']):raise ValueError('Incomplete coordinate producer')
    with Path(plan['models']).open() as f:models=[json.loads(l) for l in f]
    n=sourceplan['models_per_shard'];expected=(len(models)+n-1)//n
    if len(r['shards'])!=expected or r['models']!=len(models):raise ValueError('Incomplete shard/model universe')
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);totals=Counter();residues=0
    jobs=[]
    for i,receipt in enumerate(r['shards']):
        subset=models[i*n:(i+1)*n];mh=hashlib.sha256(json.dumps(subset,sort_keys=True).encode()).hexdigest()
        if receipt['models_sha256']!=mh or receipt['output']!=f'shard-{i:05d}.jsonl.gz':raise ValueError('Shard source mapping differs')
        jobs.append((str(folder),receipt,subset,str(out),ph))
    with ProcessPoolExecutor(max_workers=plan['workers']) as pool:
        for i,proof in enumerate(pool.map(check_shard,jobs),1):
            totals.update(proof['counts']);residues+=proof['validated_residues'];print('Checked shards',i,'/',expected,flush=True)
    verify()
    if sha(rp)!=producer_hash:raise ValueError('Producer receipt changed during audit')
    if dict(totals)!=r['counts']:raise ValueError('Full disposition counts differ')
    result=dict(status='passed_duplication_exported_ca_readback',plan_sha256=ph,producer_receipt_sha256=sha(rp),models=len(models),counts=dict(totals),validated_residues=residues,proofs={p.name:sha(p) for p in out.iterdir()},scope='Every exported accepted C-alpha sequence, coordinate and confidence value independently reconstructed from frozen CIF atom rows; full model/disposition grid and summaries checked. Shares CIF lexical parser, not producer extraction function. Rejection identities/reasons retained but rejection causes not independently adjudicated. No PAE or biological inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
