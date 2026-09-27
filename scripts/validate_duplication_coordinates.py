#!/usr/bin/env python3
"""Checkpoint raw-coordinate validation for every model in the reviewed pair queue."""
import argparse,csv,fcntl,gzip,hashlib,json,math,os,shutil,time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from extract_domain_coordinates import load_atoms
from run_ortholog_pair_guide_comparison import sha


def validate_model(model):
    sequence,residues,confidence=load_atoms(model)
    values=[confidence[i] for i in range(1,len(sequence)+1)]
    mean=sum(values)/len(values);below=sum(x<50 for x in values)/len(values)
    if not math.isclose(mean,model['mean_ca_plddt'],rel_tol=0,abs_tol=1e-9) or not math.isclose(below,model['fraction_ca_plddt_below50'],rel_tol=0,abs_tol=1e-12):raise ValueError('Catalog confidence summary differs')
    coords=[]
    for i in range(1,len(sequence)+1):
        ca=[a for a in residues[i] if a[0]=='CA']
        if len(ca)!=1:raise ValueError('Ambiguous C-alpha coordinate')
        coords.append(list(ca[0][2:5]))
    return dict(status='validated',sequence=sequence,ca_xyz=coords,ca_plddt=values,length=len(sequence),mean_ca_plddt=mean,fraction_ca_plddt_below50=below,residues_ge70=sum(x>=70 for x in values),residues_ge90=sum(x>=90 for x in values))


def run_shard(job):
    index,models,out,ph,reserve=job;out=Path(out);dest=out/f'shard-{index:05d}.jsonl.gz';rp=out/f'shard-{index:05d}.receipt.json'
    mh=hashlib.sha256(json.dumps(models,sort_keys=True).encode()).hexdigest()
    if rp.exists():
        r=json.loads(rp.read_text())
        if r['plan_sha256']!=ph or r['models_sha256']!=mh or sha(dest)!=r['output_sha256']:raise ValueError('Checkpoint mismatch')
        # Reusing a checkpoint still requires unchanged raw input bytes.
        for m in models:
            if sha(m['path'])!=m['sha256']:raise ValueError('Checkpoint source changed')
        return r
    if dest.exists() or dest.with_suffix('.partial').exists():raise FileExistsError('Unreceipted shard requires review')
    if shutil.disk_usage(out).free<reserve*2**30:raise RuntimeError('Disk reserve reached')
    counts=Counter();start=time.monotonic();temp=dest.with_suffix('.partial')
    with gzip.open(temp,'wt') as f:
        for m in models:
            row={'model_id':m['model_id'],'version':m['version'],'source_path':m['path'],'source_sha256':m['sha256'],'sequence_sha256':m['sequence_sha256']}
            # A changed/missing file is a provenance failure, not an excluded model.
            if sha(m['path'])!=m['sha256']:raise ValueError('Source coordinate changed')
            try:row.update(validate_model(m))
            except (ValueError,KeyError) as e:row.update(status='rejected_content',reason=f'{type(e).__name__}: {e}')
            counts[row['status']]+=1;f.write(json.dumps(row,separators=(',',':'))+'\n')
    temp.replace(dest);r=dict(plan_sha256=ph,models_sha256=mh,models=len(models),counts=dict(counts),output=dest.name,output_sha256=sha(dest),elapsed_seconds=time.monotonic()-start)
    rp.write_text(json.dumps(r,indent=2)+'\n');return r


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed pin '+path)
    verify();queue=Path(plan['queue']);qr=json.loads((queue/'receipt.json').read_text());source=queue/'models.jsonl'
    if qr['status']!='complete_reviewed_duplication_model_pair_queue' or sha(source)!=qr['artifacts']['models.jsonl']:raise ValueError('Queue binding mismatch')
    with source.open() as f:models=[json.loads(l) for l in f]
    if len(models)!=qr['unique_models'] or len({(m['model_id'],m['version']) for m in models})!=len(models):raise ValueError('Model universe differs')
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True);lock=(out/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    count=plan['models_per_shard'];jobs=[(i,models[start:start+count],str(out),ph,plan['minimum_free_disk_gib']) for i,start in enumerate(range(0,len(models),count))]
    totals=Counter();receipts=[]
    with ProcessPoolExecutor(max_workers=plan['workers']) as pool:
        for r in pool.map(run_shard,jobs):
            receipts.append(r);totals.update(r['counts'])
            state=dict(stage='validating_coordinates',completed_shards=len(receipts),total_shards=len(jobs),counts=dict(totals));(out/'state.json').write_text(json.dumps(state)+'\n');print(json.dumps(state),flush=True)
    verify();assert sum(totals.values())==len(models)
    result=dict(status='complete_duplication_coordinate_validation_with_dispositions',plan_sha256=ph,models=len(models),counts=dict(totals),shards=receipts,script_sha256=sha(__file__),scope='Every frozen candidate model including identical-model links. Raw byte hash, canonical polymer/atom identities, single-model/chain constraints, C-alpha coverage, finite all-atom values and catalog confidence summaries checked. Source/provenance errors fail; content rejections retained. Residue-level C-alpha coordinates and confidence retained for later alignment and masks. No PAE qualification, structural comparison, biological duplication or asymmetry inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
