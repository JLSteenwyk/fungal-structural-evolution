#!/usr/bin/env python3
"""Read back all background coordinate shards with the independently validated CIF checker."""
import argparse,hashlib,json,time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha
from readback_duplication_coordinates import check_shard

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
    if r['status']!='complete_background_coordinate_validation_with_dispositions' or r['plan_sha256']!=sha(plan['source_plan']):raise ValueError('Incomplete coordinate producer')
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
    result=dict(status='passed_background_exported_ca_readback',plan_sha256=ph,producer_receipt_sha256=sha(rp),models=len(models),counts=dict(totals),validated_residues=residues,proofs={p.name:sha(p) for p in out.iterdir()},scope='Every exported accepted C-alpha sequence, coordinate and confidence value independently reconstructed from frozen CIF atom rows; full model/disposition grid and summaries checked. Shares CIF lexical parser, not producer extraction function. Rejection identities/reasons retained but rejection causes not independently adjudicated. No PAE or biological inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
