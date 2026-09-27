#!/usr/bin/env python3
"""Validate every additional background model using checkpointed coordinate shards."""
import argparse,fcntl,json,time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha
from validate_duplication_coordinates import run_shard


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();dep=plan['dependency']
    while True:
        try:
            proc=psutil.Process(dep['pid'])
            if proc.create_time()!=dep['created'] or proc.status()==psutil.STATUS_ZOMBIE:break
            assert proc.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();root=Path(plan['inventory']);r=json.loads((root/'receipt.json').read_text());audit=json.loads(Path(plan['readback']).read_text())
    assert r['status']=='complete_background_measurement_inventory_pending_readback'
    assert audit['status']=='passed_full_background_measurement_inventory_readback' and audit['producer_receipt_sha256']==sha(root/'receipt.json')
    path=root/'additional_models.jsonl';assert sha(path)==r['artifacts'][path.name]
    models=[json.loads(line) for line in path.open()]
    assert len(models)==len({(m['model_id'],m['version']) for m in models})==r['additional_models']==audit['additional_models']
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    with (out/'run.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        n=plan['models_per_shard'];jobs=[(i,models[start:start+n],str(out),ph,plan['minimum_free_disk_gib']) for i,start in enumerate(range(0,len(models),n))]
        counts=Counter();receipts=[]
        with ProcessPoolExecutor(max_workers=plan['workers']) as pool:
            for receipt in pool.map(run_shard,jobs):
                receipts.append(receipt);counts.update(receipt['counts'])
                state=dict(stage='validating_background_coordinates',completed_shards=len(receipts),total_shards=len(jobs),counts=dict(counts));(out/'state.json').write_text(json.dumps(state)+'\n');print(json.dumps(state),flush=True)
        verify();assert sum(counts.values())==len(models)
        result=dict(status='complete_background_coordinate_validation_with_dispositions',plan_sha256=ph,models=len(models),counts=dict(counts),shards=receipts,source_inventory_receipt_sha256=sha(root/'receipt.json'),source_readback_sha256=sha(plan['readback']),scope='All additional frozen background models, including identical-model pairs. Raw hashes, polymer/atom identity, single model/chain, full C-alpha coverage, finite coordinates and confidence summaries checked with the existing pinned shard validator. Content rejections explicit; source errors fail. Independent coordinate readback, reuse of earlier model coordinates, PAE/coverage checks and alignments remain downstream. No new predictions.')
        (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
