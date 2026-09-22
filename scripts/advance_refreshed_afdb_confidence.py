#!/usr/bin/env python3
"""Bind refreshed AFDB confidence to audited coordinates and check all contexts."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from assess_small_family_output_exposure import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());pins={str(a.plan):sha(a.plan),**plan['pins']}
    def verify():
        for path,digest in pins.items():
            if sha(path)!=digest:raise ValueError('Changed pinned dependency: '+path)
    verify();out=Path(plan['output'])
    if out.exists() or any(Path(s['output']).exists() for s in plan['stages']):raise FileExistsError('Inspect existing outputs before recovery')
    out.mkdir(parents=True);start=time.time();completed=[]
    def state(status,**details):
        temp=out/'state.tmp';temp.write_text(json.dumps(dict(status=status,completed_stages=completed,elapsed_seconds=time.time()-start,**details),indent=2)+'\n');temp.replace(out/'state.json')
    try:
        state('waiting_for_exact_prefetch_and_coordinate_controllers')
        for dependency in plan['dependencies']:
            proc=Path('/proc')/str(dependency['pid'])/'stat'
            while proc.exists():
                try:fields=proc.read_text().rsplit(')',1)[1].split()
                except FileNotFoundError:break
                if fields[19]!=str(dependency['start_ticks']) or fields[0]=='Z':break
                time.sleep(20)
            receipt=json.loads(Path(dependency['receipt']).read_text())
            if receipt['status']!=dependency['expected_status'] or receipt['plan_sha256']!=sha(dependency['plan']):
                raise ValueError('Dependency did not complete successfully: '+dependency['name'])
            pins[dependency['receipt']]=sha(dependency['receipt'])
        mapping_sha=sha(Path(plan['mapping'])/'receipt.json')
        feature=json.loads(Path(plan['dependencies'][0]['receipt']).read_text())
        prefetch=json.loads((Path(plan['prefetch'])/'receipt.json').read_text())
        pref_control=json.loads(Path(plan['dependencies'][1]['receipt']).read_text())
        if (feature['mapping_receipt_sha256']!=mapping_sha or feature['models']!=plan['models']
                or pref_control['result_receipt_sha256']!=sha(Path(plan['prefetch'])/'receipt.json')
                or prefetch['mapping_receipt_sha256']!=sha(Path(plan['catalog'])/'receipt.json')
                or prefetch['models_failed'] or prefetch['models_verified']!=plan['models']
                or prefetch['models_requested']!=plan['models']):
            raise ValueError('Feature or prefetch coverage/binding differs')
        pins[str(Path(plan['mapping'])/'receipt.json')]=mapping_sha
        pins[str(Path(plan['prefetch'])/'receipt.json')]=sha(Path(plan['prefetch'])/'receipt.json')
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
        for stage in plan['stages']:
            verify()
            if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
            command=[sys.executable,*stage['arguments']];state('running_stage',stage=stage['name'],command=command)
            with (out/(stage['name']+'.log')).open('w') as log:
                subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
            result=json.loads(Path(stage['result']).read_text())
            if any(result.get(k)!=v for k,v in stage['expected'].items()):raise ValueError('Stage scope/status differs')
            if 'mapping_receipt_sha256' in result and result['mapping_receipt_sha256']!=mapping_sha:raise ValueError('Wrong mapping')
            completed.append(dict(name=stage['name'],result=stage['result'],sha256=sha(stage['result'])))
        verify();final=json.loads(Path(plan['stages'][-1]['result']).read_text())
        receipt=dict(status='complete_refreshed_afdb_mapping_bound_confidence_and_full_context_readback',models=plan['models'],
                     mapping_receipt_sha256=mapping_sha,plan_sha256=sha(a.plan),stages=completed,totals=final['totals'],
                     scope='Full mapping-bound AFDB confidence qualification and independent ordered-pair context readback. No new prediction, confidence calibration, or evolutionary inference.')
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');state(receipt['status'])
    except Exception as exc:
        state('failed_requires_review',error=repr(exc));raise


if __name__=='__main__':main()
