#!/usr/bin/env python3
"""Wait for the exact overlap producer, then require and audit its complete outputs."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True)
    a=ap.parse_args()
    config=json.loads(a.config.read_text())
    for path,digest in config['pins'].items():
        assert sha(path)==digest,path
    output=Path(config['output'])
    output.mkdir(parents=True,exist_ok=False)
    def state(status,**extra):
        (output/'state.json').write_text(json.dumps(dict(status=status,config_sha256=sha(a.config),**extra),indent=2)+'\n')
    state('waiting_for_exact_overlap_producer')
    while True:
        try:
            p=psutil.Process(config['producer']['pid'])
            alive=(abs(p.create_time()-config['producer']['created'])<.02 and p.status()!=psutil.STATUS_ZOMBIE)
        except psutil.NoSuchProcess:
            alive=False
        if not alive:
            break
        time.sleep(30)
    try:
        plan=json.loads(Path(config['plan']).read_text())
        receipt=json.loads((Path(plan['output'])/'receipt.json').read_text())
        assert receipt['status']=='complete_full_overlap_inventory_requires_independent_readback'
        assert receipt['plan_sha256']==sha(config['plan'])
        for path,digest in config['pins'].items():
            assert sha(path)==digest,path
        state('running_exhaustive_readback')
        command=['/home/bizon/anaconda3/bin/python','scripts/readback_full_domain_overlaps.py','--plan',config['plan'],'--output',str(Path(plan['output'])/'readback.json')]
        with (output/'readback.log').open('w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        result=json.loads((Path(plan['output'])/'readback.json').read_text())
        assert result['status']=='passed_complete_exhaustive_overlap_readback'
        state('complete_exhaustive_overlap_readback',readback_sha256=sha(Path(plan['output'])/'readback.json'))
    except Exception as error:
        state('failed_requires_review',error=repr(error))
        raise


if __name__=='__main__':
    main()
