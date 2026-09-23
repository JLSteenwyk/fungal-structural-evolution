#!/usr/bin/env python3
"""Export and audit full-cohort site rates after the identified paired-fit producer completes.

Derived from advance_site_rate_comparison.py; uses the same gated stage runner."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import psutil


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(spec):
    path = Path(spec['receipt']); r = json.loads(path.read_text())
    if r['status'] != spec['status'] or any(r[k] != v for k,v in spec.get('expected_counts',{}).items()):
        raise ValueError('Incomplete stage: '+str(path))
    for name,digest in r.get('artifacts',{}).items():
        if sha(path.parent/name) != digest:
            raise ValueError('Changed stage artifact: '+name)
    return sha(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--check-config',action='store_true')
    a=ap.parse_args(); c=json.loads(a.config.read_text()); digest=sha(a.config)
    def verify():
        if sha(a.config)!=digest:
            raise ValueError('Configuration changed')
        for name,pin in c['pinned_files'].items():
            if sha(Path(name))!=pin:
                raise ValueError('Pinned input changed: '+name)
    def live():
        try:
            p=psutil.Process(c['producer']['pid'])
            if p.create_time()!=c['producer']['created'] or p.status()==psutil.STATUS_ZOMBIE:
                return False
            if p.cmdline()!=c['producer']['command']:
                raise ValueError('Producer command changed')
            return True
        except psutil.NoSuchProcess:
            return False
    verify()
    if not live():
        for spec in c['prerequisites']:
            validate(spec)
    if a.check_config:
        print('Pins and producer identity/completion verified; no stages run');return
    root=Path(c['output']);root.mkdir(parents=True,exist_ok=True)
    lock=(root/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (root/'receipt.json').exists():
        raise FileExistsError('Controller already complete')
    pin=root/'config_sha256.txt'
    if pin.exists() and pin.read_text().strip()!=digest:
        raise ValueError('Configuration changed across restart')
    pin.write_text(digest+'\n');os.sched_setaffinity(0,c['cpu_affinity'])
    while live():
        print('Waiting for rate producer',c['producer']['pid'],flush=True);time.sleep(30)
    verify();prerequisites={s['receipt']:validate(s) for s in c['prerequisites']};completed=[]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for stage in c['stages']:
        verify()
        for path,pin in prerequisites.items():
            if sha(Path(path))!=pin:raise ValueError('Prerequisite changed')
        for prior in completed:
            if sha(Path(prior['receipt']))!=prior['receipt_sha256']:raise ValueError('Prior stage changed')
        checkpoint=root/(stage['name']+'.json');target=Path(stage['receipt'])
        if checkpoint.exists():
            done=json.loads(checkpoint.read_text())
            if done['config_sha256']!=digest or validate(stage)!=done['receipt_sha256']:
                raise ValueError('Checkpoint differs')
        else:
            if target.parent.exists():raise FileExistsError('Review uncheckpointed stage: '+str(target.parent))
            available=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
            if available<c['resources']['minimum_available_memory_bytes'] or shutil.disk_usage(root).free<c['resources']['minimum_free_disk_bytes']:
                raise RuntimeError('Insufficient resource headroom')
            command=[sys.executable]+stage['command'];print('Starting',stage['name'],flush=True)
            with (root/(stage['name']+'.log')).open('w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
            done={'stage':stage['name'],'config_sha256':digest,'command':command,'receipt':str(target),'receipt_sha256':validate(stage)}
            checkpoint.write_text(json.dumps(done,indent=2)+'\n')
        completed.append(done)
    verify()
    result={'status':'complete_full_cohort_site_rate_exports_and_audits','config_sha256':digest,'prerequisite_receipts':prerequisites,'stages':completed,'interpretation':'Gamma4 and FreeRate4 exports and report audits complete. Optimization sensitivity, model comparison, exposure integration and coupling inference remain separate stages.'}
    (root/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
