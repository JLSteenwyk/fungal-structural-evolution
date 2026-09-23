#!/usr/bin/env python3
"""Gate full-domain CPU clustering on the audited interval database."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time
import psutil
from catalog_whole_proteome_structures import sha


def check_partition(path, expected):
    assigned={};counts=Counter();self_rows=set()
    with Path(path).open() as f:
        for line in f:
            representative,member=line.rstrip('\n').split('\t')
            if representative not in expected or member not in expected or member in assigned:
                raise ValueError('Unknown identifier or repeated cluster member')
            assigned[member]=representative;counts[representative]+=1
            if representative==member:self_rows.add(member)
    if set(assigned)!=expected or self_rows!=set(counts):raise ValueError('Incomplete partition or missing representative self-membership')
    return counts


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',required=True,type=Path)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan_hash=sha(args.plan)
    def verify():
        if sha(args.plan)!=plan_hash:raise ValueError('Plan changed')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Pinned input changed: '+p)
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    def state(status,**kw):
        (out/'state.json').write_text(json.dumps(dict(status=status,**kw),indent=2)+'\n')
    state('waiting_for_verified_database')
    while True:
        try:
            process=psutil.Process(plan['predecessor']['pid'])
            live=process.create_time()==plan['predecessor']['create_time'] and process.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify();db=Path(plan['database']);receipt_path=db/'receipt.json'
    receipt=json.loads(receipt_path.read_text())
    if receipt['status']!='complete_domain_foldseek_database_with_full_sequence_and_coordinate_readback' or receipt['plan_sha256']!=sha(plan['database_plan']):raise ValueError('Database verification incomplete')
    for name,h in receipt['artifacts'].items():
        if sha(db/name)!=h:raise ValueError('Database artifact changed')
    expected=set()
    with (db/'domains.lookup').open() as f:
        for line in f:
            _,name,_=line.rstrip('\n').split('\t')
            if name in expected:raise ValueError('Repeated lookup name')
            expected.add(name)
    if len(expected)!=receipt['models']:raise ValueError('Database scope mismatch')
    resource=plan['resources']
    if shutil.disk_usage(out).free<resource['minimum_free_disk_gib']*2**30 or psutil.virtual_memory().available<resource['minimum_available_memory_gib']*2**30:raise ValueError('Insufficient resources')
    prefix=str(db/'domains');cluster=str(out/'clusters');commands=[
        [plan['foldseek'],'cluster',prefix,cluster,str(out/'tmp'),*plan['cluster_arguments']],
        [plan['foldseek'],'createtsv',prefix,prefix,cluster,str(out/'cluster_members.tsv'),'--threads',str(resource['cpu'])]]
    for i,command in enumerate(commands):
        state('running',stage=i,command=command)
        with (out/f'stage_{i}.log').open('w') as log:
            child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,CUDA_VISIBLE_DEVICES=''),start_new_session=True)
            try:
                while child.poll() is None:
                    if shutil.disk_usage(out).free<resource['emergency_free_disk_gib']*2**30:raise RuntimeError('Disk reserve reached')
                    time.sleep(10)
                if child.returncode:raise RuntimeError(f'Native stage {i} exited {child.returncode}')
            except BaseException:
                if child.poll() is None:
                    os.killpg(child.pid,signal.SIGTERM)
                    try:child.wait(timeout=30)
                    except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
                raise
    counts=check_partition(out/'cluster_members.tsv',expected)
    with (out/'cluster_sizes.tsv').open('w') as f:
        f.write('representative\tmodels\n')
        for representative,count in sorted(counts.items()):f.write(f'{representative}\t{count}\n')
    verify()
    # Recheck input database artifacts after native execution.
    for name,h in receipt['artifacts'].items():
        if sha(db/name)!=h:raise ValueError('Database changed during clustering')
    result={'status':'complete_domain_candidate_partition_membership_readback','intervals':len(expected),'excluded_rejected_intervals':receipt['excluded_rejected_intervals'],'exported_intervals_with_missing_backbone':receipt['exported_intervals_with_missing_backbone'],'clusters':len(counts),'singleton_clusters':sum(n==1 for n in counts.values()),'largest_cluster_models':max(counts.values()),'plan_sha256':plan_hash,'database_receipt_sha256':sha(receipt_path),'commands':commands,'artifacts':{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='state.json'},'scope':'Full exported-interval partition and representative membership verified. Both alignment and envelope intervals retained, so alternative boundaries are not independent observations. Candidate domain clustering with confidence-masked seeding, finite prefilter cap and native reassignment; not independent alignment-threshold validation, confidence qualification, orthology, remote homology confirmation or evolutionary events. Boundary and parameter sensitivity remain required; rejected intervals are excluded explicitly.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');state('complete',models=len(expected),clusters=len(counts))


if __name__=='__main__':main()
