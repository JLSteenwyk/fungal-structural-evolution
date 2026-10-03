#!/usr/bin/env python3
"""Run paired 200-draw site/block uncertainty controls over every ready input."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
import fcntl
import json
from pathlib import Path
import shutil
import psutil

from ancestral_chain_attempt import sha
from matched_predictor_resampling import load,perform,add_splits,MODES,SCHEMA,runtime_caps
from matched_predictor_branch_inputs import verify

STATUS='complete_full_matched_predictor_paired_resampling_pending_independent_readback'


def run(path):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path)
    limits=runtime_caps(plan)
    assert psutil.virtual_memory().available>=plan['resources']['memory_gib']*2**30
    source=add_splits(source,Path(plan['point_output']));root=Path(plan['output'])
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(exist_ok=False);lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    (root/'work').mkdir();(root/'cases').mkdir()
    axes=dict(schema=SCHEMA,modes=MODES,replicates_per_mode=plan['replicates_per_mode'],
        original_inputs=sorted(source['configs']),resampling_groups=source['resampling_groups'],
        native_roles=[r[0] for r in __import__('matched_predictor_branch_fits').ROLES],splits=source['splits'])
    (root/'axes.json').write_text(json.dumps(axes,indent=2)+'\n')
    jobs=iter((key,mode,rep) for key in sorted(source['configs']) for mode in MODES for rep in range(plan['replicates_per_mode']))
    manifest=[];counts=Counter();complete=slots=finite=archived=0;total=len(source['configs'])*2*plan['replicates_per_mode']
    with ThreadPoolExecutor(max_workers=plan['resources']['workers']) as pool:
        pending={}
        def submit():
            job=next(jobs,None)
            if job is not None:
                f=pool.submit(perform,plan,source,*job,root);pending[f]=job
        for _ in range(plan['resources']['workers']):submit()
        while pending:
            finished,_=wait(pending,return_when=FIRST_COMPLETED)
            for f in finished:
                key,mode,rep=pending.pop(f);value=f.result();counts.update(value['native_status_counts'])
                complete+=value['complete'];slots+=value['serialized_branch_value_slots'];finite+=value['finite_branch_values']
                archived+=len(value['archive_members'])
                folder=root/'cases'/key/mode;stem=f'{rep:04d}'
                manifest.append(dict(original_input_id=key,mode=mode,replicate=rep,case_id=value['case_id'],
                    receipt=str((folder/(stem+'.json')).relative_to(root)),receipt_sha256=sha(folder/(stem+'.json')),
                    archive=str((folder/(stem+'.tar.gz')).relative_to(root)),archive_sha256=value['archive_sha256'],
                    array=str((folder/(stem+'.npz')).relative_to(root)),array_sha256=value['array_sha256']))
                assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
                print('matched_predictor_resample',len(manifest),'/',total,key,mode,rep,'complete',value['complete'],flush=True)
                submit()
    assert len(manifest)==total==53200 and sum(counts.values())==372400
    assert not list((root/'work').iterdir());manifest.sort(key=lambda r:(r['original_input_id'],r['mode'],r['replicate']))
    (root/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    summary=dict(original_comparison_cases=8750,unique_inputs=133,replicates_per_mode=200,modes=MODES,
        resampling_cases=len(manifest),native_roles=sum(counts.values()),native_status_counts=dict(counts),
        complete_cases=complete,unresolved_cases=len(manifest)-complete,serialized_branch_value_slots=slots,
        finite_branch_values=finite,archived_native_files=archived)
    verify(bindings)
    artifacts={name:sha(root/name) for name in ['axes.json','case_manifest.json']}
    for row in manifest:
        for name in ['receipt','archive','array']:artifacts[row[name]]=row[name+'_sha256']
    receipt=dict(status=STATUS,plan_sha256=sha(path),**summary,artifacts=artifacts,source_hashes=bindings,
                 actual_cgroup_limits=limits,scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as f:f.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
