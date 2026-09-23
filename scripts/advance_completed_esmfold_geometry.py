#!/usr/bin/env python3
"""Build the full completed-cohort PAE union and paired-site geometry benchmark."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import psutil
from readback_domain_boundary_clusters import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Plan changed')
        for path,digest in plan['pins'].items():
            if sha(path)!=digest:raise ValueError('Pinned source changed: '+path)
    verify();out=Path(plan['output'])
    if out.exists() or any(Path(s['output']).exists() for s in plan['stages']):raise FileExistsError('Use fresh outputs')
    audit=json.loads(Path(plan['input_readback']).read_text())
    if (audit['status']!='passed_complete_paired_inputs_from_qualified_arrays_readback'
            or audit['source_receipt_sha256']!=sha(Path(plan['inputs'])/'receipt.json')
            or audit['markers']!=plan['markers']):raise ValueError('Incomplete input qualification')
    if psutil.virtual_memory().available<plan['resources']['minimum_available_memory_gib']*2**30:raise ValueError('Insufficient memory headroom')
    out.mkdir(parents=True);completed=[];started=time.time()
    def state(status,**extra):
        target=out/'state.tmp';target.write_text(json.dumps(dict(status=status,completed_stages=completed,elapsed_seconds=time.time()-started,**extra),indent=2)+'\n');target.replace(out/'state.json')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    try:
        for stage in plan['stages']:
            verify()
            if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk headroom')
            command=[sys.executable,*stage['arguments']];state('running_stage',stage=stage['name'],command=command)
            with (out/(stage['name']+'.log')).open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
            target=Path(stage['output'])/'receipt.json';result=json.loads(target.read_text())
            if any(result.get(k)!=v for k,v in stage['expected'].items()):raise ValueError('Stage completion scope differs')
            if stage['name']=='geometry' and result['accepted_taxon_pairs']+result['excluded_taxon_pairs']!=plan['pairs']:raise ValueError('Incomplete pair universe')
            completed.append(dict(name=stage['name'],receipt=str(target),sha256=sha(target)))
        verify();result=dict(status='complete_all_cohort_paired_geometry_and_grid_audit',plan_sha256=ph,markers=plan['markers'],pairs=plan['pairs'],stages=completed,scope='Complete paired-site geometry and independent full pair-grid/character checks, with one independent numerical geometry check per accepted marker. Not a full numerical audit, tree-path benchmark, branch-rate inference, or uncertainty analysis.')
        (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');state(result['status'])
    except Exception as exc:
        state('failed_requires_review',error=repr(exc));raise


if __name__=='__main__':main()
