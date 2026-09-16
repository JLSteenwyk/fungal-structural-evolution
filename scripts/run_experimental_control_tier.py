#!/usr/bin/env python3
"""Run a planned control tier, requiring complete predictions before audit/conversion."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args()
    plan = json.loads(a.plan.read_text())
    for path, digest in plan['pins'].items():
        if sha(path) != digest:
            raise ValueError('Changed planned source: ' + path)
    if psutil.virtual_memory().available < plan['minimum_available_memory_gib'] * 2**30:
        raise RuntimeError('Insufficient available host memory')
    if shutil.disk_usage('.').free < plan['minimum_free_disk_gib'] * 2**30:
        raise RuntimeError('Insufficient free disk')
    usage = subprocess.check_output(['nvidia-smi', '-i', plan['gpu_uuid'],
        '--query-gpu=memory.used,utilization.gpu', '--format=csv,noheader,nounits'], text=True)
    memory, utilization = map(int, usage.strip().split(','))
    if memory > 256 or utilization > 5:
        raise RuntimeError('Planned GPU is no longer idle')
    for directory in plan['output_directories']:
        if Path(directory).exists():
            raise FileExistsError('Fresh staged outputs required: ' + directory)
    os.sched_setaffinity(0, plan['cpu_affinity'])
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=plan['gpu_uuid'],
        OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='4', MKL_NUM_THREADS='1',
        PYTHONUNBUFFERED='1')
    receipt_path = Path(plan['execution_receipt'])
    receipt = dict(status='running', started_at=time.time(), plan_sha256=sha(a.plan),
                   pid=os.getpid(), stages=[])
    def save():
        tmp = receipt_path.with_suffix('.partial')
        tmp.write_text(json.dumps(receipt, indent=2) + '\n')
        tmp.replace(receipt_path)
    save()
    try:
        for i, command in enumerate(plan['commands']):
            print('Starting stage', i, command, flush=True)
            subprocess.run(command, env=env, check=True)
            if i == 0:
                chunk = json.loads(Path(plan['chunk_path']).read_text())
                if (chunk['remaining_eligible'] or chunk['oom_deferred'] or chunk['interrupted']
                    or chunk['new_predictions'] + chunk['cached_predictions'] != plan['expected_predictions']):
                    raise RuntimeError('Prediction tier incomplete; preserve outputs for review')
            receipt['stages'].append(dict(index=i, completed_at=time.time()))
            save()
        receipt['status'] = 'complete_prediction_audit_conversion'
    except Exception as exc:
        receipt['status'] = 'failed_requires_review'
        receipt['error'] = repr(exc)
        raise
    finally:
        receipt['updated_at'] = time.time()
        save()


if __name__ == '__main__':
    main()
