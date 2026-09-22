#!/usr/bin/env python3
"""Run refreshed AlphaFold native coordinate features after full residue readback."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import psutil
from audit_busco_gene_copies import ROOT, sha


def read(path):
    return json.loads(path.read_text())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args()
    plan, plan_sha = read(a.plan), sha(a.plan)

    def verify():
        if sha(a.plan) != plan_sha:
            raise ValueError('Plan changed')
        for name, digest in plan['pins'].items():
            if sha(ROOT / name) != digest:
                raise ValueError('Changed dependency: ' + name)

    verify()
    control = ROOT / plan['control']
    if control.exists() or any((ROOT / stage['output']).exists() for stage in plan['stages']):
        raise FileExistsError('Inspect existing outputs before any recovery')
    control.mkdir(parents=True)
    state = dict(status='waiting_for_exact_mapping_controller', started_unix=time.time(),
                 plan_sha256=plan_sha, completed_stages=[])

    def save():
        state['updated_unix'] = time.time()
        temporary = control / 'state.tmp'
        temporary.write_text(json.dumps(state, indent=2) + '\n')
        temporary.replace(control / 'state.json')

    save()
    try:
        while True:
            try:
                process = psutil.Process(plan['predecessor_pid'])
                live = process.create_time() == plan['predecessor_create_time'] and process.status() != psutil.STATUS_ZOMBIE
            except psutil.NoSuchProcess:
                live = False
            if not live:
                break
            time.sleep(20)
        verify()
        predecessor = read(ROOT / plan['predecessor_control'] / 'receipt.json')
        mapping_receipt = ROOT / plan['mapping'] / 'receipt.json'
        residue_receipt = ROOT / plan['residue_readback'] / 'receipt.json'
        mapping, residue = read(mapping_receipt), read(residue_receipt)
        if (predecessor['status'] != 'passed_full_refreshed_afdb_residue_mapping_readback'
                or predecessor['plan_sha256'] != sha(ROOT / plan['mapping_plan'])
                or predecessor['mapping_receipt_sha256'] != sha(mapping_receipt)
                or residue['status'] != 'passed_full_refreshed_afdb_residue_mapping_readback'
                or residue['mapping_receipt_sha256'] != sha(mapping_receipt)
                or residue['models'] != plan['models'] or mapping['distinct_models'] != plan['models']):
            raise ValueError('Complete mapping and independent residue QC required')
        mapping_sha = sha(mapping_receipt)
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', CUDA_VISIBLE_DEVICES='')
        env['PATH'] = str(Path(plan['foldseek']).parent) + os.pathsep + env.get('PATH', '')
        for stage in plan['stages']:
            verify()
            if sha(mapping_receipt) != mapping_sha:
                raise ValueError('Mapping receipt changed')
            if shutil.disk_usage(ROOT).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
                raise RuntimeError('Insufficient disk headroom')
            command = [sys.executable] + stage['arguments']
            state.update(status='running_stage', stage=stage['name'], command=command)
            save()
            print('Starting', stage['name'], flush=True)
            with (control / (stage['name'] + '.log')).open('w') as log:
                subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            result_path = ROOT / stage['result']
            result = read(result_path)
            if any(result.get(key) != value for key, value in stage['expected'].items()):
                raise ValueError('Stage output status or scope differs: ' + stage['name'])
            if 'mapping_receipt_sha256' in result and result['mapping_receipt_sha256'] != mapping_sha:
                raise ValueError('Stage belongs to a different mapping')
            if stage['name'] == 'native_extraction':
                if result['source_sha256'] != plan['native_source_sha256']:
                    raise ValueError('Native encoder source differs from the previously qualified build')
            state['completed_stages'].append(dict(name=stage['name'], result=stage['result'], sha256=sha(result_path)))
            save()
        verify()
        final = read(ROOT / plan['stages'][-1]['result'])
        result = dict(status='complete_refreshed_afdb_native_coordinate_features', models=plan['models'],
                      mapping_receipt_sha256=mapping_sha, residue_readback_receipt_sha256=sha(residue_receipt),
                      plan_sha256=plan_sha, stages=state['completed_stages'], totals=final['totals'],
                      scope='Full refreshed AlphaFold native feature extraction, export readback and coordinate reconstruction. PAE binding/qualification, biological accuracy and evolutionary analyses remain separate.')
        (control / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
        state['status'] = result['status']
        save()
    except Exception as exc:
        state.update(status='failed_requires_review', error=repr(exc))
        save()
        raise


if __name__ == '__main__':
    main()
