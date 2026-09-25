#!/usr/bin/env python3
"""Recover full geometry with durable marker checkpoints and the original audit."""
import argparse
import json
import os
import fcntl
import shutil
import psutil
from pathlib import Path
import subprocess
import sys
from inspect_focal_domain_architectures import sha
from assess_pae_sensitivity import checked_receipt


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    a = ap.parse_args(); recovery = json.loads(a.plan.read_text()); plan_hash = sha(a.plan)
    original = Path(recovery['original_plan']); plan = json.loads(original.read_text())
    assert Path('/proc/sys/kernel/random/boot_id').read_text().strip() == recovery['boot_id']
    def verify():
        assert sha(a.plan) == plan_hash
        for path, digest in {**plan['pins'], **recovery['pins']}.items():
            assert sha(path) == digest, path
    verify()
    control = Path(recovery['control']); control.mkdir(parents=True, exist_ok=True)
    lock = (control / '.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', CUDA_VISIBLE_DEVICES='')
    completed = []
    for stage in plan['stages']:
        receipt = Path(stage['output']) / 'receipt.json'
        if not receipt.exists():
            assert psutil.virtual_memory().available >= 96 * 2**30
            assert shutil.disk_usage(control).free >= 100 * 2**30
            assert stage['name'] != 'pae_union', 'Previously completed PAE union missing'
            arguments = list(stage['arguments'])
            if stage['name'] == 'geometry':
                arguments[0] = 'scripts/prepare_checkpointed_paired_site_geometry.py'
                arguments += ['--checkpoints', recovery['checkpoints']]
            command = [sys.executable, *arguments]
            (control / 'state.json').write_text(json.dumps(dict(stage=stage['name'], command=command))+'\n')
            with (control / (stage['name']+'.log')).open('a') as log:
                subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        result = checked_receipt(Path(stage['output']))
        assert all(result.get(k) == v for k, v in stage['expected'].items())
        if stage['name'] == 'geometry':
            assert result['accepted_taxon_pairs'] + result['excluded_taxon_pairs'] == plan['pairs']
        completed.append(dict(name=stage['name'], receipt=str(receipt), sha256=sha(receipt)))
        verify()
    result = dict(status='complete_all_cohort_paired_geometry_and_grid_audit',
                  plan_sha256=sha(original), recovery_plan_sha256=plan_hash,
                  markers=plan['markers'], pairs=plan['pairs'], stages=completed,
                  scope='Checkpointed geometry; full grid/character audit and one numerical pair per marker. Not full numerical validation or evolutionary inference.')
    for target in [control / 'receipt.json', Path(plan['output']) / 'receipt.json']:
        with target.open('x') as f: json.dump(result, f, indent=2); f.write('\n')


if __name__ == '__main__':
    main()
