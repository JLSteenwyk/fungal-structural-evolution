#!/usr/bin/env python3
"""Gate sequential sensitivity rate/exposure stages on the audited fit controller."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from advance_paired_fit_benchmark import identity, sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--check-config', action='store_true')
    a = p.parse_args()
    c = json.loads(a.config.read_text())
    digest = sha(a.config)

    def verify():
        if sha(a.config) != digest:
            raise ValueError('Controller configuration changed')
        for name, expected in c['pinned_files'].items():
            if sha(Path(name)) != expected:
                raise ValueError('Pinned file changed: ' + name)

    def parent_complete():
        r = json.loads(Path(c['parent_receipt']).read_text())
        if r['status'] != 'complete_fit_audit_and_descriptive_geometry_benchmark_handoff' or r['config_sha256'] != sha(Path(c['parent_config'])) or r['markers'] != c['markers']:
            raise ValueError('Parent did not complete expected audited marker grid')
        for stage in r['stages']:
            if sha(Path(stage['receipt_path'])) != stage['receipt_sha256']:
                raise ValueError('Parent stage receipt changed')
        return r

    def live():
        current = identity(c['parent_pid'])
        if current and current['start_ticks'] == c['parent_start_ticks'] and current['state'] != 'Z':
            if current['command'] != c['parent_command']:
                raise ValueError('Parent command changed')
            return True
        return False

    verify()
    if not live():
        parent_complete()
    if a.check_config:
        print('Pinned files and current parent identity/completion verified; no stages run')
        return
    root = Path(c['output'])
    root.mkdir(parents=True, exist_ok=True)
    lock = (root / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (root / 'receipt.json').exists():
        raise FileExistsError('Controller already completed')
    pin = root / 'config_sha256.txt'
    if pin.exists() and pin.read_text().strip() != digest:
        raise ValueError('Configuration changed across restart')
    pin.write_text(digest + '\n')
    while live():
        print('waiting_for_audited_fit_controller', c['parent_pid'], flush=True)
        time.sleep(30)
    verify()
    parent_complete()
    parent_hash = sha(Path(c['parent_receipt']))
    completed = []
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for stage in c['stages']:
        verify()
        if sha(Path(c['parent_receipt'])) != parent_hash:
            raise ValueError('Parent completion changed')
        for prior in completed:
            if sha(Path(prior['receipt'])) != prior['receipt_sha256']:
                raise ValueError('Earlier stage changed')
        checkpoint = root / (stage['name'] + '.json')
        target = Path(stage['receipt'])
        if checkpoint.exists():
            done = json.loads(checkpoint.read_text())
            if done['config_sha256'] != digest or sha(target) != done['receipt_sha256']:
                raise ValueError('Completed checkpoint changed')
        else:
            if target.parent.exists():
                raise FileExistsError('Uncheckpointed stage requires review: ' + str(target.parent))
            available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:'))) * 1024
            if available < c['resources']['minimum_available_memory_bytes'] or shutil.disk_usage(root).free < c['resources']['minimum_free_disk_bytes']:
                raise RuntimeError('Resource headroom insufficient; stage not started')
            command = [sys.executable] + stage['command']
            print('starting', stage['name'], flush=True)
            with (root / (stage['name'] + '.log')).open('w') as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
            result = json.loads(target.read_text())
            if result['status'] != stage['status']:
                raise ValueError('Stage completion status differs')
            for field, expected in stage.get('expected_counts', {}).items():
                if result[field] != expected:
                    raise ValueError('Stage count mismatch: ' + field)
            done = {'stage': stage['name'], 'config_sha256': digest, 'command': command,
                    'receipt': str(target), 'receipt_sha256': sha(target)}
            checkpoint.write_text(json.dumps(done, indent=2) + '\n')
        completed.append(done)
    result = {'status': 'complete_fcs_rate_and_exposure_frame_handoff', 'config_sha256': digest,
              'parent_receipt_sha256': parent_hash, 'stages': completed, 'markers': c['markers'],
              'interpretation': 'Changed-marker rate comparison, full readback and exposure frame complete. Merge with unchanged baseline markers and revised coupling inference remain separate required steps.'}
    (root / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
