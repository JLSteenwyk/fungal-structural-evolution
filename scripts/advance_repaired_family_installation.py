#!/usr/bin/env python3
"""Wait for a pinned family rebuild, independently validate, and install its tree."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import psutil


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--check-config', action='store_true')
    a = ap.parse_args()
    config_hash = sha(a.config)
    c = json.loads(a.config.read_text())

    def verify():
        if sha(a.config) != config_hash:
            raise ValueError('Controller configuration changed')
        for name, expected in c['pins'].items():
            if sha(Path(name)) != expected:
                raise ValueError('Changed pinned file: ' + name)

    def live(spec):
        try:
            p = psutil.Process(spec['pid'])
            if p.create_time() != spec['created'] or p.status() == psutil.STATUS_ZOMBIE:
                return False
            if p.cmdline() != spec['command']:
                raise ValueError('Producer command changed')
            return True
        except psutil.NoSuchProcess:
            return False

    verify()
    if a.check_config:
        print('Pinned files verified; producer live:', live(c['producer']))
        return
    out = Path(c['output'])
    out.mkdir(parents=True, exist_ok=False)
    state = dict(status='waiting_for_rebuild', pid=os.getpid(),
                 created=psutil.Process().create_time(), config_sha256=config_hash,
                 started_at=time.time())

    def save():
        p = out / 'state.partial'
        p.write_text(json.dumps(state, indent=2) + '\n')
        p.replace(out / 'state.json')

    save()
    os.sched_setaffinity(0, c['cpu_affinity'])
    try:
        while live(c['producer']):
            time.sleep(30)
        verify()
        plan_path = Path(c['rebuild_plan'])
        plan = json.loads(plan_path.read_text())
        for name, expected in plan['pinned_files'].items():
            if sha(Path(name)) != expected:
                raise ValueError('Rebuild source changed: ' + name)
        rebuild = Path(plan['output'])
        r = json.loads((rebuild / 'receipt.json').read_text())
        cfg = json.loads((rebuild / 'config.json').read_text())
        if (r['status'] != 'complete_isolated_audited_family_tree_rebuild'
                or r['family'] != plan['family']
                or r['config_sha256'] != sha(rebuild / 'config.json')
                or cfg['plan_sha256'] != sha(plan_path)
                or cfg['script_sha256'] != plan['pinned_files']['scripts/rebuild_audited_family_tree.py']):
            raise ValueError('Rebuild completion lineage mismatch')
        for name, expected in r['artifacts'].items():
            if sha(rebuild / name) != expected:
                raise ValueError('Rebuild artifact changed')
        if psutil.virtual_memory().available < c['minimum_available_memory_bytes']:
            raise RuntimeError('Insufficient validation memory headroom')
        if shutil.disk_usage(out).free < c['minimum_free_disk_bytes']:
            raise RuntimeError('Insufficient installation disk headroom')
        state['status'] = 'validating_native_tree'; save()
        readback = out / 'native_readback.json'
        command = [c['python'], 'scripts/validate_repaired_tree_native.py',
                   '--tree', str(rebuild / 'tree.nwk'), '--source', plan['source'],
                   '--output', str(readback)]
        with (out / 'native_readback.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        check = json.loads(readback.read_text())
        if (check['status'] != 'passed_native_repaired_tree_readback'
                or check['tips'] != r['tips']
                or check['tree_sha256'] != r['artifacts']['tree.nwk']
                or check['source_sha256'] != plan['pinned_files'][plan['source']]):
            raise ValueError('Independent readback differs')
        verify()
        production = Path(plan['production_tree'])
        with (out / 'installation.lock').open('w') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            if production.is_symlink() or production.stat().st_size != 0:
                raise ValueError('Production tree is no longer the empty repair target')
            old_hash = sha(production)
            if old_hash != plan['pinned_files'][str(production)]:
                raise ValueError('Production target changed')
            temporary = production.with_name(production.name + '.verified-repair.partial')
            with temporary.open('xb') as dst, (rebuild / 'tree.nwk').open('rb') as src:
                shutil.copyfileobj(src, dst)
                dst.flush(); os.fsync(dst.fileno())
            if sha(temporary) != check['tree_sha256'] or production.stat().st_size != 0:
                raise ValueError('Staged copy or production target changed')
            state['status'] = 'installing_verified_tree'; save()
            os.replace(temporary, production)
            directory = os.open(production.parent, os.O_RDONLY)
            try:
                os.fsync(directory)
            finally:
                os.close(directory)
            if sha(production) != check['tree_sha256']:
                raise ValueError('Installed tree differs')
        state.update(status='complete_verified_family_tree_installation',
                     production_tree=str(production), previous_tree_sha256=old_hash,
                     installed_tree_sha256=check['tree_sha256'], tips=check['tips'],
                     native_readback_sha256=sha(readback),
                     rebuild_receipt_sha256=sha(rebuild / 'receipt.json'),
                     interpretation='Only the audited empty family tree replaced. Native from-trees continuation and biological reconciliation validation remain separate.')
    except Exception as exc:
        state.update(status='failed_requires_review', error=repr(exc))
        raise
    finally:
        state['updated_at'] = time.time(); save()


if __name__ == '__main__':
    main()
