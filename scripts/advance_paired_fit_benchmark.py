#!/usr/bin/env python3
"""Gate complete paired-fit audits and geometry benchmarks on a pinned producer."""
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

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(pid):
    try:
        text = Path(f'/proc/{pid}/stat').read_text()
        fields = text[text.rfind(')') + 2:].split()
        command = Path(f'/proc/{pid}/cmdline').read_bytes().decode().split('\0')
    except FileNotFoundError:
        return None
    return {'start_ticks': fields[19], 'state': fields[0], 'command': [s for s in command if s]}


def verify_pins(config, root=ROOT):
    for name, expected in config['pinned_files'].items():
        if sha(root / name) != expected:
            raise ValueError('Pinned file changed: ' + name)


def completed_fit(config, root=ROOT):
    path = root / config['fits'] / 'receipt.json'
    r = json.loads(path.read_text())
    if r['status'] != 'complete_matched_topology_point_estimates' or r['config_sha256'] != sha(root / config['fits'] / 'config.json'):
        raise ValueError('Missing successful fit completion or configuration binding')
    expected = config['expected_markers']
    markers = [row['marker'] for row in r['results']]
    if len(markers) != expected or len(set(markers)) != expected:
        raise ValueError('Incomplete or duplicate marker results')
    return r


def geometry_gate(config, root=ROOT):
    r = json.loads((root / config['geometry_audit'] / 'receipt.json').read_text())
    if r['status'] != 'passed_complete_paired_grid_character_and_sampled_geometry_readback':
        raise ValueError('Geometry audit incomplete')
    for key, folder in [('comparisons', config['geometry']), ('inputs', config['inputs'])]:
        if r['source_receipts'][key] != sha(root / folder / 'receipt.json'):
            raise ValueError('Geometry audit source mismatch')
    return r


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--check-config', action='store_true')
    a = p.parse_args(); config_hash = sha(a.config); c = json.loads(a.config.read_text())
    verify_pins(c); geometry_gate(c)
    current = identity(c['producer_pid'])
    if current and current['start_ticks'] == c['producer_start_ticks'] and current['state'] != 'Z':
        if current['command'] != c['producer_command']:
            raise ValueError('Producer command changed')
    else:
        completed_fit(c)
    if a.check_config:
        print('Configuration, source bindings and live producer/completion checked; no output created')
        return
    out = ROOT / c['controller_output']; out.mkdir(parents=True, exist_ok=True)
    lock = (out / '.lock').open('w'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (out / 'receipt.json').exists():
        raise FileExistsError('Controller has completed')
    saved_config = out / 'config_sha256.txt'
    if saved_config.exists() and saved_config.read_text().strip() != config_hash:
        raise ValueError('Controller configuration changed across invocation')
    saved_config.write_text(config_hash + '\n')
    def verify():
        if sha(a.config) != config_hash:
            raise ValueError('Controller configuration changed')
        verify_pins(c)
    while True:
        current = identity(c['producer_pid'])
        if current is None or current['start_ticks'] != c['producer_start_ticks'] or current['state'] == 'Z':
            break
        if current['command'] != c['producer_command']:
            raise ValueError('Producer command changed')
        print('waiting_for_pinned_fit_producer', c['producer_pid'], flush=True)
        time.sleep(30)
    verify(); completed_fit(c); geometry = geometry_gate(c)
    available = int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))
    if available < c['resources']['minimum_available_memory_gib'] * 1024**2 or shutil.disk_usage(ROOT).free < c['resources']['output_allowance_gb'] * 10**9:
        raise RuntimeError('Insufficient current resource headroom; no stages started')
    stage_specs = [
        ('audit', 'summarize_paired_marker_fits.py', ['--inputs', c['inputs'], '--models', c['models'], '--fits', c['fits'], '--output', c['audit_output']], c['audit_output'], 'complete_paired_fit_audit'),
        ('paths', 'benchmark_paired_site_tree_points.py', ['--inputs', c['inputs'], '--fits', c['fits'], '--audit', c['audit_output'], '--geometry', c['geometry'], '--review', c['review'], '--output', c['path_output']], c['path_output'], 'complete_paired_site_tree_path_point_benchmark'),
        ('ranks', 'summarize_tree_path_geometry_ranks.py', ['--benchmark', c['path_output'], '--output', c['rank_output']], c['rank_output'], 'complete_descriptive_within_marker_rank_summary'),
        ('figure', 'plot_tree_path_geometry_ranks.py', ['--summary', c['rank_output'], '--source-label', c['source_label'], '--output', c['figure_output']], c['figure_output'], 'complete_descriptive_figure')]
    env = os.environ.copy(); env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    stages = []
    for label, script, args, folder, status in stage_specs:
        verify(); checkpoint = out / (label + '.json'); target = ROOT / folder / 'receipt.json'
        if checkpoint.exists():
            done = json.loads(checkpoint.read_text())
            if done['config_sha256'] != config_hash or sha(target) != done['receipt_sha256']:
                raise ValueError('Previously completed stage changed')
        else:
            if (ROOT / folder).exists():
                raise FileExistsError('Uncheckpointed stage output requires review: ' + folder)
            command = [sys.executable, str(ROOT / 'scripts' / script)] + args
            print('starting', label, flush=True)
            with (out / (label + '.log')).open('w') as log:
                subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            result = json.loads(target.read_text())
            if result['status'] != status:
                raise ValueError('Stage did not establish completion')
            done = {'stage': label, 'command': command, 'config_sha256': config_hash, 'receipt_path': str(target.relative_to(ROOT)), 'receipt_sha256': sha(target)}
            checkpoint.write_text(json.dumps(done, indent=2) + '\n')
        stages.append(done)
    verify()
    paths = json.loads((ROOT / c['path_output'] / 'receipt.json').read_text())
    if paths['markers'] != c['expected_markers'] or paths['accepted_pairs'] != geometry['counts']['accepted'] or paths['excluded_pairs'] != geometry['counts']['excluded']:
        raise ValueError('Final benchmark grid differs')
    r = {'status': 'complete_fit_audit_and_descriptive_geometry_benchmark_handoff', 'config_sha256': config_hash,
         'fit_receipt_sha256': sha(ROOT / c['fits'] / 'receipt.json'), 'geometry_audit_sha256': sha(ROOT / c['geometry_audit'] / 'receipt.json'),
         'stages': stages, 'markers': paths['markers'], 'accepted_pairs': paths['accepted_pairs'], 'excluded_pairs': paths['excluded_pairs'],
         'interpretation': 'Point-estimate benchmark and descriptive rank comparisons only. Branch uncertainty, phylogenetic dependence, prediction circularity and biological interpretation remain unresolved. Figure requires visual review.'}
    (out / 'receipt.json').write_text(json.dumps(r, indent=2) + '\n'); print(json.dumps(r, indent=2), flush=True)


if __name__ == '__main__':
    main()
