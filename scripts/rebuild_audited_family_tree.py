#!/usr/bin/env python3
"""Rebuild an audited empty family tree in isolation with the original FastTree."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time
from Bio import Phylo


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    a = ap.parse_args(); plan = json.loads(a.plan.read_text())
    out = Path(plan['output'])
    if out.exists():
        raise FileExistsError(out)
    for name, digest in plan['pinned_files'].items():
        if sha(Path(name)) != digest:
            raise ValueError('Changed input: '+name)
    disposition = json.loads(Path(plan['disposition']).read_text())
    matches = [r for r in disposition['required_tree_repairs'] if r['family'] == plan['family']]
    if len(matches) != 1 or matches[0]['alignment_status'] != 'content_checks_passed' or matches[0]['tree_status'] != 'empty':
        raise ValueError('Audited empty-tree repair not established')
    row = matches[0]; source = Path(plan['source']); alignment = Path(plan['alignment'])
    if sha(source) != row['source_sha256'] or sha(alignment) != row['alignment_sha256']:
        raise ValueError('Audited sequences changed')
    production = Path(plan['production_tree'])
    if production.stat().st_size != 0:
        raise ValueError('Production tree is no longer empty')
    ids = set(); count = 0
    with source.open() as handle:
        for line in handle:
            if line.startswith('>'):
                ids.add(line[1:].strip()); count += 1
    if len(ids) != count or count != int(row['source_sequences']):
        raise ValueError('Source identity mismatch')
    available = int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:'))) * 1024
    if available < plan['minimum_available_memory_bytes'] or shutil.disk_usage(out.parent).free < plan['minimum_free_disk_bytes']:
        raise RuntimeError('Resource headroom insufficient')
    os.sched_setaffinity(0, plan['cpu_affinity'])
    limit = plan['address_space_limit_bytes']
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    out.mkdir(); command = [str(Path(plan['executable']).resolve()), str(alignment.resolve())]
    config = {'plan_sha256': sha(a.plan), 'command': command, 'script_sha256': sha(Path(__file__)), 'pid': os.getpid(), 'started_at': time.time()}
    (out/'config.json').write_text(json.dumps(config, indent=2)+'\n')
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with (out/'tree.nwk').open('w') as tree, (out/'fasttree.log').open('w') as log:
        result = subprocess.run(command, stdout=tree, stderr=log, env=env)
    if result.returncode:
        (out/'failure.json').write_text(json.dumps({'returncode': result.returncode, 'time': time.time(), 'interpretation': 'Isolated repair failed; production tree untouched.'}, indent=2)+'\n')
        raise RuntimeError(f'FastTree exited {result.returncode}')
    parsed = Phylo.read(out/'tree.nwk', 'newick'); tips = []; stack = [parsed.root]
    while stack:
        node = stack.pop()
        if not node.clades:
            tips.append(node.name)
        stack.extend(node.clades)
        if node.branch_length is not None and (not math.isfinite(node.branch_length) or node.branch_length < 0):
            raise ValueError('Invalid branch length')
    if len(tips) != len(set(tips)) or set(tips) != ids:
        raise ValueError('Tree tip universe differs')
    for name, digest in plan['pinned_files'].items():
        if sha(Path(name)) != digest:
            raise ValueError('Input changed during repair: '+name)
    if sha(a.plan) != config['plan_sha256']:
        raise ValueError('Plan changed during repair')
    receipt = {'status': 'complete_isolated_audited_family_tree_rebuild', 'family': plan['family'], 'tips': len(tips), 'alignment_columns': int(row['alignment_columns']), 'elapsed_seconds': time.time()-config['started_at'], 'peak_child_rss_kib': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss, 'config_sha256': sha(out/'config.json'), 'artifacts': {n: sha(out/n) for n in ['tree.nwk', 'fasttree.log']}, 'interpretation': 'Original FastTree and unchanged audited alignment. Tip and branch checks passed; production not replaced. Installation and downstream orthology continuation remain separate.'}
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__':
    main()
