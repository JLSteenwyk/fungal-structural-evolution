#!/usr/bin/env python3
"""Run restartable exploratory alternative-copy trees on audited frozen inputs."""
import argparse
import fcntl
import hashlib
import json
import math
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from Bio import Phylo, SeqIO


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--resources', type=Path, required=True)
    a = p.parse_args()
    inputs = json.loads((a.inputs / 'receipt.json').read_text())
    assert inputs['status'] == 'complete_candidate_copy_tree_preparation'
    for name, digest in inputs['artifacts'].items():
        assert sha(a.inputs / name) == digest, name
    resources = json.loads(a.resources.read_text())
    assert resources['workers'] == 2 and resources['threads_per_worker'] == 2
    meminfo = {r.split(':')[0]: int(r.split()[1]) for r in Path('/proc/meminfo').read_text().splitlines()}
    assert meminfo['MemAvailable'] * 1024 >= resources['minimum_available_memory_bytes']
    assert shutil.disk_usage(a.inputs).free >= resources['minimum_free_disk_bytes']
    exe = Path(shutil.which('iqtree3')).resolve()
    version = subprocess.run([str(exe), '--version'], capture_output=True, text=True, check=True).stdout
    a.output.mkdir(parents=True, exist_ok=True)
    lock = (a.output / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    common = {'input_receipt_sha256': sha(a.inputs / 'receipt.json'), 'executable_sha256': sha(exe),
              'version': version, 'script_sha256': sha(Path(__file__)),
              'resource_plan_sha256': sha(a.resources), 'workers': 2,
              'support': '1000 SH-aLRT replicates, not bootstrap support',
              'interpretation': 'Exploratory alternative gene-copy trees; no accepted replacement or species-tree inference.'}
    config = a.output / 'run_config.json'
    if config.exists():
        assert json.loads(config.read_text()) == common
    else:
        config.write_text(json.dumps(common, indent=2) + '\n')
    paths = sorted(a.inputs.glob('*/input.faa'))
    assert len(paths) == inputs['markers_ready']

    def run(path):
        marker = path.parent.name
        records = list(SeqIO.parse(path, 'fasta'))
        names = [r.id for r in records]
        assert len(names) == len(set(names))
        folder = a.output / marker
        folder.mkdir(exist_ok=True)
        prefix = folder / 'tree'
        seed = int(hashlib.sha256(('fcs-copy-' + marker).encode()).hexdigest()[:8], 16) % 2147483646 + 1
        command = [str(exe), '-s', str(path.resolve()), '-st', 'AA', '-m', 'MFP', '-mset', 'LG,WAG,JTT',
                   '-mfreq', 'F', '-mrate', 'G', '--alrt', '1000', '--keep-ident', '-T', '2',
                   '--mem', '4G', '--seed', str(seed), '--prefix', str(prefix.resolve())]
        configuration = {'marker': marker, 'command': command, 'input_sha256': sha(path),
                         'taxa_or_copies': len(names), 'columns': len(records[0].seq),
                         'common_config_sha256': sha(config)}
        cp = folder / 'run_config.json'
        if cp.exists():
            assert json.loads(cp.read_text()) == configuration
        else:
            cp.write_text(json.dumps(configuration, indent=2) + '\n')
        rp = folder / 'receipt.json'
        if rp.exists():
            old = json.loads(rp.read_text())
            if old['status'] == 'inferred':
                assert old['run_config_sha256'] == sha(cp)
                for name, digest in old['artifacts'].items():
                    assert sha(folder / name) == digest
                return old
        started = time.monotonic()
        with (folder / 'stdout.log').open('a') as handle:
            subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
        tree = Phylo.read(prefix.with_suffix('.treefile'), 'newick')
        tips = [tip.name for tip in tree.get_terminals()]
        assert len(tips) == len(set(tips)) and set(tips) == set(names)
        assert all(clade.branch_length is None or (math.isfinite(clade.branch_length) and clade.branch_length >= 0)
                   for clade in tree.find_clades())
        required = ['tree.treefile', 'tree.iqtree', 'tree.log', 'run_config.json', 'stdout.log']
        result = {'status': 'inferred', 'marker': marker, 'run_config_sha256': sha(cp),
                  'elapsed_seconds': time.monotonic() - started, 'tips': len(tips),
                  'artifacts': {name: sha(folder / name) for name in required}}
        rp.write_text(json.dumps(result, indent=2) + '\n')
        print(marker, 'inferred', flush=True)
        return result

    results = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(run, path) for path in paths]
        try:
            for future in as_completed(futures):
                results.append(future.result())
        except Exception:
            for future in futures:
                future.cancel()
            raise
    result = {'status': 'complete_exploratory_copy_tree_fits', 'markers': len(results),
              'common_config_sha256': sha(config), 'results': sorted(results, key=lambda r: r['marker']),
              'interpretation': 'Tip/branch/artifact checks passed; full model/support audit and biological placement review pending.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
