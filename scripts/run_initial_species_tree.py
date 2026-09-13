#!/usr/bin/env python3
"""Run a reproducible full-matrix IQ-TREE guide analysis, not a final species tree."""
import argparse
import fcntl
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path
from Bio import Phylo, SeqIO

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    matrix = (args.matrix / 'matrix.faa').resolve()
    matrix_receipt = json.loads((args.matrix / 'receipt.json').read_text())
    if sha(matrix) != matrix_receipt['artifacts']['matrix.faa']:
        raise ValueError('Changed matrix')
    executable = shutil.which('iqtree3')
    if executable is None:
        raise FileNotFoundError('iqtree3')
    version = subprocess.run([executable, '--version'], text=True, capture_output=True, check=True).stdout
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    prefix = (args.output / 'guide').resolve()
    command = [executable, '-s', str(matrix), '-st', 'AA', '-m', 'LG+F+G4',
               '-T', '16', '--mem', '32G', '--seed', '20260913', '--prefix', str(prefix)]
    config = {'command': command, 'version': version, 'input_sha256': sha(matrix),
              'purpose': 'Initial unpartitioned homogeneous-model guide; support, mixture/partition models, discordance and sampling sensitivity remain required.'}
    config_path = args.output / 'run_config.json'
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        raise ValueError('Changed run configuration; choose a new output directory')
    config_path.write_text(json.dumps(config, indent=2) + '\n')
    started = time.monotonic()
    with (args.output / 'stdout.log').open('a') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    receipt = dict(config, returncode=result.returncode, elapsed_seconds=time.monotonic() - started)
    if result.returncode == 0:
        tree_path = prefix.with_suffix('.treefile')
        tree = Phylo.read(tree_path, 'newick')
        with matrix.open() as handle:
            expected = {record.id for record in SeqIO.parse(handle, 'fasta')}
        tips = [tip.name for tip in tree.get_terminals()]
        if len(tips) != len(set(tips)) or set(tips) != expected:
            raise ValueError('Tree taxon set differs from complete matrix')
        if any(clade.branch_length is not None and clade.branch_length < 0 for clade in tree.find_clades()):
            raise ValueError('Negative tree branch length')
        receipt.update(tree_sha256=sha(tree_path), taxa=len(tips), status='guide_inferred_not_final')
    else:
        receipt['status'] = 'failed'
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    if result.returncode:
        raise RuntimeError(f'IQ-TREE exited {result.returncode}')
    print(receipt['status'], receipt['taxa'], flush=True)


if __name__ == '__main__':
    main()
