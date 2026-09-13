#!/usr/bin/env python3
"""Infer supported individual marker trees for species-tree discordance analysis."""
import csv
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

ROOT = Path(__file__).resolve().parents[1]
AA = set('ACDEFGHIKLMNPQRSTVWY')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def select_rows(sequences):
    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) != 1:
        raise ValueError('Expected rectangular marker matrix')
    length = next(iter(lengths))
    minimum = max(50, math.ceil(.3 * length))
    included, excluded = {}, {}
    for name, seq in sequences.items():
        observed = sum(aa in AA for aa in seq)
        if observed >= minimum:
            included[name] = seq
        else:
            excluded[name] = observed
    return included, excluded, minimum


def main():
    matrix_dir = ROOT / 'results/phylogeny/profile-matrix-50-v1'
    matrix_path = matrix_dir / 'matrix.faa'
    source_receipt = json.loads((matrix_dir / 'receipt.json').read_text())
    for filename, checksum in source_receipt['artifacts'].items():
        if sha(matrix_dir / filename) != checksum:
            raise ValueError(f'Changed matrix artifact: {filename}')
    with matrix_path.open() as handle:
        matrix = {r.id: str(r.seq).upper() for r in SeqIO.parse(handle, 'fasta')}
    bounds = {}
    with (matrix_dir / 'site_mapping.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            marker, position = row['marker'], int(row['matrix_column_1based']) - 1
            bounds.setdefault(marker, []).append(position)
    for columns in bounds.values():
        if columns != list(range(columns[0], columns[-1] + 1)):
            raise ValueError('Marker columns are not contiguous')
    root = ROOT / 'results/phylogeny/marker-gene-trees-v2'
    root.mkdir(parents=True, exist_ok=True)
    lock = (root / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    executable = shutil.which('iqtree3')
    if executable is None:
        raise FileNotFoundError('iqtree3')
    version = subprocess.run([executable, '--version'], capture_output=True, text=True, check=True).stdout

    def run(marker, columns):
        sequences = {name: seq[columns[0]:columns[-1] + 1] for name, seq in matrix.items()}
        included, excluded, minimum = select_rows(sequences)
        directory = root / marker
        directory.mkdir(exist_ok=True)
        path = directory / 'input.faa'
        text = ''.join(f'>{name}\n{seq}\n' for name, seq in sorted(included.items()))
        if path.exists() and path.read_text() != text:
            raise ValueError('Changed marker gene-tree input')
        path.write_text(text)
        (directory / 'coverage_filter.json').write_text(json.dumps({
            'minimum_unambiguous_residues': minimum, 'excluded_taxa_observed_residues': excluded,
            'note': 'Per-marker coverage rule only; taxa remain in the full project.'}, indent=2) + '\n')
        if len(included) < 4:
            return {'marker': marker, 'status': 'insufficient_taxa', 'taxa': len(included)}
        seed = int(hashlib.sha256(marker.encode()).hexdigest()[:8], 16) % 2147483646 + 1
        prefix = directory / 'tree'
        command = [executable, '-s', str(path), '-st', 'AA', '-m', 'MFP', '-mset', 'LG,WAG,JTT',
                   '-mfreq', 'F', '-mrate', 'G', '--alrt', '1000', '-T', '2', '--mem', '4G',
                   '--seed', str(seed), '--prefix', str(prefix)]
        config = {'marker': marker, 'command': command, 'version': version,
                  'input_sha256': sha(path), 'matrix_receipt_sha256': sha(matrix_dir / 'receipt.json'),
                  'taxa': len(included), 'columns': len(columns), 'support': 'SH-aLRT 1000 replicates; not bootstrap support'}
        config_path = directory / 'run_config.json'
        if config_path.exists() and json.loads(config_path.read_text()) != config:
            raise ValueError('Changed gene-tree configuration')
        config_path.write_text(json.dumps(config, indent=2) + '\n')
        receipt_path = directory / 'receipt.json'
        tree_path = prefix.with_suffix('.treefile')
        if receipt_path.exists():
            old = json.loads(receipt_path.read_text())
            if old['status'] == 'inferred' and tree_path.exists() and old['tree_sha256'] == sha(tree_path):
                return old
        start = time.monotonic()
        with (directory / 'stdout.log').open('a') as log:
            process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        result = dict(config, returncode=process.returncode, elapsed_seconds=time.monotonic() - start, status='failed')
        if process.returncode == 0:
            tree = Phylo.read(tree_path, 'newick')
            tips = [tip.name for tip in tree.get_terminals()]
            if len(set(tips)) != len(tips) or set(tips) != set(included):
                raise ValueError('Gene-tree taxon mismatch')
            result.update(status='inferred', tree_path=str(tree_path.relative_to(ROOT)), tree_sha256=sha(tree_path))
        receipt_path.write_text(json.dumps(result, indent=2) + '\n')
        return result

    results = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(run, marker, columns) for marker, columns in sorted(bounds.items())]
        try:
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                print(result['marker'], result['status'], flush=True)
                if result['status'] == 'failed':
                    raise RuntimeError('Marker job failed; inspect its stdout log before resuming')
        except Exception:
            for future in futures:
                future.cancel()
            raise
    (root / 'receipt.json').write_text(json.dumps({'markers': len(results),
        'results': sorted(results, key=lambda row: row['marker'])}, indent=2) + '\n')


if __name__ == '__main__':
    main()
