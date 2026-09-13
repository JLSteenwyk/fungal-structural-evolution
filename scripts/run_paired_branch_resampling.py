#!/usr/bin/env python3
"""Paired site/block resampling of AA and 3Di branches on fixed marker topologies."""
import argparse
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
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha
from run_paired_marker_fits import tree_edges


def sampled_columns(length, block_length, seed):
    """Circular blocks give every source column equal marginal sampling weight."""
    if length < 1 or not 1 <= block_length <= length:
        raise ValueError('Invalid alignment/block length')
    rng = np.random.Generator(np.random.PCG64(seed))
    starts = rng.integers(0, length, size=math.ceil(length / block_length))
    return ((starts[:, None] + np.arange(block_length)) % length).ravel()[:length]


def draw_seed(marker, block_length, replicate):
    key = f'paired-bootstrap-v1:{marker}:{block_length}:{replicate}'
    return int(hashlib.sha256(key.encode()).hexdigest()[:16], 16)


def atomic_json(path, data):
    temporary = path.with_suffix('.partial')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['inputs', 'models', 'fits', 'audit', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--replicates', type=int, default=200)
    parser.add_argument('--blocks', type=int, nargs='+', default=[1, 10, 30])
    args = parser.parse_args()
    if args.replicates < 2 or len(set(args.blocks)) != len(args.blocks):
        raise ValueError('Invalid replicate count or duplicate block sizes')
    for name in ['inputs', 'models', 'fits', 'audit', 'output']:
        setattr(args, name, getattr(args, name).resolve())
    source = checked_receipt(args.inputs)
    checked_receipt(args.models)
    audit = checked_receipt(args.audit)
    fit_receipt = json.loads((args.fits / 'receipt.json').read_text())
    if audit['fit_receipt_sha256'] != sha(args.fits / 'receipt.json'):
        raise ValueError('Audit and fits differ')
    parent_config = json.loads((args.fits / 'config.json').read_text())
    if (fit_receipt['config_sha256'] != sha(args.fits / 'config.json')
            or parent_config['input_receipt_sha256'] != sha(args.inputs / 'receipt.json')
            or parent_config['model_receipt_sha256'] != sha(args.models / 'receipt.json')):
        raise ValueError('Fit provenance differs')
    executable = shutil.which('iqtree3')
    if not executable or sha(Path(executable)) != parent_config['executable_sha256']:
        raise ValueError('Use the audited IQ-TREE executable')
    markers = [row for row in csv.DictReader((args.inputs / 'marker_summary.tsv').open(), delimiter='\t')
               if row['status'] == 'ready_for_inference']
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    config = {'source_receipts': {name: sha(getattr(args, name) / 'receipt.json')
              for name in ['inputs', 'models', 'fits', 'audit']},
        'script_sha256': sha(Path(__file__)), 'tree_helper_sha256': sha(Path(__file__).with_name('run_paired_marker_fits.py')),
        'numpy_version': np.__version__, 'rng': 'PCG64; seed=first 64 bits of SHA256(paired-bootstrap-v1:marker:block:replicate)',
        'executable': executable, 'executable_sha256': sha(Path(executable)),
        'replicates': args.replicates, 'block_lengths': args.blocks, 'workers': 8, 'threads_per_fit': 1,
        'markers': len(markers), 'planned_paired_draws': len(markers) * args.replicates * len(args.blocks),
        'models': {'aa': 'LG+F+G4', '3di': str(args.models / 'Q.3Di.AF') + '+G4'},
        'method': 'Same resampled columns for every taxon and both alphabets; circular blocks over retained alignment columns, concatenated then truncated to original length. Refit branches and gamma (plus AA empirical frequencies) on fixed original sequence topology. Block=1 is ordinary site resampling.',
        'limits': 'Conditional sampling sensitivity only: no topology, dating, prediction, alignment or model uncertainty. Local blocks do not preserve all nonlocal 3Di feature dependencies; stationary/circular assumptions are not established for proteins. No calibrated acceleration tests or structural/AA ratios.',
        'resource_plan': 'Default: 31,200 paired draws / 62,400 fits. Provisional 0.1-5 seconds/fit yields 0.2-11 hours at eight workers; initialization and I/O can add time. Reserve 20 GB output, eight CPU threads, requested 2 GB/fit. No paid resources.'}
    config_path = args.output / 'config.json'
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        raise ValueError('Changed resampling configuration; use a new output')
    atomic_json(config_path, config)
    loaded = {}
    for row in markers:
        marker = row['marker']
        data = {label: {r.id: str(r.seq) for r in SeqIO.parse(args.inputs / marker / filename, 'fasta')}
                for label, filename in [('aa', 'aa.faa'), ('3di', '3di.faa')]}
        taxa = sorted(data['aa'])
        if set(taxa) != set(data['3di']):
            raise ValueError('Paired taxa differ')
        for taxon in taxa:
            if [x == '?' for x in data['aa'][taxon]] != [x == '?' for x in data['3di'][taxon]]:
                raise ValueError('Paired observation masks differ')
        tree_path = args.fits / marker / 'aa.treefile'
        original_receipt = json.loads((args.fits / marker / 'aa.receipt.json').read_text())
        if sha(tree_path) != original_receipt['artifacts']['aa.treefile']:
            raise ValueError('Changed reference topology')
        edges = sorted(tree_edges(tree_path, set(taxa)))
        arrays = {label: np.array([list(data[label][taxon]) for taxon in taxa]) for label in data}
        length = len(data['aa'][taxa[0]])
        if any(not 1 <= b <= length for b in args.blocks):
            raise ValueError('Block longer than alignment')
        loaded[marker] = taxa, arrays, tree_path, edges, length

    def run(marker, block):
        taxa, arrays, topology, edges, length = loaded[marker]
        folder = args.output / marker / f'block-{block}'
        folder.mkdir(parents=True, exist_ok=True)
        completed = []
        for replicate in range(args.replicates):
            work = folder / f'replicate-{replicate:04d}'
            work.mkdir(exist_ok=True)
            seed = draw_seed(marker, block, replicate)
            indices = sampled_columns(length, block, seed)
            index_text = '\n'.join(str(int(i) + 1) for i in indices) + '\n'
            input_text = {label: ''.join('>' + taxon + '\n' + ''.join(array[i, indices]) + '\n'
                                       for i, taxon in enumerate(taxa)) for label, array in arrays.items()}
            request = {'parent_config_sha256': sha(config_path), 'marker': marker, 'block_length': block,
                'replicate': replicate, 'seed': seed, 'topology_sha256': sha(topology),
                'columns_sha256': hashlib.sha256(index_text.encode()).hexdigest(),
                'alignment_sha256': {label: hashlib.sha256(text.encode()).hexdigest() for label, text in input_text.items()}}
            request_path = work / 'request.json'
            if request_path.exists() and json.loads(request_path.read_text()) != request:
                raise ValueError('Changed resampling request')
            atomic_json(request_path, request)
            receipt_path = work / 'receipt.json'
            if receipt_path.exists():
                old = json.loads(receipt_path.read_text())
                if old['request_sha256'] != sha(request_path):
                    raise ValueError('Changed completed draw request')
                for name, checksum in old['artifacts'].items():
                    if sha(work / name) != checksum:
                        raise ValueError('Changed completed draw artifact')
                completed.append(old['status'])
                continue
            (work / 'columns_1based.txt').write_text(index_text)
            for label, text in input_text.items():
                path = work / (label + '.faa')
                if path.exists() and path.read_text() != text:
                    raise ValueError('Changed partial replicate alignment')
                path.write_text(text)
            counts = (arrays['aa'][:, indices] != '?').sum(axis=1)
            reason = 'all_missing_taxon' if np.any(counts == 0) else ''
            if not reason:
                for label, array in arrays.items():
                    if not any(len(set(array[:, i]) - {'?'}) > 1 for i in indices):
                        reason = 'no_variable_columns_' + label
                        break
            result = {'request_sha256': sha(request_path), 'status': 'unestimable_draw' if reason else 'completed_draw',
                      'reason': reason, 'minimum_observed_per_taxon': int(counts.min()), 'fits': {}}
            if not reason:
                for label in ['aa', '3di']:
                    prefix = work / label
                    command = [executable, '-s', str(work / (label + '.faa')), '-st', 'AA',
                        '-m', config['models'][label], '-te', str(topology), '-keep-ident',
                        '-T', '1', '--mem', '2G', '--seed', str(seed % 2147483646 + 1), '--prefix', str(prefix)]
                    start = time.monotonic()
                    with (work / (label + '.stdout.log')).open('a') as log:
                        process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
                    if process.returncode:
                        raise RuntimeError(f'{marker}/block-{block}/replicate-{replicate}/{label} exit {process.returncode}; inspect before resuming')
                    fitted = tree_edges(work / (label + '.treefile'), set(taxa))
                    if set(fitted) != set(edges):
                        raise ValueError('Resampling fit changed topology')
                    result['fits'][label] = {'command': command, 'elapsed_seconds': time.monotonic() - start,
                                            'branch_lengths': [fitted[edge] for edge in edges]}
            result['splits'] = [list(edge) for edge in edges]
            result['artifacts'] = {p.name: sha(p) for p in work.iterdir()
                                   if p.is_file() and p.name != 'receipt.json' and not p.name.endswith('.partial')}
            atomic_json(receipt_path, result)
            completed.append(result['status'])
        summary = {'marker': marker, 'block_length': block, 'draws': len(completed),
            'completed_draws': completed.count('completed_draw'), 'unestimable_draws': completed.count('unestimable_draw'),
            'receipt_hashes': {f'replicate-{i:04d}/receipt.json': sha(folder / f'replicate-{i:04d}/receipt.json') for i in range(args.replicates)}}
        atomic_json(folder / 'receipt.json', summary)
        print(marker, block, summary['completed_draws'], summary['unestimable_draws'], flush=True)
        return {'marker': marker, 'block_length': block, 'receipt_sha256': sha(folder / 'receipt.json')}

    results = []
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = [pool.submit(run, row['marker'], block) for row in markers for block in args.blocks]
        try:
            for future in as_completed(futures):
                results.append(future.result())
        except Exception:
            for future in futures:
                future.cancel()
            raise
    atomic_json(args.output / 'receipt.json', {'status': 'complete_paired_resampling_execution',
        'config_sha256': sha(config_path), 'results': sorted(results, key=lambda r: (r['marker'], r['block_length'])),
        'interpretation': config['limits']})
    print('All resampling batches complete', flush=True)


if __name__ == '__main__':
    main()
