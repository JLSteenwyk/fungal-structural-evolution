#!/usr/bin/env python3
"""Validate paired resampling artifacts and summarize conditional branch intervals."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha
from prepare_paired_phylogenetic_inputs import write_table
from run_paired_branch_resampling import sampled_columns, draw_seed
from run_paired_marker_fits import tree_edges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['inputs', 'fits', 'resampling', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    checked_receipt(args.inputs)
    config_path = args.resampling / 'config.json'
    config = json.loads(config_path.read_text())
    receipt = json.loads((args.resampling / 'receipt.json').read_text())
    if (receipt['status'] != 'complete_paired_resampling_execution'
            or receipt['config_sha256'] != sha(config_path)
            or config['source_receipts']['inputs'] != sha(args.inputs / 'receipt.json')
            or config['source_receipts']['fits'] != sha(args.fits / 'receipt.json')):
        raise ValueError('Resampling provenance mismatch')
    expected_markers = {r['marker'] for r in csv.DictReader((args.inputs / 'marker_summary.tsv').open(), delimiter='\t')
                        if r['status'] == 'ready_for_inference'}
    expected_batches = {(marker, block) for marker in expected_markers for block in config['block_lengths']}
    if len(receipt['results']) != len(expected_batches) or {(r['marker'], r['block_length']) for r in receipt['results']} != expected_batches:
        raise ValueError('Incomplete resampling batches')
    if args.output.exists():
        raise FileExistsError('Use a new immutable output')
    summaries, batches = [], []
    total_draws = unestimable = total_fits = fits_with_warnings = 0
    for batch in receipt['results']:
        marker, block = batch['marker'], batch['block_length']
        folder = args.resampling / marker / f'block-{block}'
        if sha(folder / 'receipt.json') != batch['receipt_sha256']:
            raise ValueError('Changed batch receipt')
        br = json.loads((folder / 'receipt.json').read_text())
        count = config['replicates']
        if set(br['receipt_hashes']) != {f'replicate-{i:04d}/receipt.json' for i in range(count)}:
            raise ValueError('Missing replicate receipts')
        original = {label: {r.id: str(r.seq) for r in SeqIO.parse(args.inputs / marker / filename, 'fasta')}
                    for label, filename in [('aa', 'aa.faa'), ('3di', '3di.faa')]}
        taxa = sorted(original['aa'])
        length = len(original['aa'][taxa[0]])
        topology = args.fits / marker / 'aa.treefile'
        splits = sorted(tree_edges(topology, set(taxa)))
        point = {label: tree_edges(args.fits / marker / filename, set(taxa))
                 for label, filename in [('aa', 'aa.treefile'), ('3di', '3di_af.treefile')]}
        values = {'aa': [], '3di': []}
        invalid = 0
        for replicate in range(count):
            work = folder / f'replicate-{replicate:04d}'
            rp = work / 'receipt.json'
            if sha(rp) != br['receipt_hashes'][f'replicate-{replicate:04d}/receipt.json']:
                raise ValueError('Changed replicate receipt')
            r = json.loads(rp.read_text())
            request_path = work / 'request.json'
            request = json.loads(request_path.read_text())
            if (r['request_sha256'] != sha(request_path) or request['parent_config_sha256'] != sha(config_path)
                    or request['topology_sha256'] != sha(topology)):
                raise ValueError('Changed replicate request')
            for name, checksum in r['artifacts'].items():
                if sha(work / name) != checksum:
                    raise ValueError('Changed resampling artifact')
            seed = draw_seed(marker, block, replicate)
            indices = sampled_columns(length, block, seed)
            expected_columns = '\n'.join(str(int(i) + 1) for i in indices) + '\n'
            if (work / 'columns_1based.txt').read_text() != expected_columns:
                raise ValueError('Resampled columns differ from deterministic design')
            draw_sequences = {}
            for label in original:
                draw_sequences[label] = {taxon: ''.join(original[label][taxon][i] for i in indices) for taxon in taxa}
                actual = {r.id: str(r.seq) for r in SeqIO.parse(work / (label + '.faa'), 'fasta')}
                if actual != draw_sequences[label] or sha(work / (label + '.faa')) != request['alignment_sha256'][label]:
                    raise ValueError('Resampled FASTA differs from original columns')
            reason = 'all_missing_taxon' if any(set(s) == {'?'} for s in draw_sequences['aa'].values()) else ''
            if not reason:
                for label in ['aa', '3di']:
                    if not any(len(set(column) - {'?'}) > 1 for column in zip(*draw_sequences[label].values())):
                        reason = 'no_variable_columns_' + label
                        break
            if r['reason'] != reason:
                raise ValueError('Unestimable draw reason differs')
            total_draws += 1
            if reason:
                if r['status'] != 'unestimable_draw' or r['fits']:
                    raise ValueError('Unexpected invalid-draw state')
                invalid += 1
                unestimable += 1
                continue
            if r['status'] != 'completed_draw' or r['splits'] != [list(s) for s in splits]:
                raise ValueError('Unexpected completed draw state/splits')
            for label in ['aa', '3di']:
                report = (work / (label + '.iqtree')).read_text()
                if 'Model of substitution: ' + config['models'][label] + '\n' not in report:
                    raise ValueError('Replicate model differs')
                fitted = tree_edges(work / (label + '.treefile'), set(taxa))
                if set(fitted) != set(splits):
                    raise ValueError('Replicate topology differs')
                lengths = np.array([fitted[s] for s in splits])
                if not np.allclose(lengths, r['fits'][label]['branch_lengths'], rtol=0, atol=1e-12):
                    raise ValueError('Recorded branch values differ from trees')
                values[label].append(lengths)
                total_fits += 1
                fits_with_warnings += 'WARNING:' in report or 'WARNING:' in (work / (label + '.log')).read_text()
        valid = count - invalid
        if valid != br['completed_draws'] or invalid != br['unestimable_draws']:
            raise ValueError('Batch completion counts differ')
        batches.append({'marker': marker, 'block_length': block, 'draws': count, 'valid_draws': valid,
                        'unestimable_draws': invalid, 'branches': len(splits)})
        enough = valid >= max(2, math.ceil(.9 * count))
        for label in values:
            values[label] = np.asarray(values[label])
        for i, split in enumerate(splits):
            covariance = float(np.cov(values['aa'][:, i], values['3di'][:, i], ddof=1)[0, 1]) if enough else ''
            for label in ['aa', '3di']:
                samples = values[label][:, i] if valid else np.array([])
                quantiles = np.quantile(samples, [.025, .5, .975], method='linear') if enough else ['', '', '']
                summaries.append({'marker': marker, 'block_length': block, 'split_taxa': ','.join(split),
                    'fit': label, 'point_estimate': point[label][split], 'planned_draws': count, 'valid_draws': valid,
                    'interval_status': 'conditional_percentile_sensitivity' if enough else 'insufficient_estimable_draws',
                    'percentile_2_5': quantiles[0], 'median': quantiles[1], 'percentile_97_5': quantiles[2],
                    'fraction_le_1e_minus5': float(np.mean(samples <= 1e-5)) if valid else '',
                    'sampling_sd': float(np.std(samples, ddof=1)) if enough else '',
                    'paired_sampling_covariance': covariance})
        print(marker, block, valid, invalid, flush=True)
    args.output.mkdir(parents=True)
    write_table(args.output / 'conditional_intervals.tsv', summaries)
    write_table(args.output / 'batch_summary.tsv', batches)
    result = {'status': 'complete_paired_resampling_audit', 'markers': len(expected_markers),
        'batches': len(batches), 'attempted_paired_draws': total_draws, 'unestimable_draws': unestimable,
        'validated_fits': total_fits, 'fits_with_warnings': fits_with_warnings, 'branch_interval_rows': len(summaries),
        'resampling_receipt_sha256': sha(args.resampling / 'receipt.json'),
        'script_sha256': sha(Path(__file__)), 'resampling_helper_sha256': sha(Path(__file__).with_name('run_paired_branch_resampling.py')),
        'numpy_version': np.__version__, 'quantile_method': 'linear; 2.5 and 97.5 percentiles, at least 90% draws estimable',
        'interpretation': 'Conditional site/circular-block sampling sensitivity on fixed sequence marker topologies. Intervals are not calibrated for topology/model/prediction error or all spatial feature dependence. Paired sampling covariance describes estimation error, not evolutionary coupling. No branch ratios, p-values or acceleration rankings.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
