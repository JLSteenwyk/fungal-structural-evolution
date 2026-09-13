#!/usr/bin/env python3
"""Append conditional joint-draw path summaries to the exact-site geometry benchmark."""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import read_table, sha
from prepare_paired_phylogenetic_inputs import write_table
from run_paired_marker_fits import tree_edges


def joint_path_statistics(aa, structure, masks):
    aa, structure, masks = map(np.asarray, (aa, structure, masks))
    if (aa.ndim != 2 or structure.shape != aa.shape or masks.ndim != 2
            or masks.shape[1] != aa.shape[1] or len(aa) < 2
            or masks.dtype != bool or not np.isfinite(aa).all()
            or not np.isfinite(structure).all() or (aa < 0).any() or (structure < 0).any()):
        raise ValueError('Invalid paired branch draws or path masks')
    paths = {'aa': masks.astype(float) @ aa.T, '3di_af': masks.astype(float) @ structure.T}
    result = {}
    for label, values in paths.items():
        q = np.quantile(values, [.025, .5, .975], axis=1, method='linear')
        for name, array in zip(['p025', 'median', 'p975'], q):
            result[label + '_' + name] = array
        result[label + '_sd'] = values.std(axis=1, ddof=1)
    x, y = paths['aa'], paths['3di_af']
    result['paired_sampling_covariance'] = ((x - x.mean(axis=1, keepdims=True)) *
                                           (y - y.mean(axis=1, keepdims=True))).sum(axis=1) / (x.shape[1] - 1)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['points', 'inputs', 'fits', 'resampling', 'resampling-audit', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable output')
    points = checked_receipt(a.points); audit = checked_receipt(a.resampling_audit); checked_receipt(a.inputs)
    receipt = json.loads((a.resampling / 'receipt.json').read_text())
    config = json.loads((a.resampling / 'config.json').read_text())
    if (audit['status'] != 'complete_paired_resampling_audit'
            or audit['resampling_receipt_sha256'] != sha(a.resampling / 'receipt.json')
            or receipt['config_sha256'] != sha(a.resampling / 'config.json')
            or any(points['source_receipts'][k] != sha(getattr(a, k) / 'receipt.json')
                   or config['source_receipts'][k] != sha(getattr(a, k) / 'receipt.json') for k in ['inputs', 'fits'])):
        raise ValueError('Mismatched point/resampling lineage')
    batches = {(r['marker'], r['block_length']): r for r in receipt['results']}
    audited = {(r['marker'], int(r['block_length'])): r for r in read_table(a.resampling_audit / 'batch_summary.tsv')}
    grouped = defaultdict(list)
    for accepted, filename in [(True, 'path_geometry_points.tsv'), (False, 'geometry_exclusions_with_paths.tsv')]:
        for row in read_table(a.points / filename):
            grouped[row['marker']].append((accepted, row))
    expected = {(m, b) for m in grouped for b in config['block_lengths']}
    if set(batches) != expected or set(audited) != expected:
        raise ValueError('Incomplete marker/block universe')
    summaries = []; outputs = {True: [], False: []}
    for marker, pairs in sorted(grouped.items()):
        taxa = sorted(x.id for x in SeqIO.parse(a.inputs / marker / 'aa.faa', 'fasta'))
        index = {t: i for i, t in enumerate(taxa)}
        splits = sorted(tree_edges(a.fits / marker / 'aa.treefile', set(taxa)))
        membership = np.array([[t in split for split in splits] for t in taxa], dtype=bool)
        identities = [(r['taxon_a'], r['taxon_b']) for _, r in pairs]
        if len(set(identities)) != len(identities) or len(identities) != len(taxa) * (len(taxa) - 1) // 2 or any(x >= y or x not in index or y not in index for x, y in identities):
            raise ValueError('Invalid pair grid')
        for block in config['block_lengths']:
            folder = a.resampling / marker / f'block-{block}'
            if sha(folder / 'receipt.json') != batches[marker, block]['receipt_sha256']:
                raise ValueError('Batch changed')
            batch = json.loads((folder / 'receipt.json').read_text()); draws = {'aa': [], '3di': []}; invalid = 0
            if set(batch['receipt_hashes']) != {f'replicate-{i:04d}/receipt.json' for i in range(config['replicates'])}:
                raise ValueError('Incomplete draw identities')
            for name, digest in sorted(batch['receipt_hashes'].items()):
                path = folder / name
                if sha(path) != digest:
                    raise ValueError('Joint draw changed')
                draw = json.loads(path.read_text())
                if draw['status'] == 'unestimable_draw':
                    if draw['fits']:
                        raise ValueError('Unexpected fitted unestimable draw')
                    invalid += 1; continue
                if draw['status'] != 'completed_draw' or draw['splits'] != [list(s) for s in splits]:
                    raise ValueError('Draw topology differs')
                for label in draws:
                    draws[label].append(draw['fits'][label]['branch_lengths'])
            n = len(draws['aa']); ar = audited[marker, block]
            if n != batch['completed_draws'] or invalid != batch['unestimable_draws'] or n != int(ar['valid_draws']) or invalid != int(ar['unestimable_draws']):
                raise ValueError('Audited draw accounting differs')
            enough = n >= max(2, math.ceil(.9 * config['replicates']))
            prefix = f'block{block}_'
            for offset in range(0, len(pairs), 512):
                chunk = pairs[offset:offset + 512]
                masks = np.array([membership[index[r['taxon_a']]] ^ membership[index[r['taxon_b']]] for _, r in chunk])
                stats = joint_path_statistics(draws['aa'], draws['3di'], masks) if enough else {}
                for j, (_, row) in enumerate(chunk):
                    row[prefix + 'valid_draws'] = n; row[prefix + 'unestimable_draws'] = invalid
                    row[prefix + 'interval_status'] = 'conditional_percentile_sensitivity' if enough else 'insufficient_estimable_draws'
                    for key in ['aa_p025', 'aa_median', 'aa_p975', 'aa_sd', '3di_af_p025', '3di_af_median', '3di_af_p975', '3di_af_sd', 'paired_sampling_covariance']:
                        row[prefix + key] = float(stats[key][j]) if enough else ''
            summaries.append({'marker': marker, 'block_length': block, 'pairs': len(pairs), 'valid_draws': n, 'unestimable_draws': invalid})
        for accepted, row in pairs:
            row['source_point_uncertainty_status'] = row['uncertainty_status']
            row['uncertainty_status'] = 'conditional_AA_and_3Di_AF_paths_other_models_point_only'
            outputs[accepted].append(row)
        print(marker, 'joint path summaries complete', flush=True)
    if len(outputs[True]) != points['accepted_pairs'] or len(outputs[False]) != points['excluded_pairs']:
        raise ValueError('Changed geometry cohort')
    a.output.mkdir(parents=True)
    for filename, rows in [('path_geometry_uncertainty.tsv', outputs[True]), ('geometry_exclusions_with_uncertainty.tsv', outputs[False]), ('marker_block_summary.tsv', summaries)]:
        write_table(a.output / filename, rows)
    result = {'status': 'complete_joint_path_sampling_sensitivity', 'markers': len(grouped),
              'accepted_pairs': len(outputs[True]), 'excluded_pairs': len(outputs[False]),
              'source_receipts': {k: sha(getattr(a, k) / 'receipt.json') for k in ['points', 'inputs', 'fits', 'resampling', 'resampling_audit']},
              'script_sha256': sha(Path(__file__)),
              'interpretation': 'Same pair cohorts and geometry fields as the point benchmark. AA and 3Di AF path percentiles sum branches within each joint draw; covariance is sampling covariance, not evolutionary coupling. Circular blocks 1/10/30 do not cover all spatial dependencies. Conditional sensitivity only, not calibrated confidence, physical displacement, selection or an acceleration test. Other 3Di models retain point estimates only.',
              'artifacts': {p.name: sha(p) for p in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
