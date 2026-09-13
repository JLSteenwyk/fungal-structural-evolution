#!/usr/bin/env python3
"""Read back all path rows and independently traverse one hashed pair per marker."""
import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from Bio import Phylo


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def checked(path):
    receipt = json.loads((path / 'receipt.json').read_text())
    for name, digest in receipt.get('artifacts', {}).items():
        if sha(path / name) != digest:
            raise ValueError(f'Changed artifact: {path / name}')
    return receipt


def rows(path):
    with path.open() as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def key(row):
    return tuple(row[k] for k in ('marker', 'taxon_a', 'taxon_b'))


def quantile(values, probability):
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lo, hi = math.floor(position), math.ceil(position)
    return ordered[lo] + (position - lo) * (ordered[hi] - ordered[lo])


def close(actual, expected):
    if not math.isclose(float(actual), expected, rel_tol=1e-9, abs_tol=1e-10):
        raise ValueError(f'Numerical mismatch: {actual} != {expected}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('paths', 'points', 'resampling', 'resampling-audit', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    path_receipt = checked(args.paths)
    point_receipt = checked(args.points)
    audit_receipt = checked(args.resampling_audit)
    if (path_receipt['status'] != 'complete_joint_path_sampling_sensitivity'
            or audit_receipt['status'] != 'complete_paired_resampling_audit'):
        raise ValueError('Incomplete source')
    for source in ('points', 'resampling', 'resampling_audit'):
        if path_receipt['source_receipts'][source] != sha(getattr(args, source) / 'receipt.json'):
            raise ValueError('Source lineage mismatch')
    if audit_receipt['resampling_receipt_sha256'] != sha(args.resampling / 'receipt.json'):
        raise ValueError('Resampling audit lineage mismatch')
    config = json.loads((args.resampling / 'config.json').read_text())
    batches = {(r['marker'], int(r['block_length'])): r
               for r in rows(args.resampling_audit / 'batch_summary.tsv')}
    counts, chosen, seen = {}, {}, set()
    for kind, old_name, new_name in (
            ('accepted', 'path_geometry_points.tsv', 'path_geometry_uncertainty.tsv'),
            ('excluded', 'geometry_exclusions_with_paths.tsv', 'geometry_exclusions_with_uncertainty.tsv')):
        original = {key(r): r for r in rows(args.points / old_name)}
        counts[kind] = 0
        for row in rows(args.paths / new_name):
            identity = key(row)
            if identity in seen or identity not in original:
                raise ValueError('Duplicate or unexpected pair')
            seen.add(identity)
            old = original.pop(identity)
            for field, value in old.items():
                target = 'source_point_uncertainty_status' if field == 'uncertainty_status' else field
                if row[target] != value:
                    raise ValueError(f'Changed inherited field: {identity} {field}')
            if row['uncertainty_status'] != 'conditional_AA_and_3Di_AF_paths_other_models_point_only':
                raise ValueError('Incorrect interpretation status')
            for block in config['block_lengths']:
                prefix = f'block{block}_'
                batch = batches[identity[0], block]
                n = int(batch['valid_draws'])
                if (int(row[prefix + 'valid_draws']) != n or
                        int(row[prefix + 'unestimable_draws']) != int(batch['unestimable_draws'])):
                    raise ValueError('Draw count mismatch')
                enough = n >= max(2, math.ceil(.9 * config['replicates']))
                expected = 'conditional_percentile_sensitivity' if enough else 'insufficient_estimable_draws'
                if row[prefix + 'interval_status'] != expected:
                    raise ValueError('Interval availability mismatch')
                if not enough:
                    for field in ('aa_p025', 'aa_median', 'aa_p975', 'aa_sd',
                                  '3di_af_p025', '3di_af_median', '3di_af_p975', '3di_af_sd',
                                  'paired_sampling_covariance'):
                        if row[prefix + field] != '':
                            raise ValueError('Unexpected insufficient-draw statistic')
                    continue
                for label in ('aa', '3di_af'):
                    values = [float(row[prefix + label + '_' + field])
                              for field in ('p025', 'median', 'p975', 'sd')]
                    if not all(math.isfinite(v) and v >= 0 for v in values) or not values[0] <= values[1] <= values[2]:
                        raise ValueError('Invalid path interval')
                covariance = float(row[prefix + 'paired_sampling_covariance'])
                bound = float(row[prefix + 'aa_sd']) * float(row[prefix + '3di_af_sd'])
                if not math.isfinite(covariance) or abs(covariance) > bound + 1e-10:
                    raise ValueError('Invalid covariance')
            digest = hashlib.sha256('\t'.join(identity).encode()).hexdigest()
            marker = identity[0]
            if marker not in chosen or digest < chosen[marker][0]:
                chosen[marker] = (digest, row)
            counts[kind] += 1
        if original or counts[kind] != point_receipt[kind + '_pairs']:
            raise ValueError('Incomplete cohort')
    summary = {(r['marker'], int(r['block_length'])): r
               for r in rows(args.paths / 'marker_block_summary.tsv')}
    if set(summary) != set(batches):
        raise ValueError('Incomplete marker/block summary')
    pair_counts = {}
    for marker, _, _ in seen:
        pair_counts[marker] = pair_counts.get(marker, 0) + 1
    for identity, row in summary.items():
        if int(row['pairs']) != pair_counts[identity[0]]:
            raise ValueError('Summary pair count mismatch')
        for field in ('valid_draws', 'unestimable_draws'):
            if int(row[field]) != int(batches[identity][field]):
                raise ValueError('Summary draw count mismatch')
    traversals, checks = 0, []
    for marker, (digest, row) in sorted(chosen.items()):
        for block in config['block_lengths']:
            folder = args.resampling / marker / f'block-{block}'
            batch = json.loads((folder / 'receipt.json').read_text())
            paths = {'aa': [], '3di_af': []}
            for name, receipt_hash in sorted(batch['receipt_hashes'].items()):
                path = folder / name
                if sha(path) != receipt_hash:
                    raise ValueError('Changed draw receipt')
                draw = json.loads(path.read_text())
                if draw['status'] == 'unestimable_draw':
                    continue
                if draw['status'] != 'completed_draw':
                    raise ValueError('Unexpected draw status')
                for label, stem in (('aa', 'aa'), ('3di_af', '3di')):
                    tree_path = path.parent / (stem + '.treefile')
                    if sha(tree_path) != draw['artifacts'][tree_path.name]:
                        raise ValueError('Changed tree')
                    tree = Phylo.read(tree_path, 'newick')
                    paths[label].append(tree.distance(row['taxon_a'], row['taxon_b']))
                    traversals += 1
            prefix = f'block{block}_'
            n = len(paths['aa'])
            if n != int(row[prefix + 'valid_draws']):
                raise ValueError('Independent draw count mismatch')
            if row[prefix + 'interval_status'] == 'conditional_percentile_sensitivity':
                for label, values in paths.items():
                    for field, probability in (('p025', .025), ('median', .5), ('p975', .975)):
                        close(row[prefix + label + '_' + field], quantile(values, probability))
                    close(row[prefix + label + '_sd'], statistics.stdev(values))
                close(row[prefix + 'paired_sampling_covariance'], statistics.covariance(paths['aa'], paths['3di_af']))
            checks.append({'marker': marker, 'block_length': block, 'taxon_a': row['taxon_a'],
                           'taxon_b': row['taxon_b'], 'selection_sha256': digest, 'valid_draws': n})
        print(marker, 'independent tree paths checked', flush=True)
    result = {'status': 'complete_joint_path_readback', 'pair_counts': counts,
              'markers': len(chosen), 'independent_pair_block_checks': len(checks),
              'independent_tree_traversals': traversals,
              'scope': 'All inherited fields, cohort identities, summary counts, interval bounds and covariance bounds; one SHA256-selected pair per marker independently traversed in every estimable draw and block with Bio.Phylo.distance and Python statistics. Other pair numerical summaries are not independently recomputed.',
              'source_receipts': {k: sha(getattr(args, k) / 'receipt.json')
                                  for k in ('paths', 'points', 'resampling', 'resampling_audit')},
              'script_sha256': sha(Path(__file__)), 'checks': checks}
    args.output.mkdir(parents=True)
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
