#!/usr/bin/env python3
"""Quantify alignment/envelope boundary differences for every candidate domain hit."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt_path = args.manifest / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    if receipt['status'] != 'complete_all_candidate_domain_extraction_manifest':
        raise ValueError('Manifest is incomplete')
    pins = {str(receipt_path): sha(receipt_path)}
    pins.update({str(args.manifest / name): digest for name, digest in receipt['artifacts'].items()})
    for path, digest in pins.items():
        if sha(path) != digest:
            raise ValueError('Changed manifest artifact: ' + path)
    intervals = pd.read_csv(args.manifest / 'intervals.tsv', sep='\t')
    links = pd.read_csv(args.manifest / 'boundary_links.tsv.gz', sep='\t')
    if intervals.interval_id.duplicated().any() or links.duplicated(['model_key', 'hit_id', 'boundary']).any():
        raise ValueError('Duplicate interval or boundary identity')
    if set(links.boundary) != {'alignment', 'envelope'} or set(links.interval_id) != set(intervals.interval_id):
        raise ValueError('Boundary or interval universe differs')
    joined = links.merge(intervals[['interval_id', 'model_key', 'start', 'end', 'residues']],
                         on=['interval_id', 'model_key'], how='left', validate='many_to_one')
    if joined.isna().any().any():
        raise ValueError('Unmatched source interval')
    pairs = joined.pivot(index=['model_key', 'hit_id'], columns='boundary', values=['start', 'end', 'residues'])
    if pairs.isna().any().any() or len(pairs) != receipt['candidate_model_hit_pairs']:
        raise ValueError('Incomplete candidate pair grid')
    result = pairs.index.to_frame(index=False)
    for boundary in ['alignment', 'envelope']:
        for field in ['start', 'end', 'residues']:
            result[boundary + '_' + field] = pairs[field, boundary].to_numpy()
    result['n_terminal_extension'] = result.alignment_start - result.envelope_start
    result['c_terminal_extension'] = result.envelope_end - result.alignment_end
    result['added_residues'] = result.envelope_residues - result.alignment_residues
    if ((result[['n_terminal_extension', 'c_terminal_extension']] < 0).any().any()
            or not (result.added_residues == result.n_terminal_extension + result.c_terminal_extension).all()):
        raise ValueError('Invalid alignment/envelope containment')
    for boundary in ['alignment', 'envelope']:
        if not (result[boundary + '_residues'] == result[boundary + '_end'] - result[boundary + '_start'] + 1).all():
            raise ValueError('Interval length mismatch')
    result['added_fraction_of_alignment'] = result.added_residues / result.alignment_residues
    if (len(intervals) != receipt['unique_intervals'] or len(links) != receipt['boundary_links']
            or result.model_key.nunique() != receipt['source_models']):
        raise ValueError('Manifest scope mismatch')
    summary = {
        'status': 'complete_full_candidate_domain_boundary_sensitivity',
        'candidate_model_hit_pairs': len(result),
        'source_models': int(result.model_key.nunique()),
        'identical_boundary_pairs': int((result.added_residues == 0).sum()),
        'different_boundary_pairs': int((result.added_residues > 0).sum()),
        'pairs_with_at_least_10_added_residues': int((result.added_residues >= 10).sum()),
        'pairs_with_at_least_20_percent_added_residues': int((result.added_fraction_of_alignment >= .2).sum()),
        'added_residues_quantiles': {str(q): float(result.added_residues.quantile(q)) for q in [0, .25, .5, .75, .9, .95, .99, 1]},
        'added_fraction_quantiles': {str(q): float(result.added_fraction_of_alignment.quantile(q)) for q in [0, .5, .9, .95, .99, 1]},
        'source_hashes': pins,
        'script_sha256': sha(__file__),
        'scope': 'Every candidate model/hit in the four-policy union, counted once with both alignment and envelope boundaries. Descriptive boundary sensitivity only; thresholds are reporting bins, not validated exclusions. No confidence, PAE, structural-distance or evolutionary-effect validation.'
    }
    args.output.mkdir(parents=True, exist_ok=False)
    target = args.output / 'boundary_sensitivity.tsv.gz'
    result.to_csv(target, sep='\t', index=False, compression={'method': 'gzip', 'mtime': 0})
    for path, digest in pins.items():
        if sha(path) != digest:
            raise ValueError('Source changed during analysis')
    summary['artifacts'] = {target.name: sha(target)}
    (args.output / 'receipt.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
