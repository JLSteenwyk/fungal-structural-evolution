#!/usr/bin/env python3
"""Measure shared-coordinate feature overlap in the paired phylogenetic inputs."""
import argparse
import csv
import gzip
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha
from prepare_paired_phylogenetic_inputs import write_table


def shared_feature_edges(features):
    by_residue = defaultdict(list)
    for column, residues in features.items():
        for residue in set(residues):
            by_residue[residue].append(column)
    return {tuple(sorted(pair)) for columns in by_residue.values() for pair in combinations(columns, 2)}


def component_sizes(nodes, edges):
    parent = {n: n for n in nodes}
    def root(n):
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n
    for a, b in edges:
        a, b = root(a), root(b)
        parent[a] = b
    counts = defaultdict(int)
    for n in nodes:
        counts[root(n)] += 1
    return sorted(counts.values(), reverse=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['inputs', 'encodings', 'snapshot', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    receipts = {name: checked_receipt(getattr(args, name)) for name in ['inputs', 'encodings', 'snapshot']}
    for name in ['encodings', 'snapshot']:
        if receipts['inputs']['source_receipts'][name]['sha256'] != sha(getattr(args, name) / 'receipt.json'):
            raise ValueError('Paired inputs use a different source snapshot')
    if args.output.exists():
        raise FileExistsError('Use a new immutable output')
    encoded = {}
    for row in csv.DictReader((args.encodings / 'model_summary.tsv').open(), delimiter='\t'):
        path = ROOT / row['encoding_path']
        if sha(path) != row['encoding_sha256']:
            raise ValueError('Changed encoding')
        with np.load(path, allow_pickle=False) as source:
            encoded[row['model_name']] = {k: source[k].copy() for k in source.files}
    links = {(r['marker'], r['taxon_id']): Path(r['model_path']).stem
             for r in csv.DictReader((args.snapshot / 'marker_structure_links.tsv').open(), delimiter='\t')}
    mapped = defaultdict(dict)
    with gzip.open(args.snapshot / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            mapped[row['marker'], row['taxon_id']][int(row['matrix_column_1based'])] = int(row['protein_residue_1based'])
    results = []
    for folder in sorted(args.inputs.iterdir()):
        if not folder.is_dir():
            continue
        marker = folder.name
        sequences = {r.id: str(r.seq) for r in SeqIO.parse(folder / '3di.faa', 'fasta')}
        columns = {int(r['paired_column_1based']): int(r['matrix_column_1based'])
                   for r in csv.DictReader((folder / 'columns.tsv').open(), delimiter='\t')}
        for taxon, states in sequences.items():
            key = marker, taxon
            data = encoded[links[key]]
            full_states = str(data['states'])
            features = {}
            partner_separations = []
            for column, state in enumerate(states, 1):
                if state == '?':
                    continue
                position = mapped[key][columns[column]]
                i = position - 1
                if (not data['valid'][i] or data['feature_min_plddt'][i] < 70
                        or data['feature_max_pae'][i] > 10 or full_states[i] != state):
                    raise ValueError('State or feature confidence differs from paired input')
                partner = int(data['partner_residue_1based'][i])
                residues = {position - 1, position, position + 1, partner - 1, partner, partner + 1}
                if min(residues) < 1 or max(residues) > len(full_states):
                    raise ValueError('Feature context outside protein')
                features[column] = residues
                partner_separations.append(abs(partner - position))
            edges = shared_feature_edges(features)
            sizes = component_sizes(features, edges)
            circular_separations = [min(b - a, len(states) - (b - a)) for a, b in edges]
            row = {'marker': marker, 'taxon_id': taxon, 'model_name': links[key],
                'alignment_columns': len(states), 'observed_features': len(features),
                'shared_coordinate_feature_pairs': len(edges), 'overlap_components': len(sizes),
                'largest_overlap_component': max(sizes),
                'median_partner_sequence_separation': float(np.median(partner_separations)),
                'max_partner_sequence_separation': max(partner_separations)}
            for block in [10, 30]:
                count = sum(separation >= block for separation in circular_separations)
                row[f'overlap_pairs_beyond_circular_block_{block}'] = count
                row[f'fraction_overlap_pairs_beyond_circular_block_{block}'] = count / len(edges) if edges else ''
            results.append(row)
    args.output.mkdir(parents=True)
    write_table(args.output / 'feature_overlap.tsv', results)
    total = sum(r['shared_coordinate_feature_pairs'] for r in results)
    receipt = {'status': 'complete_paired_feature_overlap_audit', 'taxon_marker_alignments': len(results),
        'markers': len({r['marker'] for r in results}),
        'observed_features': sum(r['observed_features'] for r in results),
        'shared_coordinate_feature_pairs': total,
        'overlap_pairs_beyond_circular_blocks': {str(b): sum(r[f'overlap_pairs_beyond_circular_block_{b}'] for r in results) for b in [10, 30]},
        'source_receipts': {name: {'path': str(getattr(args, name)), 'sha256': sha(getattr(args, name) / 'receipt.json')} for name in receipts},
        'script_sha256': sha(Path(__file__)),
        'interpretation': 'A feature-overlap edge means two observed 3Di features use at least one identical physical residue within the same model. It is a potential dependency, not measured statistical covariance or an effective sample size. Circular separation >= block length means no sampled block of that length can contain both columns. Local block sensitivity does not capture all feature-sharing dependencies. Counts across related taxa are not independent observations.',
        'artifacts': {'feature_overlap.tsv': sha(args.output / 'feature_overlap.tsv')}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
