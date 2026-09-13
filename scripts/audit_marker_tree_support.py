#!/usr/bin/env python3
"""Audit completed marker trees and record root-independent supported splits."""
import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from Bio import Phylo, SeqIO
from build_species_matrix import ROOT, sha


def supported_splits(tree):
    names = [n.name for n in tree.get_terminals()]
    if any(n is None for n in names) or len(set(names)) != len(names):
        raise ValueError('Missing or duplicate terminal identities')
    universe = frozenset(names)
    if len(tree.root.clades) != 3:
        raise ValueError('Expected a trifurcating representation of an unrooted IQ-TREE tree')
    rows, seen = [], set()
    for node in tree.find_clades():
        if node is tree.root:
            continue
        if node.branch_length is None or not math.isfinite(node.branch_length) or node.branch_length < 0:
            raise ValueError('Missing, nonfinite or negative branch length')
        if node.is_terminal():
            continue
        a = tuple(sorted(n.name for n in node.get_terminals()))
        b = tuple(sorted(universe - set(a)))
        side = min([a, b], key=lambda s: (len(s), s))
        if len(side) < 2 or side in seen:
            raise ValueError('Trivial or duplicate internal split')
        seen.add(side)
        support = node.confidence
        if support is not None and (not math.isfinite(support) or not 0 <= support <= 100):
            raise ValueError('Invalid SH-aLRT support')
        split_id = hashlib.sha256(json.dumps([sorted(universe), side]).encode()).hexdigest()
        rows.append({'split_sha256': split_id, 'smaller_side_taxa': ';'.join(side),
                     'smaller_side_size': len(side), 'sh_alrt_percent': support if support is not None else '',
                     'support_status': 'reported' if support is not None else 'not_reported',
                     'branch_length_substitutions_per_site': node.branch_length})
    if len(rows) != len(names) - 3:
        raise ValueError('Expected n-3 internal edges for a resolved unrooted tree')
    return sorted(rows, key=lambda r: r['split_sha256'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trees', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--allow-incomplete', action='store_true')
    args = parser.parse_args()
    matrix = ROOT / 'results/phylogeny/profile-matrix-50-v1'
    mr = json.loads((matrix / 'receipt.json').read_text())
    if sha(matrix / 'site_mapping.tsv') != mr['artifacts']['site_mapping.tsv']:
        raise ValueError('Changed expected marker map')
    with (matrix / 'site_mapping.tsv').open() as handle:
        expected = {r['marker'] for r in csv.DictReader(handle, delimiter='\t')}
    completed = {}
    for path in sorted(args.trees.glob('*/receipt.json')):
        r = json.loads(path.read_text())
        if r.get('status') == 'inferred':
            if r['marker'] not in expected or r['marker'] in completed:
                raise ValueError('Unexpected or duplicate completed marker')
            completed[r['marker']] = (path, r)
    pending = sorted(expected - completed.keys())
    if pending and not args.allow_incomplete:
        raise ValueError('Full marker-tree batch is not complete')
    if not completed:
        raise ValueError('No completed trees to audit')
    if args.output.exists():
        raise FileExistsError('Use a new immutable audit snapshot')
    summaries, branches, inputs = [], [], []
    for marker, (path, receipt) in sorted(completed.items()):
        if receipt['returncode'] != 0 or receipt['matrix_receipt_sha256'] != sha(matrix / 'receipt.json'):
            raise ValueError('Tree completion or source provenance mismatch')
        command = receipt['command']
        if command[command.index('--alrt') + 1] != '1000':
            raise ValueError('Unexpected support procedure')
        tree_path, fasta = ROOT / receipt['tree_path'], path.parent / 'input.faa'
        if sha(tree_path) != receipt['tree_sha256'] or sha(fasta) != receipt['input_sha256']:
            raise ValueError('Changed tree or alignment')
        tree = Phylo.read(tree_path, 'newick')
        taxa = [r.id for r in SeqIO.parse(fasta, 'fasta')]
        if len(set(taxa)) != len(taxa) or set(taxa) != {n.name for n in tree.get_terminals()} or len(taxa) != receipt['taxa']:
            raise ValueError('Tree taxa differ from inference input')
        splits = supported_splits(tree)
        supports = [r['sh_alrt_percent'] for r in splits if r['support_status'] == 'reported']
        row = {'marker': marker, 'taxa': len(taxa), 'alignment_columns': receipt['columns'],
               'internal_edges': len(splits), 'edges_with_reported_support': len(supports),
               'edges_without_reported_support': len(splits) - len(supports),
               'median_reported_sh_alrt_percent': statistics.median(supports) if supports else ''}
        for cutoff in [50, 80, 95]:
            row[f'edges_sh_alrt_ge{cutoff}'] = sum(s >= cutoff for s in supports)
        summaries.append(row)
        branches.extend(dict(marker=marker, **r) for r in splits)
        inputs.append({'marker': marker, 'receipt_sha256': sha(path), 'tree_sha256': sha(tree_path)})
    args.output.mkdir(parents=True)
    for name, rows in [('marker_support.tsv', summaries), ('branch_support.tsv', branches)]:
        with (args.output / name).open('w') as handle:
            writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    receipt = {'status': 'incomplete_snapshot' if pending else 'complete_support_audit',
        'planned_markers': len(expected), 'completed_markers': len(completed), 'pending_markers': pending,
        'internal_edges_audited': len(branches), 'inputs': inputs, 'script_sha256': sha(Path(__file__)),
        'interpretation': 'SH-aLRT from 1000 replicates, not bootstrap percentages or posterior probabilities. Split identities include each tree taxon universe. No rooting, time calibration, species-tree reconciliation or independent transition inference.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['pending_markers','inputs']}, indent=2))


if __name__ == '__main__':
    main()
