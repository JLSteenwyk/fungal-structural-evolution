#!/usr/bin/env python3
"""Rebuild completed marker inputs and verify support splits by graph edge removal."""
import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from Bio import Phylo, SeqIO


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['snapshot', 'trees', 'matrix', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    sr = json.loads((a.snapshot / 'receipt.json').read_text())
    mr = json.loads((a.matrix / 'receipt.json').read_text())
    for folder, receipt in [(a.snapshot, sr), (a.matrix, mr)]:
        for name, digest in receipt['artifacts'].items():
            if sha(folder / name) != digest:
                raise ValueError('Changed source artifact')
    matrix = {r.id: str(r.seq) for r in SeqIO.parse(a.matrix / 'matrix.faa', 'fasta')}
    columns = defaultdict(list)
    for row in table(a.matrix / 'site_mapping.tsv'):
        columns[row['marker']].append(int(row['matrix_column_1based']) - 1)
    branches = defaultdict(dict)
    for row in table(a.snapshot / 'branch_support.tsv'):
        side = tuple(row['smaller_side_taxa'].split(';'))
        if side in branches[row['marker']]:
            raise ValueError('Duplicate snapshot split')
        branches[row['marker']][side] = row
    summaries = []
    characters = edges = 0
    for source in sr['inputs']:
        marker = source['marker']
        folder = a.trees / marker
        rp = folder / 'receipt.json'
        if sha(rp) != source['receipt_sha256']:
            raise ValueError('Tree receipt changed')
        receipt = json.loads(rp.read_text())
        threshold = max(50, (3 * len(columns[marker]) + 9) // 10)
        expected, excluded = {}, {}
        for taxon, seq in matrix.items():
            projected = ''.join(seq[i] for i in columns[marker])
            observed = sum(projected.count(c) for c in 'ACDEFGHIKLMNPQRSTVWY')
            if observed >= threshold:
                expected[taxon] = projected
            else:
                excluded[taxon] = observed
        actual = list(SeqIO.parse(folder / 'input.faa', 'fasta'))
        if len(actual) != len(expected) or {r.id: str(r.seq) for r in actual} != expected:
            raise ValueError('Input not exact frozen matrix projection under original coverage rule')
        if sha(folder / 'input.faa') != receipt['input_sha256']:
            raise ValueError('Input hash mismatch')
        cf = json.loads((folder / 'coverage_filter.json').read_text())
        if cf['minimum_unambiguous_residues'] != threshold or cf['excluded_taxa_observed_residues'] != excluded:
            raise ValueError('Coverage dispositions differ')
        command = receipt['command']
        for option, value in [('-m', 'MFP'), ('-mset', 'LG,WAG,JTT'), ('-mfreq', 'F'), ('-mrate', 'G'), ('--alrt', '1000')]:
            if command[command.index(option) + 1] != value:
                raise ValueError('Model/support specification differs')
        report = (folder / 'tree.iqtree').read_text()
        models = re.findall(r'^Model of substitution: (.+)$', report, re.M)
        if len(models) != 1 or models[0] not in {'LG+F+G4', 'WAG+F+G4', 'JTT+F+G4'}:
            raise ValueError('Unexpected fitted model')
        path = folder / 'tree.treefile'
        if sha(path) != source['tree_sha256']:
            raise ValueError('Tree changed')
        tree = Phylo.read(path, 'newick')
        adjacency = defaultdict(list)
        for node in tree.find_clades():
            for child in node.clades:
                adjacency[node].append(child)
                adjacency[child].append(node)
        universe = set(expected)
        if {x.name for x in tree.get_terminals()} != universe:
            raise ValueError('Terminal universe differs')
        checked = set()
        for parent in tree.find_clades():
            for child in parent.clades:
                if child.is_terminal():
                    continue
                # Traverse the undirected tree after deleting this edge.
                stack = [child]
                visited = {parent}
                names = set()
                while stack:
                    node = stack.pop()
                    if node in visited:
                        continue
                    visited.add(node)
                    if node.is_terminal():
                        names.add(node.name)
                    stack.extend(n for n in adjacency[node] if n not in visited)
                side = min([tuple(sorted(names)), tuple(sorted(universe - names))], key=lambda x: (len(x), x))
                row = branches[marker][side]
                digest = hashlib.sha256(json.dumps([sorted(universe), side]).encode()).hexdigest()
                if row['split_sha256'] != digest or int(row['smaller_side_size']) != len(side):
                    raise ValueError('Graph split identity differs')
                support = '' if child.confidence is None else str(child.confidence)
                if row['sh_alrt_percent'] != support or float(row['branch_length_substitutions_per_site']) != child.branch_length:
                    raise ValueError('Support/length association differs')
                checked.add(side)
        if checked != set(branches[marker]):
            raise ValueError('Incomplete split grid')
        edges += len(checked)
        characters += sum(map(len, expected.values()))
        summaries.append({'marker': marker, 'taxa': len(expected), 'columns': len(columns[marker]),
                          'excluded_taxa': len(excluded), 'model': models[0], 'internal_splits_checked': len(checked),
                          'report_warning_lines': sum('WARNING' in line for line in report.splitlines()),
                          'report_sha256': sha(folder / 'tree.iqtree'),
                          'coverage_filter_sha256': sha(folder / 'coverage_filter.json')})
    if len(summaries) != sr['completed_markers'] or edges != sr['internal_edges_audited']:
        raise ValueError('Snapshot totals differ')
    a.output.mkdir(parents=True)
    path = a.output / 'marker_readback.tsv'
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(summaries[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(summaries)
    result = {'status': 'passed_all_snapshot_input_and_graph_split_readbacks',
              'snapshot_receipt_sha256': sha(a.snapshot / 'receipt.json'),
              'matrix_receipt_sha256': sha(a.matrix / 'receipt.json'),
              'script_sha256': sha(Path(__file__)), 'markers': len(summaries),
              'retained_alignment_characters_checked': characters, 'internal_splits_checked': edges,
              'artifacts': {path.name: sha(path)},
              'interpretation': 'Every completed snapshot input reconstructed from original matrix and coverage exclusions; every internal support split independently recovered by graph-edge removal. ' + ('Full planned marker batch; no model adequacy or biological discordance-cause claim.' if sr['completed_markers'] == sr['planned_markers'] and not sr.get('pending_markers') else 'Incomplete, completion-order-biased marker subset; no full-batch discordance inference or model adequacy claim.')}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
