#!/usr/bin/env python3
"""Independently audit rooted guide inputs by deleting edges in undirected graphs."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
from Bio import Phylo


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def graph_edges(path, mapping=None):
    tree = Phylo.read(path, 'newick')
    nodes = list(tree.find_clades())
    index = {n: i for i, n in enumerate(nodes)}
    leaves = {index[n]: mapping[n.name] if mapping else n.name for n in tree.get_terminals()}
    if len(leaves) != 526 or len(set(leaves.values())) != 526:
        raise ValueError('Taxon identity count mismatch')
    universe = set(leaves.values())
    adjacent = defaultdict(list)
    edges = []
    for node in nodes:
        for child in node.clades:
            u, v = index[node], index[child]
            length = child.branch_length
            if length is None or not math.isfinite(length) or length < 0:
                raise ValueError('Invalid branch length')
            adjacent[u].append(v); adjacent[v].append(u)
            edges.append((u, v, length))
    result = defaultdict(float)
    for u, v, length in edges:
        visited = {u}; pending = [v]; side = set()
        while pending:
            x = pending.pop()
            if x in visited:
                continue
            visited.add(x)
            if x in leaves:
                side.add(leaves[x])
            pending.extend(adjacent[x])
        key = min([tuple(sorted(side)), tuple(sorted(universe - side))], key=lambda x: (len(x), x))
        result[key] += length
    return tree, universe, dict(result)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    plan = json.loads(a.plan.read_text())
    folder = Path(plan['output'])
    receipt = json.loads((folder / 'receipt.json').read_text())
    assert receipt['plan_sha256'] == sha(a.plan)
    assert all(sha(Path(p)) == h for p, h in plan['pins'].items())
    with Path(plan['manifest']).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    roles = {row['taxon_id']: row['study_role'] for row in rows}
    expected_root = {frozenset(k for k, v in roles.items() if v == role) for role in ['ingroup', 'outgroup']}
    mapping = {}
    for line in Path(plan['species_ids']).read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            native, filename = line.split(': ', 1)
            assert native not in mapping
            mapping[native] = filename.rsplit('.', 1)[0]
    assert len(mapping) == len(set(mapping.values())) == 526
    assert json.loads((folder / 'species_id_mapping.json').read_text()) == mapping
    assert sha(folder / 'species_id_mapping.json') == receipt['mapping_sha256']
    checks = []
    for item in receipt['guides']:
        guide = item['guide']
        _, universe, original = graph_edges(Path(plan['guides'][guide]))
        assert universe == set(roles) and len(original) == 1049
        for filename in ['species_tree_taxa.nwk', 'species_tree_ids.nwk', 'native_converter_check.nwk']:
            path = folder / guide / filename
            assert sha(path) == item['artifacts'][filename]
            ids = None if filename == 'species_tree_taxa.nwk' else mapping
            tree, taxa, current = graph_edges(path, ids)
            assert taxa == universe and current.keys() == original.keys()
            assert all(len(n.clades) in [0, 2] for n in tree.find_clades())
            root_partition = {frozenset(ids[t.name] if ids else t.name for t in n.get_terminals()) for n in tree.root.clades}
            assert root_partition == expected_root
            errors = [abs(current[k] - original[k]) for k in original]
            if filename == 'native_converter_check.nwk':
                assert all(abs(current[k] - original[k]) <= 1e-5 * max(original[k], 1e-6) for k in original)
            else:
                assert max(errors) < 1e-12
                assert all(abs(n.branch_length - item['root_edge_original_length'] / 2) < 1e-12 for n in tree.root.clades)
            checks.append(dict(guide=guide, file=filename, taxa_checked=526,
                               canonical_edges_checked=len(original), max_edge_length_error=max(errors)))
    assert len(checks) == 6
    result = dict(status='passed_all_rooted_guide_graph_and_id_readbacks', checks=checks,
                  source_receipt_sha256=sha(folder / 'receipt.json'), script_sha256=sha(Path(__file__)),
                  scope='All six output trees: independent undirected edge deletion, exact native-ID mapping, 501/25 root partition, full bifurcation and all 6294 canonical edge comparisons. Intended named/IDs inputs preserve lengths to 1e-12; native converter rounding is separately quantified. Does not establish a final supported species tree or biological root.')
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print('Verified six rooted trees, exact species IDs and 6294 canonical edge comparisons')


if __name__ == '__main__':
    main()
