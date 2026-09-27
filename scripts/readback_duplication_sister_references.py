#!/usr/bin/env python3
"""Reconstruct all sister-reference choices with leaf intervals and upward path sums."""
import argparse
import csv
import json
import math
import sqlite3
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path
from Bio import Phylo
from run_ortholog_pair_guide_comparison import sha


def index_tree(newick):
    tree = Phylo.read(StringIO(newick), 'newick')
    nodes, parents, spans, leaves = {}, {}, {}, []
    stack = [(tree.root, None, False)]
    while stack:
        node, parent, exit_node = stack.pop()
        if exit_node:
            spans[node.name][1] = len(leaves)
            continue
        if not node.name or node.name in nodes:
            raise ValueError('Missing or duplicate node name')
        if parent is not None and (node.branch_length is None or
                not math.isfinite(node.branch_length) or node.branch_length < 0):
            raise ValueError('Invalid branch length')
        nodes[node.name] = node
        parents[node.name] = parent
        spans[node.name] = [len(leaves), None]
        stack.append((node, parent, True))
        if node.is_terminal():
            leaves.append(node.name)
        else:
            stack.extend((child, node.name, False) for child in reversed(node.clades))
    return nodes, parents, spans, leaves


def reconstruct(row, index, duplications, models):
    nodes, parents, spans, leaves = index
    name = row['gene_node']
    node = nodes[name]
    a, b = row['gene_a'], row['gene_b']
    if set(leaves[slice(*spans[name])]) != {a, b} or len(node.clades) != 2:
        raise ValueError('Candidate is not an exact pair')
    if parents[a] != name or parents[b] != name:
        raise ValueError('Candidate tips are not immediate children')
    parent = parents[name]
    result = dict(parent_node=parent or '', parent_reported_duplication=int(parent in duplications),
                  parent_children=0, sister_genes=0, sister_taxa=0,
                  sister_focal_taxon_genes=0, modeled_nonfocal_sister_genes=0,
                  modeled_nonfocal_sister_taxa=0, nearest_reference_genes=[],
                  chosen_reference_gene='', reference_model='', reference_version='',
                  reference_distance_from_duplicate_node='', distance_a_to_reference='',
                  distance_b_to_reference='',
                  duplicate_pair_sequence_distance=math.fsum([nodes[a].branch_length, nodes[b].branch_length]),
                  status='no_parent')
    if parent is None:
        return result
    # A leaf interval for the parent minus the duplicate interval, independent of
    # the producer's downward sister traversal and path accumulation.
    lo, hi = spans[parent]
    da, db = spans[name]
    sisters = leaves[lo:da] + leaves[db:hi]
    focal = row['taxon_id']
    taxon = {gene: gene.split('_', 1)[0] for gene in sisters}
    available = [gene for gene in sisters if taxon[gene] != focal and gene in models]
    result.update(parent_children=len(nodes[parent].clades), sister_genes=len(sisters),
                  sister_taxa=len(set(taxon.values())),
                  sister_focal_taxon_genes=sum(t == focal for t in taxon.values()),
                  modeled_nonfocal_sister_genes=len(available),
                  modeled_nonfocal_sister_taxa=len({taxon[g] for g in available}))
    if available:
        distances = {}
        for gene in available:
            edges = [node.branch_length]
            current = gene
            while current != parent:
                if current is None:
                    raise ValueError('Reference is outside parent clade')
                edges.append(nodes[current].branch_length)
                current = parents[current]
            distances[gene] = math.fsum(edges)
        best = min(distances.values())
        ties = sorted(g for g, d in distances.items() if abs(d-best) <= 1e-12)
        chosen = ties[0]
        distance = distances[chosen]
        model, version = models[chosen]
        result.update(nearest_reference_genes=ties, chosen_reference_gene=chosen,
                      reference_model=model, reference_version=version,
                      reference_distance_from_duplicate_node=distance,
                      distance_a_to_reference=math.fsum([distance, nodes[a].branch_length]),
                      distance_b_to_reference=math.fsum([distance, nodes[b].branch_length]))
    if parent in duplications:
        result['status'] = 'parent_reported_duplication'
    elif result['parent_children'] != 2:
        result['status'] = 'parent_not_bifurcating'
    elif result['sister_focal_taxon_genes']:
        result['status'] = 'focal_taxon_in_sister_clade'
    elif not available:
        result['status'] = 'no_modeled_nonfocal_sister'
    else:
        result['status'] = 'provisional_reference_available'
    return result


def check_row(actual, expected):
    max_delta = 0.0
    for name, value in expected.items():
        if isinstance(value, list):
            if json.loads(actual[name]) != value:
                raise ValueError('Reference ties differ')
        elif isinstance(value, float):
            observed = float(actual[name])
            if not math.isfinite(observed) or not math.isclose(observed, value, rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError('Sequence path differs: ' + name)
            max_delta = max(max_delta, abs(observed-value))
        elif str(value) != actual[name]:
            raise ValueError('Reference field differs: ' + name)
    return max_delta


def key(row):
    return row['family'], row['gene_node'], row['gene_a'], row['gene_b']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)
    def verify():
        if sha(args.plan) != plan_hash:
            raise ValueError('Changed readback plan')
        for path, digest in plan['pins'].items():
            if sha(path) != digest:
                raise ValueError('Changed source: ' + path)
    verify()
    sourceplan = json.loads(Path(plan['source_plan']).read_text())
    reviewplan = json.loads(Path(sourceplan['review_plan']).read_text())
    source = Path(sourceplan['output'])
    receipt = json.loads((source/'receipt.json').read_text())
    if receipt['status'] != 'complete_duplication_sister_reference_inventory' or receipt['plan_sha256'] != sha(plan['source_plan']):
        raise ValueError('Unbound source receipt')
    review = Path(sourceplan['review'])
    review_receipt = json.loads((review/'receipt.json').read_text())
    with sqlite3.connect('file:'+str(Path(sourceplan['bridge']).resolve())+'?mode=ro', uri=True) as conn:
        models = {gene: (model, version) for gene, model, version in conn.execute(
            "SELECT taxon_id || '_' || protein_id,model_id,version FROM structures")}
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    reports = []
    for entry in reviewplan['guides']:
        guide = entry['guide']
        original = review/(guide+'_candidate_tree_checks.tsv')
        observed = source/(guide+'_sister_references.tsv')
        if sha(original) != review_receipt['artifacts'][original.name] or sha(observed) != receipt['artifacts'][observed.name]:
            raise ValueError('Changed candidate/reference table')
        with original.open() as handle:
            candidates = list(csv.DictReader(handle, delimiter='\t'))
        with observed.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        mapping = {key(row): row for row in rows}
        if len(mapping) != len(rows) or len(rows) != len(candidates) or set(mapping) != {key(r) for r in candidates}:
            raise ValueError('Candidate universe differs')
        groups = defaultdict(list)
        for candidate in candidates:
            actual = mapping[key(candidate)]
            if any(actual[k] != v for k, v in candidate.items()):
                raise ValueError('Original candidate field changed')
            groups[candidate['family']].append(actual)
        del candidates, rows, mapping
        duplications = defaultdict(set)
        with Path(entry['events']).open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if row['Orthogroup'] in groups:
                    duplications[row['Orthogroup']].add(row['Gene Tree Node'])
        seen = set()
        counts = Counter()
        delta = 0.0
        with Path(entry['trees']).open() as handle:
            for line in handle:
                family, newick = line.rstrip().split(': ', 1)
                if family not in groups:
                    continue
                if family in seen:
                    raise ValueError('Repeated family tree')
                index = index_tree(newick)
                for actual in groups[family]:
                    expected = reconstruct(actual, index, duplications[family], models)
                    delta = max(delta, check_row(actual, expected))
                    counts[expected['status']] += 1
                seen.add(family)
                if len(seen) % 1000 == 0:
                    print(guide, 'checked families', len(seen), flush=True)
        if seen != set(groups):
            raise ValueError('Missing family trees')
        expected_counts = next(r for r in receipt['guides'] if r['guide'] == guide)
        if dict(counts) != expected_counts['counts'] or sum(counts.values()) != expected_counts['candidates']:
            raise ValueError('Disposition totals differ')
        report = dict(guide=guide, candidates=sum(counts.values()), families=len(seen),
                      counts=dict(counts), maximum_sequence_distance_difference=delta)
        (out/(guide+'.json')).write_text(json.dumps(report, indent=2)+'\n')
        reports.append(report)
        print(json.dumps(report), flush=True)
    verify()
    result = dict(status='passed_full_duplication_sister_reference_readback',
                  plan_sha256=plan_hash, producer_receipt_sha256=sha(source/'receipt.json'),
                  guides=reports, artifacts={p.name: sha(p) for p in out.iterdir()},
                  scope='Every candidate and all exported sister/reference fields checked. '
                  'Leaf-interval subtraction and upward edge sums independently reconstruct '
                  'sister membership, nearest ties, representative/model and sequence distances. '
                  'Shares Bio.Phylo Newick parser, not producer selection functions. '
                  'Fixed gene trees and native duplication calls are inputs, not biological '
                  'truth. No independent orthology, rooting, event timing or structural inference.')
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    csv.field_size_limit(32*1024*1024)
    main()
