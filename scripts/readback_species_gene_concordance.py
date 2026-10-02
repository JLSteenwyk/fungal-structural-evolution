#!/usr/bin/env python3
"""Independently reconstruct every native gCF cell from DendroPy gene splits.

The native producer uses IQ-TREE, not this classification algorithm. The four
incident clades determine decisiveness; observed restricted gene bipartitions
determine concordance, either NNI alternative, or residual discordance. Native
NNI names may differ by orientation, but one consistent bijection is required
for every gene around each reference branch. Missing evidence remains NA.
"""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path

import dendropy
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha
from species_gcf_sources import load


STAT_FIELDS = ['ID', 'gCF', 'gCF_N', 'gDF1', 'gDF1_N', 'gDF2',
               'gDF2_N', 'gDFP', 'gDFP_N', 'gN', 'Label', 'Length']
CELL_FIELDS = ['ID', 'TreeID', 'gC', 'gD1', 'gD2']
PERCENT_FIELDS = [('gCF', 'gCF_N'), ('gDF1', 'gDF1_N'),
                  ('gDF2', 'gDF2_N'), ('gDFP', 'gDFP_N')]


def table(path, fields):
    with Path(path).open() as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith('#')),
                                delimiter='\t')
        assert reader.fieldnames == fields, (path, reader.fieldnames)
        rows = list(reader)
        assert all(set(row) == set(fields) and None not in row.values() for row in rows)
        return rows


def canonical(mask, universe):
    assert mask & universe == mask
    return min(mask, universe ^ mask)


def parse_tree(text, index):
    tree = dendropy.Tree.get(data=text, schema='newick', rooting='force-unrooted',
                             preserve_underscores=True)
    nodes = list(tree.preorder_node_iter())
    graph = {node: [] for node in nodes}
    observed = [node.taxon.label for node in tree.leaf_node_iter()]
    assert len(observed) == len(set(observed)) >= 4 and set(observed) <= set(index)
    universe = sum(1 << index[label] for label in observed)
    desc = {}
    for node in tree.postorder_node_iter():
        desc[node] = (1 << index[node.taxon.label] if node.is_leaf()
                      else sum(desc[child] for child in node.child_node_iter()))
        for child in node.child_node_iter():
            graph[node].append(child)
            graph[child].append(node)
    assert all(len(neighbors) == (1 if node.is_leaf() else 3)
               for node, neighbors in graph.items()), 'Expected binary unrooted tree'
    edges = {}
    for node in nodes:
        if node.parent_node is None:
            continue
        key = canonical(desc[node], universe)
        assert key not in edges
        edges[key] = dict(node=node, parent=node.parent_node, label=node.label,
                          length=Decimal(str(node.edge_length)))
    assert len(edges) == 2 * len(observed) - 3
    return dict(tree=tree, graph=graph, edges=edges, universe=universe,
                internal={key: value for key, value in edges.items()
                          if key.bit_count() > 1 and (universe ^ key).bit_count() > 1})


def component(graph, node, forbidden, index):
    stack = [(node, forbidden)]
    mask = 0
    while stack:
        current, parent = stack.pop()
        if current.is_leaf():
            mask |= 1 << index[current.taxon.label]
        stack.extend((neighbor, current) for neighbor in graph[current] if neighbor is not parent)
    return mask


def incident_clades(reference, edge, index):
    left, right = edge['node'], edge['parent']
    graph = reference['graph']
    sides = [sorted(component(graph, neighbor, endpoint, index)
                    for neighbor in graph[endpoint] if neighbor is not opposite)
             for endpoint, opposite in [(left, right), (right, left)]]
    assert all(len(side) == 2 for side in sides)
    sides.sort(key=lambda side: side[0] | side[1])
    clades = sides[0] + sides[1]
    assert sum(clades) == reference['universe'] and all(clades)
    return clades


def percentage(value, numerator, denominator):
    if denominator == 0:
        assert numerator == 0 and value == 'NA'
    else:
        expected = Decimal(100) * Decimal(numerator) / Decimal(denominator)
        # IQ-TREE prints concordance percentages to two decimal places.
        assert value != 'NA' and abs(Decimal(value) - expected) <= Decimal('0.00501'), (value, expected)


def validate_run(folder, reference_text, gene_texts, taxon_names, allow_reference_labels=False):
    """Validate all native exports, yielding independently classified cells."""
    folder = Path(folder)
    index = {label: i for i, label in enumerate(sorted(taxon_names))}
    reference = parse_tree(reference_text, index)
    assert reference['universe'].bit_count() == len(index)
    if not allow_reference_labels:
        assert all(edge['label'] in [None, ''] for edge in reference['internal'].values())
    native_branches = parse_tree((folder / 'gcf.cf.branch').read_text(), index)
    annotated = parse_tree((folder / 'gcf.cf.tree').read_text(), index)
    for native in [native_branches, annotated]:
        assert native['universe'] == reference['universe'] and set(native['edges']) == set(reference['edges'])
        assert all(abs(native['edges'][key]['length'] - edge['length']) <= Decimal('0.0000000001')
                   for key, edge in reference['edges'].items())
    ids = {int(edge['label']): key for key, edge in native_branches['internal'].items()}
    assert len(ids) == len(reference['internal'])
    stat_rows = table(folder / 'gcf.cf.stat', STAT_FIELDS)
    stats = {int(row['ID']): row for row in stat_rows}
    assert len(stats) == len(stat_rows) and set(stats) == set(ids)
    raw_cells = table(folder / 'gcf.cf.stat_tree', CELL_FIELDS)
    cells = {(int(row['ID']), int(row['TreeID'])): row for row in raw_cells}
    assert len(cells) == len(raw_cells) == len(ids) * len(gene_texts)
    assert set(cells) == {(branch, i) for branch in ids for i in range(1, len(gene_texts) + 1)}
    genes = [parse_tree(text, index) for text in gene_texts]
    ledger, summaries = [], []
    for branch, key in sorted(ids.items()):
        clades = incident_clades(reference, reference['internal'][key], index)
        assert canonical(clades[0] | clades[1], reference['universe']) == key
        alternatives = [canonical(clades[0] | clades[i], reference['universe']) for i in [2, 3]]
        assert len(set(alternatives + [key])) == 3
        mappings = [(0, 1), (1, 0)]
        counts = Counter()
        for i, gene in enumerate(genes, 1):
            covered = [mask & gene['universe'] for mask in clades]
            native = cells[branch, i]
            native_values = [native[field] for field in ['gC', 'gD1', 'gD2']]
            if not all(covered):
                state = 'not_decisive_missing_incident_clade'
                assert native_values == ['NA'] * 3, (branch, i, native_values)
            else:
                restricted = [canonical(split & gene['universe'], gene['universe'])
                              for split in [key] + alternatives]
                assert len(set(restricted)) == 3
                found = [split in gene['internal'] for split in restricted]
                assert sum(found) <= 1
                state = ['concordant', 'alternative_1', 'alternative_2'][found.index(True)] if any(found) else 'residual_discordance'
                assert set(native_values) <= {'0', '1'} and sum(map(int, native_values)) <= 1
                assert int(native['gC']) == int(found[0]), (branch, i, state, native_values)
                mappings = [mapping for mapping in mappings
                            if int(native['gD1']) == int(found[1 + mapping[0]])
                            and int(native['gD2']) == int(found[1 + mapping[1]])]
                assert mappings, ('Native NNI orientation inconsistent', branch, i, state, native_values)
            counts[state] += 1
            ledger.append(dict(branch_id=branch, canonical_split_mask=hex(key), tree_id=i,
                               clade_taxa_counts=[mask.bit_count() for mask in covered], state=state,
                               native_gC=native['gC'], native_gD1=native['gD1'], native_gD2=native['gD2']))
        row = stats[branch]
        n = len(genes) - counts['not_decisive_missing_incident_clade']
        mapping = mappings[0]
        if n == 0:
            # This IQ-TREE build prints NA even for all count fields when no
            # gene is decisive, and leaves the Newick branch unlabelled.
            assert all(row[field] == 'NA' for field in STAT_FIELDS[1:10])
        else:
            assert int(row['gN']) == n and int(row['gCF_N']) == counts['concordant'] and int(row['gDFP_N']) == counts['residual_discordance']
            assert int(row['gDF1_N']) == counts['alternative_' + str(mapping[0] + 1)]
            assert int(row['gDF2_N']) == counts['alternative_' + str(mapping[1] + 1)]
            assert sum(int(row[count]) for _, count in PERCENT_FIELDS) == n
            for pct, count in PERCENT_FIELDS:
                percentage(row[pct], int(row[count]), n)
        if not allow_reference_labels:
            assert row['Label'] == ''
            assert annotated['internal'][key]['label'] == (row['gCF'] if n else None)
        assert abs(Decimal(row['Length']) - reference['internal'][key]['length']) <= Decimal('0.0000000001')
        split_taxa = [label for label, bit in index.items() if key & (1 << bit)]
        split_id = hashlib.sha256(json.dumps(split_taxa, separators=(',', ':')).encode()).hexdigest()
        summaries.append(dict(branch_id=branch, canonical_split_id=split_id, canonical_split_mask=hex(key),
                              split_taxa=split_taxa, four_clades=[[label for label, bit in index.items() if mask & (1 << bit)] for mask in clades],
                              canonical_alternative_masks=list(map(hex, alternatives)),
                              native_nni_mapping=[v + 1 for v in mapping],
                              native_nni_mapping_identified=len(mappings) == 1, marker_count=len(genes), decisive=n,
                              concordant=counts['concordant'], alternative_1=counts['alternative_1'],
                              alternative_2=counts['alternative_2'], residual_discordance=counts['residual_discordance'],
                              not_decisive=counts['not_decisive_missing_incident_clade'],
                              native_statistics=row))
    return summaries, ledger


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    sources, universe, markers, bindings = load(plan, args.plan)
    root = Path(plan['output'])
    receipt_path = root / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'complete_full_native_species_gene_concordance_batch_pending_independent_readback'
    assert receipt['plan_sha256'] == sha(args.plan) and receipt['scientific_eligibility'] is False
    assert receipt['references'] == 8 and receipt['native_runs'] == 16 and receipt['markers_per_alignment'] == 125
    assert receipt['marker_alignments'] == 2 and receipt['branch_marker_cells'] == 1046000
    for path, digest in receipt['source_hashes'].items(): bind(bindings, path, digest)
    bind(bindings, receipt_path)
    bind(bindings, root / 'config.json')
    config = json.loads((root / 'config.json').read_text())
    assert config['plan_sha256'] == sha(args.plan)
    for path, digest in config['source_hashes'].items(): bind(bindings, path, digest)
    expected_grid = {(run, kind, alignment) for run in sources for kind in ['ml', 'consensus'] for alignment in markers}
    runs = {(run['species_run'], run['tree_type'], run['marker_alignment']): run for run in receipt['runs']}
    assert len(runs) == len(receipt['runs']) == 16 and set(runs) == expected_grid
    gene_texts = {alignment: [Path(item['path']).read_text().strip() for item in items] for alignment, items in markers.items()}
    for alignment, items in markers.items():
        assert (root / 'inputs' / (alignment + '.genes.tree')).read_text() == ''.join(text + '\n' for text in gene_texts[alignment])
        assert json.loads((root / 'inputs' / (alignment + '.gene_order.json')).read_text()) == items
    verify(bindings)
    summary_path = root / 'independent_branch_concordance.jsonl.gz'
    cell_path = root / 'independent_branch_marker_states.jsonl.gz'
    total_cells = total_branches = 0
    totals = Counter()
    run_summaries = []
    with gzip.open(summary_path, 'xt') as sh, gzip.open(cell_path, 'xt') as ch:
        for grid, run in runs.items():
            label, kind, alignment = grid
            name = label + '-' + kind + '-' + alignment
            assert run['label'] == name and run['returncode'] == 0 and run['marker_count'] == 125
            assert run['taxa'] == 526 and run['expected_internal_branches'] == 523
            folder = root / name
            reference_path = root / 'inputs' / (label + '.' + kind + '.reference.tree')
            assert run['reference'] == str(reference_path)
            assert run['genes'] == str(root / 'inputs' / (alignment + '.genes.tree'))
            assert run['gene_order'] == str(root / 'inputs' / (alignment + '.gene_order.json'))
            original = Path(sources[label]['spec']['run']) / ('pmsf.treefile' if kind == 'ml' else 'pmsf.contree')
            index = {taxon: i for i, taxon in enumerate(sorted(universe))}
            raw = parse_tree(original.read_text(), index)
            prepared = parse_tree(reference_path.read_text(), index)
            assert set(raw['edges']) == set(prepared['edges'])
            assert all(abs(raw['edges'][key]['length'] - edge['length']) <= Decimal('0.0000000001')
                       for key, edge in prepared['edges'].items())
            command = [plan['executable'], '-t', str(reference_path), '--gcf', run['genes'], '--cf-verbose',
                       '-T', str(plan['resources']['cpu']), '--mem', '8G', '--seed', str(plan['seed']), '--prefix', str(folder / 'gcf')]
            assert run['native_command'] == command
            launched = json.loads((folder / 'native_launch.json').read_text())
            assert launched['command'] == launched['cmdline'] == command and launched['plan_sha256'] == sha(args.plan)
            rp = folder / 'receipt.json'
            r = json.loads(rp.read_text())
            assert r == dict(status='complete_native_gene_concordance_pending_full_independent_readback', **run)
            assert set(run['artifacts']) == {'stdout.log', 'native_launch.json', 'gcf.cf.stat', 'gcf.cf.stat_tree', 'gcf.cf.branch', 'gcf.cf.tree', 'gcf.log'}
            for file, digest in run['artifacts'].items(): bind(bindings, folder / file, digest)
            bind(bindings, rp)
            summaries, cells = validate_run(folder, reference_path.read_text(), gene_texts[alignment], universe)
            assert len(summaries) == 523 and len(cells) == 65375
            for row in summaries:
                key = int(row['canonical_split_mask'], 16)
                row.update(run=name, species_run=label, tree_type=kind, marker_alignment=alignment,
                           original_support_label=raw['internal'][key]['label'],
                           original_branch_length=str(raw['internal'][key]['length']))
                sh.write(json.dumps(row, separators=(',', ':')) + '\n')
            counts = Counter(row['state'] for row in cells)
            for row in cells:
                row.update(run=name, marker=markers[alignment][row['tree_id'] - 1]['marker'])
                ch.write(json.dumps(row, separators=(',', ':')) + '\n')
            run_summaries.append(dict(run=name, branches=len(summaries), branch_marker_cells=len(cells),
                                      state_counts=dict(counts), zero_decisive_branches=sum(row['decisive'] == 0 for row in summaries)))
            totals.update(counts)
            total_branches += len(summaries)
            total_cells += len(cells)
            print('independent_native_gcf_readback_complete', name, len(run_summaries), '/16', flush=True)
    assert total_branches == 8368 and total_cells == 1046000
    for path in [summary_path, cell_path]: bind(bindings, path)
    verify(bindings)
    result = dict(status='passed_full_native_species_gene_concordance_independent_readback',
                  plan_sha256=sha(args.plan), producer_receipt_sha256=sha(receipt_path),
                  references=8, marker_alignments=2, markers_per_alignment=125, native_runs=16,
                  branch_marker_cells=total_cells, branch_summary_rows=total_branches,
                  state_counts=dict(totals), run_summaries=run_summaries,
                  artifacts={path.name: sha(path) for path in [summary_path, cell_path]},
                  source_hashes=bindings, scientific_eligibility=False, scope=plan['scope'])
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'native_runs', 'branch_marker_cells', 'branch_summary_rows', 'state_counts']}), flush=True)


if __name__ == '__main__':
    main()
