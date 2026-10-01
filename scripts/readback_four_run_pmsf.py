#!/usr/bin/env python3
"""Independently exhaust all four crossed PMSF split/support/conflict exports."""
import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import dendropy
from four_run_pmsf_sources import load_sources, PAIR_FIELDS, PRESENCE_FIELDS, CONFLICT_FIELDS, BOUNDARY_FIELDS, SUMMARY_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def read(path, expected_fields):
    with Path(path).open() as f:
        reader = csv.DictReader(f, delimiter='\t'); assert reader.fieldnames == expected_fields
        return list(reader)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text()); sources, manifest, bindings = load_sources(plan, args.plan)
    root = Path(plan['output']); rp = root / 'receipt.json'; r = json.loads(rp.read_text()); rh = sha(rp)
    assert r['status'] == 'complete_full_four_run_pmsf_ML_and_consensus_sensitivity_pending_readback' and r['plan_sha256'] == sha(args.plan)
    assert r['source_hashes'] == bindings
    for name, digest in r['artifacts'].items(): assert sha(root / name) == digest
    universe = {row['taxon_id'] for row in manifest}; ns = dendropy.TaxonNamespace(); edges = {}; supports = {}; kinds = {}; labels = {}
    def key(side):
        a, b = tuple(sorted(side)), tuple(sorted(universe.difference(side)))
        return a if (len(a), a) <= (len(b), b) else b
    for label, source in sources.items():
        for kind, file in [('ml', 'pmsf.treefile'), ('consensus', 'pmsf.contree')]:
            name = label + ':' + kind
            tree = dendropy.Tree.get(path=str(Path(source['spec']['run']) / file), schema='newick', rooting='force-unrooted', preserve_underscores=True, taxon_namespace=ns)
            assert len(list(tree.leaf_node_iter())) == len(universe) and {n.taxon.label for n in tree.leaf_node_iter()} == universe
            tree.encode_bipartitions(); lookup = {}
            for node in tree.preorder_node_iter():
                if node is tree.seed_node or node.is_leaf(): continue
                split = key({n.taxon.label for n in node.leaf_iter()}); assert split not in lookup; lookup[split] = node
            assert len(lookup) == 523
            edges[name] = lookup; supports[name] = {}; kinds[name] = kind; labels[name] = label
            for row in source['rows']:
                if row['tree'] != kind: continue
                split = key(set(json.loads(row['split_taxa_json']))); assert split not in supports[name]
                node = lookup[split]
                if kind == 'ml':
                    alrt, ufb = map(float, node.label.split('/')); assert alrt == float(row['sh_alrt_percent'])
                else:
                    ufb = float(node.label); assert row['sh_alrt_percent'] == ''
                assert abs(ufb - float(row['empirical_ufboot_percent'])) <= .500001
                assert float(row['reported_ufboot_percent']) == ufb and float(row['branch_length']) == node.edge.length
                supports[name][split] = row
            assert set(supports[name]) == set(lookup)
    union = set().union(*(set(v) for v in edges.values())); presence = read(root / 'split_presence.tsv', PRESENCE_FIELDS); seen = set()
    for row in presence:
        name = row['view']; split = tuple(json.loads(row['split_taxa_json'])); assert name in edges and split in union
        assert (name, split) not in seen; seen.add((name, split))
        assert row['run'] == labels[name] and row['tree_type'] == kinds[name]
        present = split in edges[name]; assert row['present'] == str(present)
        for field, source_field in [('branch_length', 'branch_length'), ('sh_alrt', 'sh_alrt_percent'), ('empirical_ufb', 'empirical_ufboot_percent')]:
            assert row[field] == (supports[name][split][source_field] if present else '')
    assert seen == {(name, split) for name in edges for split in union}
    rule = lambda name: 'SH_aLRT80_and_empirical_UFB95' if kinds[name] == 'ml' else 'empirical_UFB95_SH_aLRT_unavailable'
    high = lambda name, split: float(supports[name][split]['empirical_ufboot_percent']) >= 95 and (kinds[name] == 'consensus' or float(supports[name][split]['sh_alrt_percent']) >= 80)
    conflicts = read(root / 'conflicts.tsv', CONFLICT_FIELDS)
    indexed = {(v['view_a'], v['view_b'], tuple(json.loads(v['split_a_taxa_json'])), tuple(json.loads(v['split_b_taxa_json']))): v for v in conflicts}
    assert len(indexed) == len(conflicts)
    pair_grid = [(a, b) for a, b in itertools.combinations(edges, 2) if kinds[a] == kinds[b] or labels[a] == labels[b]]
    pairs = read(root / 'comparisons.tsv', PAIR_FIELDS); assert len(pairs) == len(pair_grid) == 16
    assert [(v['view_a'], v['view_b']) for v in pairs] == pair_grid
    expected_conflicts = set(); expected_pairs = []
    for actual, (a, b) in zip(pairs, pair_grid):
        left, right = edges[a], edges[b]; common = set(left) & set(right); nconf = nhigh = 0
        for x, y in itertools.product(set(left) - common, set(right) - common):
            if left[x].bipartition.is_compatible_with(right[y].bipartition): continue
            state = a, b, x, y; expected_conflicts.add(state); row = indexed[state]; nconf += 1
            assert row['support_rule_a'] == rule(a) and row['support_rule_b'] == rule(b)
            for suffix, view, split in [('a', a, x), ('b', b, y)]:
                assert row['sh_alrt_' + suffix] == supports[view][split]['sh_alrt_percent']
                assert row['empirical_ufb_' + suffix] == supports[view][split]['empirical_ufboot_percent']
            supported = high(a, x) and high(b, y); nhigh += supported; assert row['both_support_criteria_met'] == str(supported)
            quartet = row['witness_quartet'].split(';'); assert len(set(quartet)) == 4 and set(quartet) <= universe
            assert {(taxon in x, taxon in y) for taxon in quartet} == {(False, False), (False, True), (True, False), (True, True)}
            assert quartet == [min(set(x) & set(y)), min(set(x) - set(y)), min(set(y) - set(x)), min(universe - (set(x) | set(y)))]
        distance = len(left) + len(right) - 2 * len(common)
        wanted = dict(view_a=a, view_b=b, comparison_scope='within_run_ML_vs_consensus' if labels[a] == labels[b] else 'cross_run_' + kinds[a],
                      support_rule_a=rule(a), support_rule_b=rule(b), shared_internal_splits=len(common), unique_internal_splits_a=len(left) - len(common), unique_internal_splits_b=len(right) - len(common),
                      rf_distance=distance, normalized_rf=distance / (len(left) + len(right)), incompatible_split_pairs=nconf, both_supported_incompatible_pairs=nhigh)
        for field, value in wanted.items():
            if isinstance(value, float): assert math.isclose(float(actual[field]), value, rel_tol=1e-12, abs_tol=1e-12)
            else: assert actual[field] == str(value), field
        expected_pairs.append(wanted)
    assert expected_conflicts == set(indexed) and r['comparisons'] == expected_pairs
    role_split = key({row['taxon_id'] for row in manifest if row['study_role'] == 'outgroup'})
    boundary = read(root / 'rooting_boundary.tsv', BOUNDARY_FIELDS); assert len(boundary) == 8 and [v['view'] for v in boundary] == list(edges)
    present_views = 0
    for row in boundary:
        name = row['view']; present = role_split in edges[name]; present_views += present
        assert row['ingroup_count'] == '501' and row['outgroup_count'] == '25' and row['boundary_split_present'] == str(present)
        assert tuple(json.loads(row['boundary_taxa_json'])) == role_split
        assert row['sh_alrt'] == (supports[name][role_split]['sh_alrt_percent'] if present else '')
        assert row['empirical_ufb'] == (supports[name][role_split]['empirical_ufboot_percent'] if present else '')
    summary = dict(taxa=526, runs=4, tree_views=8, internal_splits_per_view=523, comparison_rows=len(pairs), split_presence_rows=len(presence), incompatible_pairs=len(expected_conflicts),
                   shared_all_ml=len(set.intersection(*(set(v) for name, v in edges.items() if kinds[name] == 'ml'))),
                   shared_all_consensus=len(set.intersection(*(set(v) for name, v in edges.items() if kinds[name] == 'consensus'))),
                   shared_all_eight=len(set.intersection(*(set(v) for v in edges.values()))), boundary_views_with_role_split=present_views)
    assert all(r[field] == summary[field] for field in SUMMARY_FIELDS)
    verify(bindings); assert sha(rp) == rh
    for name, digest in r['artifacts'].items(): assert sha(root / name) == digest
    proof = dict(status='passed_full_four_run_pmsf_ML_consensus_split_support_readback', plan_sha256=sha(args.plan), producer_receipt_sha256=rh,
                 **summary, source_hashes=bindings, dendropy_version=dendropy.__version__, scientific_eligibility=False,
                 scope='Independent DendroPy raw-tree parsing, node labels/lengths and bipartition compatibility reconstruct all eight523split views, union/missing cells,16comparison rows, every incompatible pair and canonical quartet, tree-specific support criteria and full501/25role split. Bound upstream1000bootstrap/profile audits; no SH-aLRT rerun, root assignment, adequacy/biological discordance or accepted species framework.')
    with args.output.open('x') as f: f.write(json.dumps(proof, indent=2) + '\n')
    print(json.dumps({k: v for k, v in proof.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
