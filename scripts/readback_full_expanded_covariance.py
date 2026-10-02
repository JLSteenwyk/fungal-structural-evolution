#!/usr/bin/env python3
"""Rebuild every expanded reuse link, family component and tree covariance."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import dendropy
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from full_expanded_covariance_sources import load, identity, CASE_FIELDS, INCIDENCE_FIELDS, PATTERN_FIELDS, FAMILY_FIELDS, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(path, output):
    path = Path(path); plan = json.loads(path.read_text()); source, bindings = load(plan, path); original = dict(bindings)
    root = Path(plan['output']); rp = root / 'receipt.json'; receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_full_expanded_covariance_pending_independent_readback'
    assert receipt['plan_sha256'] == sha(path) and receipt['scientific_eligibility'] is False and receipt['source_hashes'] == original
    names = ['stage_plan.json', 'family_components.tsv', 'patterns.tsv', 'species_contrast_design.npz', 'taxa.json',
        'case_covariance_index.tsv.gz', 'entity_incidence.tsv.gz', 'pattern_basis.npy', 'species_projection.npy'] + [t + '.npz' for t in plan['trees']]
    assert set(receipt['artifacts']) == set(names); bind(bindings, rp)
    for name, digest in receipt['artifacts'].items(): bind(bindings, root / name, digest)
    verify(bindings)
    assert json.loads((root / 'stage_plan.json').read_text()) == dict(plan_sha256=sha(path), schema='expanded-covariance-v1')
    taxa = source['taxa']; index = {t: i for i, t in enumerate(taxa)}
    assert json.loads((root / 'taxa.json').read_text()) == taxa
    patterns = {}; pattern_rows = {}
    with (root / 'patterns.tsv').open() as handle:
        rows = csv.DictReader(handle, delimiter='\t'); assert rows.fieldnames == PATTERN_FIELDS
        for position, row in enumerate(rows):
            labels, weights = json.loads(row['taxa']), json.loads(row['twice_weights'])
            assert len(labels) == len(weights) and len(labels) == len(set(labels))
            assert labels == sorted(labels) and all(t in index for t in labels)
            assert all(type(w) is int and w != 0 for w in weights) and sum(weights) == 0
            pattern = list(zip(labels, weights)); pid = identity('species-pattern', pattern)
            assert row['species_pattern_id'] == pid and int(row['row_index']) == position and pid not in patterns
            patterns[pid] = pattern; pattern_rows[pid] = position
    assert list(patterns) == sorted(patterns) and len(patterns) == plan['expected']['patterns']
    design = sparse.load_npz(root / 'species_contrast_design.npz'); dense = design.toarray()
    expected_design = np.zeros((len(patterns), len(taxa)))
    for pid, pairs in patterns.items():
        for taxon, weight in pairs: expected_design[pattern_rows[pid], index[taxon]] = weight / 2
    np.testing.assert_array_equal(dense, expected_design)
    seen_patterns = set(); family_links = defaultdict(set); entity_counts = Counter(); entity_owners = {}
    exported_cases = []; occurrences = 0; selection_uses = 0
    with gzip.open(root / 'case_covariance_index.tsv.gz', 'rt') as handle, gzip.open(root / 'entity_incidence.tsv.gz', 'rt') as incidence:
        rows = csv.DictReader(handle, delimiter='\t'); links = csv.DictReader(incidence, delimiter='\t')
        assert rows.fieldnames == CASE_FIELDS and links.fieldnames == INCIDENCE_FIELDS
        for case in source['cases']:
            observed = next(rows, None); assert observed is not None
            weights = {}
            for taxon, coefficient in [(case['focal_taxon'], 2), (case['background_taxon_a'], -1), (case['background_taxon_b'], -1)]:
                weights[taxon] = weights.get(taxon, 0) + coefficient
            wanted = sorted((t, w) for t, w in weights.items() if w); pid = identity('species-pattern', wanted)
            assert patterns[pid] == wanted; seen_patterns.add(pid)
            expected = {k: case[k] for k in CASE_FIELDS if k in case}
            expected.update(species_pattern_id=pid, species_pattern_row=str(pattern_rows[pid]),
                focal_index=str(index[case['focal_taxon']]), background_index_a=str(index[case['background_taxon_a']]),
                background_index_b=str(index[case['background_taxon_b']]))
            assert set(observed) == set(expected) | {'family_component'}
            assert all(observed[k] == v for k, v in expected.items()); exported_cases.append(observed)
            a, b = case['target_family'], case['background_family']; family_links[a].add(b); family_links[b].add(a)
            selection_uses += int(case['selection_records'])
            # Reconstruct the 14 original endpoint occurrences separately. No
            # gene/model collapsing or signed cancellation is permitted here.
            for side in ['target', 'background']:
                node = source['nodes'][side][case[side + '_id']]; sign = 1 if side == 'target' else -1
                expected_occurrences = [(side + '_node', [node['node_id']], '', sign, 1),
                    ('model_pair', [node['pair_key']], '', sign, 1), ('family', [node['family']], '', sign, 1)]
                for endpoint in ['a', 'b']:
                    expected_occurrences.append(('gene', [node['gene_' + endpoint]], endpoint, sign / 2, .5))
                    expected_occurrences.append(('model', [node['model_id_' + endpoint], int(node['version_' + endpoint]), node['sha256_' + endpoint]], endpoint, sign / 2, .5))
                for kind, token, endpoint, signed, unsigned in expected_occurrences:
                    actual = next(links, None); assert actual is not None; eid = identity(kind, token)
                    wanted_link = dict(case_id=case['case_id'], entity_kind=kind, entity_id=eid,
                        entity_value=json.dumps(token, separators=(',', ':')), side=side, endpoint=endpoint)
                    assert all(actual[k] == v for k, v in wanted_link.items())
                    assert float(actual['signed_loading']) == signed and float(actual['unsigned_loading']) == unsigned
                    entity_counts[kind, eid] += 1; occurrences += 1
                    key = kind, eid; family = node['family']
                    if key in entity_owners:
                        other = entity_owners[key]; family_links[family].add(other); family_links[other].add(family)
                    else: entity_owners[key] = family
        assert next(rows, None) is None and next(links, None) is None
    assert seen_patterns == set(patterns) and occurrences == 14 * len(exported_cases)
    families = sorted(family_links); family_rows = {f: i for i, f in enumerate(families)}
    edge_rows, edge_columns = [], []
    for f, neighbors in family_links.items():
        for other in neighbors: edge_rows.append(family_rows[f]); edge_columns.append(family_rows[other])
    adjacency = sparse.csr_matrix((np.ones(len(edge_rows)), (edge_rows, edge_columns)), shape=(len(families), len(families)))
    ncomponents, labels = connected_components(adjacency, directed=False)
    groups = defaultdict(list)
    for f, label in zip(families, labels): groups[int(label)].append(f)
    family_expected = {}
    for members in groups.values():
        component = identity('family-component', members)
        for family in members: family_expected[family] = dict(family=family, family_component=component, component_families=str(len(members)))
    seen_families = set()
    with (root / 'family_components.tsv').open() as handle:
        rows = csv.DictReader(handle, delimiter='\t'); assert rows.fieldnames == FAMILY_FIELDS
        for row in rows:
            assert row['family'] not in seen_families and row == family_expected[row['family']]; seen_families.add(row['family'])
    assert seen_families == set(families)
    for row in exported_cases:
        assert row['family_component'] == family_expected[row['target_family']]['family_component'] == family_expected[row['background_family']]['family_component']
    q, projection = np.load(root / 'pattern_basis.npy'), np.load(root / 'species_projection.npy')
    singular = np.linalg.svd(dense, compute_uv=False)
    tolerance = max(dense.shape) * np.finfo(float).eps * (singular[0] if len(singular) else 0)
    rank = int((singular > tolerance).sum()); assert q.shape == (len(patterns), rank) and projection.shape == (rank, len(taxa))
    assert np.isclose(receipt['svd_rank_tolerance'], tolerance, rtol=1e-12, atol=0)
    np.testing.assert_allclose(q.T @ q, np.eye(rank), atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(q @ projection, dense, atol=1e-12, rtol=1e-10)
    tree_checks = {}
    for label in plan['trees']:
        spec = source['kernel_plan']['trees'][label]
        tree = dendropy.Tree.get(path=spec['tree'], schema='newick', preserve_underscores=True)
        assert {leaf.taxon.label for leaf in tree.leaf_node_iter()} == set(taxa)
        edges = [node for node in tree.preorder_node_iter() if node.parent_node is not None]
        features = np.zeros((len(taxa), len(edges)))
        for column, node in enumerate(edges):
            length = float(node.edge_length); assert np.isfinite(length) and length >= 0
            for leaf in node.leaf_iter(): features[index[leaf.taxon.label], column] = np.sqrt(length)
        paths = design @ features
        with np.load(root / (label + '.npz')) as arrays:
            factor, core, lower = arrays['factor'], arrays['reduced_covariance'], arrays['cholesky']
            assert set(arrays.files) == {'factor', 'reduced_covariance', 'cholesky'}
        assert factor.shape == (len(patterns), rank) and lower.shape == core.shape == (rank, rank)
        np.testing.assert_array_equal(np.triu(lower, 1), 0); assert np.all(np.diag(lower) > 0)
        np.testing.assert_allclose(lower @ lower.T, core, atol=1e-11, rtol=1e-10)
        np.testing.assert_allclose(factor, q @ lower, atol=1e-12, rtol=1e-10)
        projected_paths = projection @ features
        np.testing.assert_allclose(projected_paths @ projected_paths.T, core, atol=1e-10, rtol=1e-10)
        minimum = float(np.linalg.eigvalsh(core)[0]) if rank else None
        claim = receipt['tree_checks'][label]; assert claim['covariance_entries_checked'] == len(patterns) ** 2
        if rank: assert np.isclose(minimum, claim['minimum_core_eigenvalue'], rtol=1e-8, atol=1e-12)
        else: assert claim['minimum_core_eigenvalue'] is None
        maximum = 0.; checked = 0
        for start in range(0, len(patterns), plan['covariance_block_rows']):
            end = min(start + plan['covariance_block_rows'], len(patterns))
            expected = paths[start:end] @ paths.T; actual = factor[start:end] @ factor.T
            np.testing.assert_allclose(actual, expected, atol=1e-11, rtol=1e-10)
            maximum = max(maximum, float(np.max(abs(actual - expected), initial=0))); checked += expected.size
        tree_checks[label] = dict(covariance_entries_checked=checked, maximum_tree_edge_error=maximum)
        print('independent_full_expanded_tree_covariance', label, checked, flush=True)
    summary = dict(logical_cases=len(exported_cases), physical_cases=len({r['physical_case_id'] for r in source['cases']}),
        selected_records=selection_uses, species_columns=len(taxa), patterns=len(patterns), zero_patterns=sum(not p for p in patterns.values()),
        rank=rank, families=len(families), family_components=ncomponents, entity_occurrences=occurrences,
        unique_entities=dict(Counter(k for k, e in entity_counts)), trees=plan['trees'], covariance_entries_per_tree=len(patterns) ** 2)
    assert all(summary[k] == receipt[k] for k in SUMMARY_FIELDS); verify(bindings)
    result = dict(status='passed_full_expanded_covariance_entity_and_tree_readback', plan_sha256=sha(path),
        producer_receipt_sha256=sha(rp), **summary, tree_checks=tree_checks, source_hashes=bindings,
        scientific_eligibility=False, scope=plan['scope'])
    with Path(output).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'scope']}, indent=2), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); run(a.plan, a.output)
