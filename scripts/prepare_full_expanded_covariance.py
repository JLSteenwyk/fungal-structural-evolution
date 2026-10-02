#!/usr/bin/env python3
"""Preserve complete entity reuse and factor all expanded species contrasts."""
import argparse
from collections import Counter, defaultdict
import csv
import fcntl
import gzip
import json
from pathlib import Path
import shutil
import numpy as np
from scipy import sparse
from full_expanded_covariance_sources import load, identity, CASE_FIELDS, INCIDENCE_FIELDS, PATTERN_FIELDS, FAMILY_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def occurrences(case, nodes):
    for side, sign in [('target', 1), ('background', -1)]:
        node = nodes[side][case[side + '_id']]
        for kind, token in [(side + '_node', [node['node_id']]),
            ('model_pair', [node['pair_key']]), ('family', [node['family']])]:
            yield kind, token, side, '', sign, 1
        for end in ['a', 'b']:
            yield 'gene', [node['gene_' + end]], side, end, sign / 2, .5
            yield 'model', [node['model_id_' + end], int(node['version_' + end]), node['sha256_' + end]], side, end, sign / 2, .5


def run(path):
    path = Path(path); plan = json.loads(path.read_text()); source, bindings = load(plan, path)
    root = Path(plan['output']); assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    root.mkdir(exist_ok=True); lock = (root / 'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (root / 'receipt.json').exists(), 'Completed stage cannot be restarted'
    state = dict(plan_sha256=sha(path), schema='expanded-covariance-v1')
    marker = root / 'stage_plan.json'
    if marker.exists(): assert json.loads(marker.read_text()) == state
    else: marker.write_text(json.dumps(state, indent=2) + '\n')
    cases, nodes, taxa = source['cases'], source['nodes'], source['taxa']; taxa_index = {t: i for i, t in enumerate(taxa)}
    parent = {}; entity_family = {}; patterns = {}; assignments = {}; entity_counts = Counter()
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        a, b = find(a), find(b)
        if a != b: parent[max(a, b)] = min(a, b)
    incidence_tmp = root / 'entity_incidence.tsv.gz.tmp'
    with gzip.open(incidence_tmp, 'wt', compresslevel=1) as handle:
        writer = csv.DictWriter(handle, INCIDENCE_FIELDS, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for case in cases:
            family = case['target_family']; union(family, case['background_family'])
            for kind, token, side, end, signed, unsigned in occurrences(case, nodes):
                eid = identity(kind, token); entity_counts[kind, eid] += 1
                if (kind, eid) in entity_family: union(family, entity_family[kind, eid])
                else: entity_family[kind, eid] = family
                writer.writerow(dict(case_id=case['case_id'], entity_kind=kind, entity_id=eid,
                    entity_value=json.dumps(token, separators=(',', ':')), side=side, endpoint=end,
                    signed_loading=signed, unsigned_loading=unsigned))
            weights = Counter()
            for field, coefficient in [('focal_taxon', 2), ('background_taxon_a', -1), ('background_taxon_b', -1)]:
                weights[case[field]] += coefficient
            pattern = tuple(sorted((t, w) for t, w in weights.items() if w)); assert sum(w for t, w in pattern) == 0
            pid = identity('species-pattern', pattern); assert pid not in patterns or patterns[pid] == pattern
            patterns[pid] = pattern; assignments[case['case_id']] = pid
    components = defaultdict(list)
    for family in sorted(parent): components[find(family)].append(family)
    family_index = {}
    with (root / 'family_components.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, FAMILY_FIELDS, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for members in sorted(components.values()):
            component = identity('family-component', members)
            for family in members:
                family_index[family] = component
                writer.writerow(dict(family=family, family_component=component, component_families=len(members)))
    labels = sorted(patterns); row_index = {pid: i for i, pid in enumerate(labels)}
    rows, columns, values = [], [], []
    with (root / 'patterns.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, PATTERN_FIELDS, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for i, pid in enumerate(labels):
            pattern = patterns[pid]
            writer.writerow(dict(species_pattern_id=pid, row_index=i,
                taxa=json.dumps([t for t, w in pattern], separators=(',', ':')),
                twice_weights=json.dumps([w for t, w in pattern], separators=(',', ':'))))
            for t, w in pattern: rows.append(i); columns.append(taxa_index[t]); values.append(w / 2)
    design = sparse.csr_matrix((values, (rows, columns)), shape=(len(labels), len(taxa)))
    assert np.all(np.asarray(design.sum(axis=1)) == 0)
    sparse.save_npz(root / 'species_contrast_design.npz', design)
    (root / 'taxa.json').write_text(json.dumps(taxa, indent=2) + '\n')
    with gzip.open(root / 'case_covariance_index.tsv.gz.tmp', 'wt', compresslevel=1) as handle:
        writer = csv.DictWriter(handle, CASE_FIELDS, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for case in cases:
            result = {k: case[k] for k in CASE_FIELDS if k in case}; pid = assignments[case['case_id']]
            assert family_index[case['target_family']] == family_index[case['background_family']]
            result.update(family_component=family_index[case['target_family']], species_pattern_id=pid,
                species_pattern_row=row_index[pid], focal_index=taxa_index[case['focal_taxon']],
                background_index_a=taxa_index[case['background_taxon_a']], background_index_b=taxa_index[case['background_taxon_b']])
            writer.writerow(result)
    dense = design.toarray(); u, singular, vt = np.linalg.svd(dense, full_matrices=False)
    tolerance = max(dense.shape) * np.finfo(float).eps * (singular[0] if len(singular) else 0)
    rank = int((singular > tolerance).sum()); q = u[:, :rank]; projection = singular[:rank, None] * vt[:rank]
    np.testing.assert_allclose(q @ projection, dense, atol=1e-12, rtol=1e-10)
    np.save(root / 'pattern_basis.npy', q); np.save(root / 'species_projection.npy', projection)
    tree_summaries = {}
    for tree in plan['trees']:
        with np.load(source['kernel'] / (tree + '.npz')) as arrays: kernel = arrays['centered_kernel']
        core = projection @ kernel @ projection.T
        assert np.max(abs(core - core.T), initial=0) < 1e-10; core = (core + core.T) / 2
        lower = np.linalg.cholesky(core); factor = q @ lower
        np.savez_compressed(root / (tree + '.npz'), factor=factor, reduced_covariance=core, cholesky=lower)
        maximum = 0.; checked = 0
        for start in range(0, len(labels), plan['covariance_block_rows']):
            end = min(start + plan['covariance_block_rows'], len(labels))
            expected = (design[start:end] @ kernel) @ design.T; actual = factor[start:end] @ factor.T
            np.testing.assert_allclose(actual, expected, atol=1e-11, rtol=1e-10)
            maximum = max(maximum, float(np.max(abs(actual - expected), initial=0))); checked += expected.size
        tree_summaries[tree] = dict(covariance_entries_checked=checked, maximum_absolute_error=maximum,
            minimum_core_eigenvalue=float(np.linalg.eigvalsh(core)[0]) if rank else None)
        print('full_expanded_covariance_tree', tree, 'patterns', len(labels), 'rank', rank, flush=True)
    verify(bindings)
    incidence_tmp.replace(root / 'entity_incidence.tsv.gz')
    (root / 'case_covariance_index.tsv.gz.tmp').replace(root / 'case_covariance_index.tsv.gz')
    summary = dict(logical_cases=len(cases), physical_cases=len({r['physical_case_id'] for r in cases}),
        selected_records=sum(int(r['selection_records']) for r in cases), species_columns=len(taxa), patterns=len(labels),
        zero_patterns=sum(not p for p in patterns.values()), rank=rank, families=len(parent), family_components=len(components),
        entity_occurrences=sum(entity_counts.values()), unique_entities=dict(Counter(k for k, e in entity_counts)),
        trees=plan['trees'], covariance_entries_per_tree=len(labels) ** 2)
    assert summary['patterns'] == plan['expected']['patterns'] and summary['entity_occurrences'] == 14 * len(cases)
    names = ['stage_plan.json', 'family_components.tsv', 'patterns.tsv', 'species_contrast_design.npz', 'taxa.json',
        'case_covariance_index.tsv.gz', 'entity_incidence.tsv.gz', 'pattern_basis.npy', 'species_projection.npy'] + [t + '.npz' for t in plan['trees']]
    receipt = dict(status='complete_full_expanded_covariance_pending_independent_readback', plan_sha256=sha(path),
        **summary, svd_rank_tolerance=float(tolerance), tree_checks=tree_summaries, source_hashes=bindings,
        artifacts={n: sha(root / n) for n in names}, scientific_eligibility=False, scope=plan['scope'])
    with (root / 'receipt.json').open('x') as handle: handle.write(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k not in ['source_hashes', 'scope']}, indent=2), flush=True)
    return receipt


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); run(p.parse_args().plan)
