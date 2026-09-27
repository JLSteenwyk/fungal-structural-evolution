#!/usr/bin/env python3
"""Reconstruct all tree-comparison rows using DendroPy and direct tree distances."""
import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import dendropy
import pandas as pd
import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path, keys):
    data = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    if data.duplicated(keys).any():
        raise ValueError('Duplicate output keys: ' + str(path))
    return data.set_index(keys)


def close(actual, expected):
    if expected == '':
        assert actual == '', (actual, expected)
    elif isinstance(expected, (int, float)):
        assert math.isfinite(float(actual)) and math.isclose(float(actual), expected, rel_tol=1e-10, abs_tol=1e-10), (actual, expected)
    else:
        assert actual == str(expected), (actual, expected)


def source(trees, audit):
    proof = json.loads((audit / 'receipt.json').read_text())
    assert proof['status'] == 'passed_full_genus_tree_audit' and not proof['pending_cases']
    case_hashes = {r['case_id']: r['receipt_sha256'] for r in proof['source_case_receipts']}
    for name, digest in proof['artifacts'].items():
        assert sha(audit / name) == digest
    summary = table(audit / 'case_audit.tsv', ['case_id'])
    assert set(summary.index) == set(case_hashes)
    parsed = {}
    for case, digest in case_hashes.items():
        folder = trees / case
        assert sha(folder / 'receipt.json') == digest
        receipt = json.loads((folder / 'receipt.json').read_text())
        assert sha(folder / 'tree.treefile') == receipt['artifacts']['tree.treefile']
        tree = dendropy.Tree.get(path=str(folder / 'tree.treefile'), schema='newick',
                                 preserve_underscores=True, rooting='force-unrooted')
        taxa = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
        edges = {}
        for node in tree.preorder_node_iter():
            if node.parent_node is None:
                continue
            descendants = sorted(leaf.taxon.label for leaf in node.leaf_iter())
            complement = sorted(taxa - set(descendants))
            selected = sorted([descendants, complement], key=lambda x: (len(x), x))[0]
            key = ';'.join(selected)
            assert key and key not in edges
            if len(selected) > 1:
                alrt, ufb = map(float, node.label.split('/'))
            else:
                alrt = ufb = ''
            edges[key] = dict(length=node.edge_length, sh_alrt=alrt, ufb=ufb)
        distances = tree.phylogenetic_distance_matrix()
        taxon_objects = {t.label: t for t in tree.taxon_namespace}
        pairs = {(a, b): distances.distance(taxon_objects[a], taxon_objects[b])
                 for a, b in itertools.combinations(sorted(taxa), 2)}
        parsed[case] = taxa, edges, pairs
    return summary, parsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison-plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--wait-for', type=Path)
    args = parser.parse_args()
    pinned = {str(p): sha(p) for p in [args.comparison_plan, Path(__file__)]}
    plan = json.loads(args.comparison_plan.read_text())
    if args.wait_for:
        pinned[str(args.wait_for)] = sha(args.wait_for)
        identity = json.loads(args.wait_for.read_text())
        while True:
            try:
                process = psutil.Process(identity['pid'])
                if process.create_time() != identity['created'] or process.status() == psutil.STATUS_ZOMBIE:
                    break
                if process.cmdline() != identity['cmdline']:
                    raise ValueError('Comparison controller changed identity')
            except psutil.NoSuchProcess:
                break
            print('waiting_for_comparison', identity['pid'], flush=True)
            time.sleep(30)
    paths = {k: Path(v) for k, v in plan['arguments'].items()}
    root = paths['output']
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_audited_codon_alignment_tree_comparison_pending_readback'
    for path, digest in receipt['sources'].items():
        assert sha(path) == digest
    for path, digest in receipt['artifacts'].items():
        assert sha(root / path) == digest
    ledger = table(paths['ledger'], ['case_id'])
    cases = table(root / 'cases.tsv', ['case_id'])
    splits = table(root / 'splits.tsv', ['case_id', 'smaller_side_taxa'])
    pairs = table(root / 'pairs.tsv', ['case_id', 'taxon_a', 'taxon_b'])
    assert set(cases.index) == set(ledger.index)
    old_summary, old = source(paths['original-trees'], paths['original-audit'])
    new_summary, new = source(paths['local-trees'], paths['local-audit'])
    assert (set(old) | set(new)) <= set(ledger.index)
    matched = set(old) & set(new)
    assert set(splits.index.get_level_values(0)) == set(pairs.index.get_level_values(0)) == matched
    changed = 0
    statuses = Counter()
    for case, row in cases.iterrows():
        status = 'matched' if case in matched else 'original_only' if case in old else 'local_only' if case in new else 'neither_fitted'
        assert row.comparison_status == status
        statuses[status] += 1
        for dst, src in [('historical_review_flags', 'case_specific_review_flags'), ('copy_caveat', 'copy_caveat'),
                         ('fcs_omission_status', 'fcs_omission_status'), ('alignment_median_pair_retention', 'alignment_median_pair_retention')]:
            assert row[dst] == ledger.loc[case, src]
        for label, summary in [('original', old_summary), ('local', new_summary)]:
            for suffix, src in [('nucleotide_columns', 'nucleotide_columns'), ('total_tree_length', 'total_tree_length')]:
                expected = float(summary.loc[case, src]) if case in summary.index else ''
                close(row[label + '_' + suffix], expected)
        metrics = ['taxa', 'shared_internal_splits', 'original_only_internal_splits', 'local_only_internal_splits',
                   'rf_distance', 'normalized_rf', 'total_tree_length_delta', 'mean_absolute_pair_distance_delta', 'maximum_absolute_pair_distance_delta']
        if case not in matched:
            assert all(row[k] == '' for k in metrics)
            continue
        taxa, old_edges, old_paths = old[case]
        new_taxa, new_edges, new_paths = new[case]
        assert taxa == new_taxa
        left, right = {s for s in old_edges if ';' in s}, {s for s in new_edges if ';' in s}
        rf = len(left - right) + len(right - left)
        changed += rf > 0
        absolute = [abs(new_paths[p] - v) for p, v in old_paths.items()]
        expected = [len(taxa), len(left & right), len(left - right), len(right - left), rf,
                    rf / (2 * (len(taxa) - 3)),
                    sum(v['length'] for v in new_edges.values()) - sum(v['length'] for v in old_edges.values()),
                    sum(absolute) / len(absolute), max(absolute)]
        for key, value in zip(metrics, expected):
            close(row[key], value)
        case_splits = splits.loc[case]
        assert set(case_splits.index) == set(old_edges) | set(new_edges)
        for split, observed in case_splits.iterrows():
            a, b = old_edges.get(split), new_edges.get(split)
            assert observed['internal'] == str(';' in split)
            assert observed['status'] == ('shared' if a and b else 'original_only' if a else 'local_only')
            for field in ['length', 'sh_alrt', 'ufb']:
                close(observed['original_' + field], a[field] if a else '')
                close(observed['local_' + field], b[field] if b else '')
                close(observed[field + '_delta'], b[field] - a[field] if a and b and a[field] != '' else '')
        case_pairs = pairs.loc[case]
        assert set(case_pairs.index) == set(old_paths)
        for pair, observed in case_pairs.iterrows():
            close(observed.original_distance, old_paths[pair])
            close(observed.local_distance, new_paths[pair])
            close(observed.distance_delta, new_paths[pair] - old_paths[pair])
    assert dict(statuses) == receipt['dispositions']
    assert len(cases) == receipt['cases'] and len(splits) == receipt['split_rows'] and len(pairs) == receipt['pair_rows']
    assert changed == receipt['matched_cases_with_changed_topology']
    for path, digest in pinned.items():
        assert sha(path) == digest
    result = dict(status='passed_full_codon_tree_comparison_independent_readback',
                  cases=len(cases), dispositions=dict(statuses), split_rows=len(splits), pair_rows=len(pairs),
                  changed_topology_cases=changed, source_receipt_sha256=sha(root / 'receipt.json'),
                  pins=pinned, dendropy_version=dendropy.__version__,
                  scope='Independent DendroPy parser and direct phylogenetic distance matrix reconstruct every split, support/length delta and taxon-pair distance. Exact case, flag, missing-field and full row grids checked; numerical tolerance 1e-10 absolute/relative. Does not rerun phylogenetic inference or establish model adequacy.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
