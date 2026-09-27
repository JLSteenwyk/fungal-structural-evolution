#!/usr/bin/env python3
"""Independently reconstruct normalized distances using DendroPy path matrices."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import dendropy

METRICS = {'ds_equal_alternative': 'dS_upstream_equal_alternative_convention',
           'dn_equal_alternative': 'dN_upstream_equal_alternative_convention'}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    with Path(p).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def check_rows(path, expected, keys):
    actual = read(path)
    indexed = {tuple(r[k] for k in keys): r for r in actual}
    assert len(indexed) == len(actual) == len(expected), path
    for row in expected:
        observed = indexed.pop(tuple(row[k] for k in keys))
        assert observed.keys() == row.keys(), path
        for key, value in row.items():
            if isinstance(value, (int, float)):
                assert math.isfinite(float(observed[key]))
                assert math.isclose(float(observed[key]), value, rel_tol=1e-10, abs_tol=1e-10), (path, key, value, observed[key])
            else:
                assert observed[key] == value, (path, key)
    assert not indexed


def load(spec):
    fits, norm, audit = (Path(spec[k]) for k in ('fits', 'normalization', 'audit'))
    fr, nr, ar = (json.loads((p / 'receipt.json').read_text()) for p in (fits, norm, audit))
    assert nr['status'] == 'complete_independently_checked_mg94_opportunity_normalization'
    assert ar['status'] == 'passed_saved_fit_readback' and ar['scope'] == 'all_executed_cases'
    assert nr['source_fit_receipt_sha256'] == sha(fits / 'receipt.json')
    assert nr['source_audit_sha256'] == sha(audit / 'receipt.json')
    expected = {r['case_id']: r['receipt_sha256'] for r in fr['case_receipts']}
    assert len(expected) == len(fr['case_receipts'])
    assert expected == {r['case_id']: r['source_receipt_sha256'] for r in ar['source_cases']}
    assert expected == {r['case_id']: r['source_fit_receipt_sha256'] for r in nr['checks']}
    for name, digest in nr['artifacts'].items():
        assert sha(norm / name) == digest
    grouped = defaultdict(dict)
    for r in read(norm / 'normalized_branches.tsv'):
        assert r['node'] not in grouped[r['case_id']]
        grouped[r['case_id']][r['node']] = r
    assert set(grouped) == set(expected)
    result = {}
    for case, digest in expected.items():
        folder = fits / case
        assert sha(folder / 'receipt.json') == digest
        receipt = json.loads((folder / 'receipt.json').read_text())
        config = json.loads((folder / 'config.json').read_text())
        assert receipt['config_sha256'] == sha(folder / 'config.json')
        assert config['tree_sha256'] == sha(folder / 'tree.nwk')
        tree = dendropy.Tree.get(path=str(folder / 'tree.nwk'), schema='newick', rooting='force-unrooted', preserve_underscores=True)
        taxa = sorted(n.taxon.label for n in tree.leaf_node_iter())
        assert len(taxa) == len(set(taxa))
        edges, used = {}, set()
        for node in tree.preorder_node_iter():
            if node is tree.seed_node:
                node.edge.length = 0.0
                continue
            label = node.taxon.label if node.taxon else node.label
            assert label not in used
            used.add(label)
            r = grouped[case][label]
            assert r['translation_table'] == config['translation_table']
            assert float(r['fitted_global_omega']) == receipt['reported_omega']
            a = sorted(n.taxon.label for n in node.leaf_iter())
            b = sorted(set(taxa).difference(a))
            split = tuple(a if (len(a), a) < (len(b), b) else b)
            assert split and split not in edges
            edges[split] = {m: float(r[c]) for m, c in METRICS.items()}
            assert all(math.isfinite(v) and v >= 0 for v in edges[split].values())
        assert used == set(grouped[case])
        paths = {}
        for metric, column in METRICS.items():
            for node in tree.preorder_node_iter():
                if node is not tree.seed_node:
                    label = node.taxon.label if node.taxon else node.label
                    node.edge.length = float(grouped[case][label][column])
            matrix = tree.phylogenetic_distance_matrix()
            paths[metric] = {(a, b): matrix.distance(tree.taxon_namespace.get_taxon(label=a), tree.taxon_namespace.get_taxon(label=b))
                             for a, b in itertools.combinations(taxa, 2)}
        result[case] = dict(taxa=taxa, edges=edges, paths=paths, code=config['translation_table'], omega=receipt['reported_omega'])
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--proof', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['plan_sha256'] == sha(args.plan)
    assert receipt['status'] == 'complete_normalized_codon_alignment_divergence_comparison_pending_readback'
    pins = dict(receipt['sources'])
    pins.update({str(root / name): digest for name, digest in receipt['artifacts'].items()})
    pins[str(args.plan)] = sha(args.plan)
    for p, digest in pins.items():
        assert sha(p) == digest
    old, new = load(plan['original']), load(plan['local'])
    ledger = read(plan['ledger'])
    assert len({r['case_id'] for r in ledger}) == len(ledger)
    assert set(old) | set(new) <= {r['case_id'] for r in ledger}
    cases, splits, pairs = [], [], []
    for flags in ledger:
        case = flags['case_id']
        a, b = old.get(case), new.get(case)
        status = 'matched' if a and b else 'original_only' if a else 'local_only' if b else 'neither_fitted'
        row = dict(case_id=case, comparison_status=status, original_global_omega=float(a['omega']) if a else '',
                   local_global_omega=float(b['omega']) if b else '', global_omega_delta=b['omega']-a['omega'] if a and b else '',
                   historical_review_flags=flags['case_specific_review_flags'], copy_caveat=flags['copy_caveat'], selection_eligibility='not_established')
        for m in METRICS:
            for label, fit in [('original', a), ('local', b)]:
                row[label+'_tree_'+m] = math.fsum(e[m] for e in fit['edges'].values()) if fit else ''
            row['tree_'+m+'_delta'] = row['local_tree_'+m]-row['original_tree_'+m] if a and b else ''
        cases.append(row)
        if not (a and b):
            continue
        assert a['taxa'] == b['taxa'] and a['code'] == b['code']
        for key in sorted(a['edges'].keys() | b['edges'].keys()):
            x, y = a['edges'].get(key), b['edges'].get(key)
            r = dict(case_id=case, smaller_side_taxa=';'.join(key), split_status='shared' if x and y else 'original_only' if x else 'local_only')
            for m in METRICS:
                r.update({'original_'+m: x[m] if x else '', 'local_'+m: y[m] if y else '', m+'_delta': y[m]-x[m] if x and y else ''})
            splits.append(r)
        for x, y in itertools.combinations(a['taxa'], 2):
            r = dict(case_id=case, taxon_a=x, taxon_b=y)
            for m in METRICS:
                left, right = a['paths'][m][x, y], b['paths'][m][x, y]
                r.update({'original_'+m: left, 'local_'+m: right, m+'_delta': right-left})
            pairs.append(r)
    check_rows(root/'cases.tsv', cases, ['case_id'])
    check_rows(root/'splits.tsv', splits, ['case_id', 'smaller_side_taxa'])
    check_rows(root/'pairs.tsv', pairs, ['case_id', 'taxon_a', 'taxon_b'])
    dispositions = dict(Counter(r['comparison_status'] for r in cases))
    for key, val in [('cases',len(cases)), ('split_rows',len(splits)), ('pair_rows',len(pairs)), ('dispositions',dispositions)]:
        assert receipt[key] == val
    for p, digest in pins.items():
        assert sha(p) == digest
    proof = dict(status='passed_full_normalized_codon_divergence_independent_readback', cases=len(cases), split_rows=len(splits), pair_rows=len(pairs), dispositions=dispositions,
                 source_receipt_sha256=sha(root/'receipt.json'), plan_sha256=sha(args.plan), script_sha256=sha(__file__), dendropy_version=dendropy.__version__,
                 scope='Independent DendroPy trees and direct path matrices check every normalized branch, taxon-pair distance, total, omega delta, flag, missing cell and ledger disposition; tolerance 1e-10 absolute/relative. Uses audited normalized lengths, does not independently refit models or establish selection eligibility.')
    args.proof.write_text(json.dumps(proof, indent=2)+'\n')
    print(json.dumps(proof, indent=2))


if __name__ == '__main__':
    main()
