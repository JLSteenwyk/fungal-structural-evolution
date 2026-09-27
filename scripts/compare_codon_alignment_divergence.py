#!/usr/bin/env python3
"""Compare fully audited normalized codon fits by unrooted split and taxon pair."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

from Bio import Phylo
import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def load(spec):
    fits, norm, audit = [Path(spec[k]) for k in ['fits', 'normalization', 'audit']]
    fit_receipt, receipt, proof = [json.loads((p / 'receipt.json').read_text()) for p in [fits, norm, audit]]
    if (receipt['status'] != 'complete_independently_checked_mg94_opportunity_normalization'
            or receipt['source_fit_receipt_sha256'] != sha(fits / 'receipt.json')
            or receipt['source_audit_sha256'] != sha(audit / 'receipt.json')
            or proof['status'] != 'passed_saved_fit_readback' or proof['scope'] != 'all_executed_cases'):
        raise ValueError('Full verified normalization and clean numerical audit required')
    expected = {r['case_id']: r['receipt_sha256'] for r in fit_receipt['case_receipts']}
    if expected != {r['case_id']: r['source_receipt_sha256'] for r in proof['source_cases']}:
        raise ValueError('Fit/audit grids differ')
    if expected != {r['case_id']: r['source_fit_receipt_sha256'] for r in receipt['checks']}:
        raise ValueError('Normalization fit grid differs')
    for name, digest in receipt['artifacts'].items():
        if sha(norm / name) != digest:
            raise ValueError('Normalization artifact changed')
    grouped = defaultdict(dict)
    for r in table(norm / 'normalized_branches.tsv'):
        if r['node'] in grouped[r['case_id']]:
            raise ValueError('Duplicate normalized node')
        grouped[r['case_id']][r['node']] = r
    if set(grouped) != set(expected):
        raise ValueError('Normalized branch cases differ')
    out = {}
    for case, digest in expected.items():
        folder = fits / case
        if sha(folder / 'receipt.json') != digest:
            raise ValueError('Fit receipt changed')
        r = json.loads((folder / 'receipt.json').read_text())
        config = json.loads((folder / 'config.json').read_text())
        if r['config_sha256'] != sha(folder / 'config.json') or config['tree_sha256'] != sha(folder / 'tree.nwk'):
            raise ValueError('Fit tree binding changed')
        tree = Phylo.read(folder / 'tree.nwk', 'newick')
        names = [t.name for t in tree.get_terminals()]
        taxa = set(names)
        nodes = [n for n in tree.find_clades() if n is not tree.root]
        if len(taxa) != len(names) or set(grouped[case]) != {n.name for n in nodes} or len(nodes) != len(grouped[case]):
            raise ValueError('Node grid differs')
        edges = {}
        for node in nodes:
            side = {t.name for t in node.get_terminals()}
            key = min(tuple(sorted(side)), tuple(sorted(taxa - side)), key=lambda s: (len(s), s))
            if not key or key in edges:
                raise ValueError('Ambiguous split identity')
            row = grouped[case][node.name]
            if row['translation_table'] != config['translation_table'] or float(row['fitted_global_omega']) != r['reported_omega']:
                raise ValueError('Code or omega mismatch')
            values = {name: float(row[col]) for name, col in
                      [('ds_equal_alternative', 'dS_upstream_equal_alternative_convention'),
                       ('dn_equal_alternative', 'dN_upstream_equal_alternative_convention')]}
            if any(not math.isfinite(v) or v < 0 for v in values.values()):
                raise ValueError('Invalid normalized length')
            edges[key] = values
        out[case] = dict(taxa=taxa, edges=edges, omega=r['reported_omega'], code=config['translation_table'])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--wait-for', type=Path)
    args = parser.parse_args()
    plan_hash, script_hash = sha(args.plan), sha(__file__)
    plan = json.loads(args.plan.read_text())
    if args.wait_for:
        identity = json.loads(args.wait_for.read_text())
        while True:
            try:
                process = psutil.Process(identity['pid'])
                if process.create_time() != identity['created'] or process.status() == psutil.STATUS_ZOMBIE:
                    break
                if process.cmdline() != identity['cmdline']:
                    raise ValueError('Normalization controller identity changed')
            except psutil.NoSuchProcess:
                break
            print('waiting_for_complete_normalization', identity['pid'], flush=True)
            time.sleep(30)
    sources = [Path(plan[label][k]) / 'receipt.json' for label in ['original', 'local'] for k in ['fits', 'normalization', 'audit']]
    sources.append(Path(plan['ledger']))
    hashes = {str(p): sha(p) for p in sources}
    old, new = load(plan['original']), load(plan['local'])
    ledger_rows = table(plan['ledger'])
    ledger = {r['case_id']: r for r in ledger_rows}
    if len(ledger) != len(ledger_rows) or (set(old) | set(new)) - set(ledger):
        raise ValueError('Full ledger grid differs')
    case_rows, edge_rows, pair_rows = [], [], []
    metrics = ['ds_equal_alternative', 'dn_equal_alternative']
    for case, flags in ledger.items():
        a, b = old.get(case), new.get(case)
        status = 'matched' if a and b else 'original_only' if a else 'local_only' if b else 'neither_fitted'
        row = dict(case_id=case, comparison_status=status, original_global_omega=a['omega'] if a else '',
                   local_global_omega=b['omega'] if b else '', global_omega_delta=b['omega'] - a['omega'] if a and b else '',
                   historical_review_flags=flags['case_specific_review_flags'], copy_caveat=flags['copy_caveat'],
                   selection_eligibility='not_established')
        if a and b:
            if a['taxa'] != b['taxa'] or a['code'] != b['code']:
                raise ValueError('Cannot compare changed taxa or genetic code')
            for split in sorted(set(a['edges']) | set(b['edges'])):
                left, right = a['edges'].get(split), b['edges'].get(split)
                edge = dict(case_id=case, smaller_side_taxa=';'.join(split),
                            split_status='shared' if left and right else 'original_only' if left else 'local_only')
                for metric in metrics:
                    edge['original_' + metric] = left[metric] if left else ''
                    edge['local_' + metric] = right[metric] if right else ''
                    edge[metric + '_delta'] = right[metric] - left[metric] if left and right else ''
                edge_rows.append(edge)
            for x, y in itertools.combinations(sorted(a['taxa']), 2):
                pair = dict(case_id=case, taxon_a=x, taxon_b=y)
                for metric in metrics:
                    left = sum(v[metric] for split, v in a['edges'].items() if (x in split) != (y in split))
                    right = sum(v[metric] for split, v in b['edges'].items() if (x in split) != (y in split))
                    pair.update({'original_' + metric: left, 'local_' + metric: right, metric + '_delta': right-left})
                pair_rows.append(pair)
        for metric in metrics:
            for label, fit in [('original', a), ('local', b)]:
                row[label + '_tree_' + metric] = sum(v[metric] for v in fit['edges'].values()) if fit else ''
            row['tree_' + metric + '_delta'] = row['local_tree_' + metric] - row['original_tree_' + metric] if a and b else ''
        case_rows.append(row)
    for p in sources:
        if sha(p) != hashes[str(p)]:
            raise ValueError('Source changed during comparison')
    if sha(args.plan) != plan_hash or sha(__file__) != script_hash:
        raise ValueError('Plan/script changed')
    output = Path(plan['output'])
    output.mkdir(parents=True, exist_ok=False)
    for name, rows in [('cases.tsv', case_rows), ('splits.tsv', edge_rows), ('pairs.tsv', pair_rows)]:
        with (output / name).open('x') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    receipt = dict(status='complete_normalized_codon_alignment_divergence_comparison_pending_readback',
                   cases=len(case_rows), dispositions=dict(Counter(r['comparison_status'] for r in case_rows)),
                   split_rows=len(edge_rows), pair_rows=len(pair_rows), sources=hashes,
                   plan_sha256=plan_hash, script_sha256=script_hash,
                   artifacts={p.name: sha(p) for p in output.iterdir()},
                   scope='Same-taxon/code normalized-distance sensitivity. Branches matched by unrooted splits, never internal labels. Missing branches are absent, not zero. Original means archived baseline fits, not later multistart optima; historical identifiability flags retained. Pair distances are correlated, not independent samples. No cross-alignment likelihood test, branch omega inference, saturation clearance or selection claim.')
    (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['cases', 'dispositions', 'split_rows', 'pair_rows']}), flush=True)


if __name__ == '__main__':
    main()
