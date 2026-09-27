#!/usr/bin/env python3
"""Compare audited unrooted nucleotide trees on exactly matched case/taxon grids."""
import argparse
from collections import Counter
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

from Bio import Phylo


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def keyed(rows):
    result = {r['case_id']: r for r in rows}
    if len(result) != len(rows):
        raise ValueError('Duplicate cases')
    return result


def tree_map(tree):
    tips = [n.name for n in tree.get_terminals()]
    taxa = set(tips)
    if len(taxa) != len(tips) or len(taxa) < 4:
        raise ValueError('Invalid tree tip grid')
    edges = {}
    for node in tree.find_clades():
        if node is tree.root:
            continue
        length = node.branch_length
        if length is None or not math.isfinite(length) or length < 0:
            raise ValueError('Invalid branch length')
        side = {n.name for n in node.get_terminals()}
        key = min(tuple(sorted(side)), tuple(sorted(taxa - side)), key=lambda s: (len(s), s))
        if not key or key in edges:
            raise ValueError('Empty/duplicate unrooted split; unsupported degree-two root')
        if len(key) > 1:
            supports = str(node.name).split('/')
            if len(supports) != 2:
                raise ValueError('Missing internal support labels')
            alrt, ufb = map(float, supports)
            if any(not math.isfinite(s) or not 0 <= s <= 100 for s in [alrt, ufb]):
                raise ValueError('Invalid support')
        else:
            alrt = ufb = ''
        edges[key] = dict(length=length, sh_alrt=alrt, ufb=ufb)
    if len(edges) != 2 * len(taxa) - 3 or sum(len(s) > 1 for s in edges) != len(taxa) - 3:
        raise ValueError('Expected resolved unrooted tree')
    return taxa, edges


def pair_distances(taxa, edges):
    return {(a, b): sum(v['length'] for split, v in edges.items() if (a in split) != (b in split))
            for a, b in itertools.combinations(sorted(taxa), 2)}


def load_source(trees, audit):
    receipt = json.loads((audit / 'receipt.json').read_text())
    batch = json.loads((trees / 'receipt.json').read_text())
    if receipt['status'] != 'passed_full_genus_tree_audit' or receipt['pending_cases']:
        raise ValueError('Full completed audit required')
    if receipt['source_config_sha256'] != sha(trees / 'config.json') or batch['config_sha256'] != receipt['source_config_sha256']:
        raise ValueError('Audit/config binding differs')
    source = {r['case_id']: r['receipt_sha256'] for r in receipt['source_case_receipts']}
    if source != {r['case_id']: r['receipt_sha256'] for r in batch['case_receipts']}:
        raise ValueError('Audit/producer case grids differ')
    for name, expected in receipt['artifacts'].items():
        if sha(audit / name) != expected:
            raise ValueError('Changed audited artifact: ' + name)
    summaries = keyed(read_table(audit / 'case_audit.tsv'))
    if set(summaries) != set(source) or len(source) != receipt['audited_cases']:
        raise ValueError('Audit summary grid differs')
    parsed = {}
    for case, expected in source.items():
        folder = trees / case
        if sha(folder / 'receipt.json') != expected:
            raise ValueError('Changed case receipt: ' + case)
        case_receipt = json.loads((folder / 'receipt.json').read_text())
        if sha(folder / 'tree.treefile') != case_receipt['artifacts']['tree.treefile']:
            raise ValueError('Changed audited tree: ' + case)
        tree = Phylo.read(folder / 'tree.treefile', 'newick')
        taxa, edges = tree_map(tree)
        if len(taxa) != int(summaries[case]['taxa']):
            raise ValueError('Tree/audit taxon count differs')
        parsed[case] = taxa, edges
    return summaries, parsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['original-trees', 'original-audit', 'local-trees', 'local-audit', 'ledger', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    ledger = keyed(read_table(args.ledger))
    sources = [args.original_trees / 'receipt.json', args.original_audit / 'receipt.json',
               args.local_trees / 'receipt.json', args.local_audit / 'receipt.json', args.ledger]
    hashes = {str(p): sha(p) for p in sources}
    before_summary, before = load_source(args.original_trees, args.original_audit)
    after_summary, after = load_source(args.local_trees, args.local_audit)
    if (set(before) | set(after)) - set(ledger):
        raise ValueError('Trees outside full ledger')
    cases, splits, pairs = [], [], []
    for case, flags in ledger.items():
        matched = case in before and case in after
        status = 'matched' if matched else 'original_only' if case in before else 'local_only' if case in after else 'neither_fitted'
        row = dict(case_id=case, comparison_status=status, historical_review_flags=flags['case_specific_review_flags'],
                   copy_caveat=flags['copy_caveat'], fcs_omission_status=flags['fcs_omission_status'],
                   alignment_median_pair_retention=flags['alignment_median_pair_retention'],
                   taxa='', original_nucleotide_columns='', local_nucleotide_columns='',
                   original_total_tree_length='', local_total_tree_length='', total_tree_length_delta='',
                   shared_internal_splits='', original_only_internal_splits='', local_only_internal_splits='',
                   rf_distance='', normalized_rf='', mean_absolute_pair_distance_delta='', maximum_absolute_pair_distance_delta='')
        for label, summary in [('original', before_summary), ('local', after_summary)]:
            if case in summary:
                row[label + '_nucleotide_columns'] = summary[case]['nucleotide_columns']
                row[label + '_total_tree_length'] = summary[case]['total_tree_length']
        if matched:
            taxa, old_edges = before[case]
            new_taxa, new_edges = after[case]
            if taxa != new_taxa:
                raise ValueError('Exact taxon grid differs: ' + case)
            old_internal = {s for s in old_edges if len(s) > 1}
            new_internal = {s for s in new_edges if len(s) > 1}
            rf = len(old_internal ^ new_internal)
            old_paths, new_paths = pair_distances(taxa, old_edges), pair_distances(taxa, new_edges)
            delta = [abs(new_paths[k] - v) for k, v in old_paths.items()]
            row.update(taxa=len(taxa), shared_internal_splits=len(old_internal & new_internal),
                       original_only_internal_splits=len(old_internal - new_internal),
                       local_only_internal_splits=len(new_internal - old_internal), rf_distance=rf,
                       normalized_rf=rf / (2 * (len(taxa) - 3)),
                       total_tree_length_delta=sum(v['length'] for v in new_edges.values()) - sum(v['length'] for v in old_edges.values()),
                       mean_absolute_pair_distance_delta=sum(delta) / len(delta), maximum_absolute_pair_distance_delta=max(delta))
            for split in sorted(set(old_edges) | set(new_edges)):
                old, new = old_edges.get(split), new_edges.get(split)
                entry = dict(case_id=case, smaller_side_taxa=';'.join(split), internal=len(split) > 1,
                             status='shared' if old and new else 'original_only' if old else 'local_only')
                for key in ['length', 'sh_alrt', 'ufb']:
                    entry['original_' + key] = old[key] if old else ''
                    entry['local_' + key] = new[key] if new else ''
                    entry[key + '_delta'] = new[key] - old[key] if old and new and old[key] != '' else ''
                splits.append(entry)
            for pair, distance in old_paths.items():
                pairs.append(dict(case_id=case, taxon_a=pair[0], taxon_b=pair[1], original_distance=distance,
                                  local_distance=new_paths[pair], distance_delta=new_paths[pair] - distance))
        cases.append(row)
    for p in sources:
        if sha(p) != hashes[str(p)]:
            raise ValueError('Source changed during comparison')
    args.output.mkdir(parents=True, exist_ok=False)
    for name, rows in [('cases.tsv', cases), ('splits.tsv', splits), ('pairs.tsv', pairs)]:
        if not rows:
            raise ValueError('Comparison requires matched cases')
        with (args.output / name).open('x') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    receipt = dict(status='complete_audited_codon_alignment_tree_comparison_pending_readback',
                   cases=len(cases), dispositions=dict(Counter(r['comparison_status'] for r in cases)),
                   matched_cases_with_changed_topology=sum(r['rf_distance'] != '' and r['rf_distance'] > 0 for r in cases),
                   split_rows=len(splits), pair_rows=len(pairs), sources=hashes, script_sha256=sha(__file__),
                   artifacts={p.name: sha(p) for p in args.output.iterdir()},
                   scope='Descriptive same-taxon unrooted topology and nucleotide-length sensitivity. Zero-length edges retained as inferred resolutions. Missing splits are absent, not zero-length. Likelihoods across different alignment columns are not compared. Pair distances and cases are not independent observations. No selection or homology-correctness claim.')
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['cases', 'dispositions', 'matched_cases_with_changed_topology', 'pair_rows']}))


if __name__ == '__main__':
    main()
