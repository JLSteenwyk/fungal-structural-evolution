#!/usr/bin/env python3
"""Count exact merged-family tree tasks without inferring or installing trees."""
import argparse
import csv
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def clusters(path):
    with path.open() as handle:
        if next(handle).strip() != '(mclmatrix' or next(handle).strip() != 'begin':
            raise ValueError('Unexpected cluster header')
        for number, line in enumerate(handle):
            if line.strip() == ')':
                if handle.read().strip():
                    raise ValueError('Content after closing delimiter')
                return
            fields = line.split()
            if fields[0] != str(number) or fields[-1] != '$':
                raise ValueError('Unexpected family index or terminator')
            genes = fields[1:-1]
            if not genes or len(set(genes)) != len(genes):
                raise ValueError('Empty family or repeated gene')
            yield f'OG{number:07d}', genes
    raise ValueError('Missing closing delimiter')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--merged', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    receipt = json.loads((args.merged / 'receipt.json').read_text())
    audit = json.loads((args.merged / 'readback.json').read_text())
    if audit['status'] != 'passed_complete_guide_discovery_merge_readback':
        raise ValueError('Independent membership readback required')
    if audit['merge_receipt_sha256'] != sha(args.merged / 'receipt.json'):
        raise ValueError('Changed merge receipt')
    summaries = []
    tasks = {}
    for guide in receipt['guides']:
        for relative, digest in guide['artifacts'].items():
            if sha(args.merged / relative) != digest:
                raise ValueError('Changed merged artifact: ' + relative)
        counts = Counter()
        folder = args.merged / guide['guide']
        with (folder / 'family_sources.tsv').open() as handle:
            rows = csv.DictReader(handle, delimiter='\t')
            for row, pair in itertools.zip_longest(rows, clusters(folder / 'clusters_id_pairs.txt')):
                if row is None or pair is None:
                    raise ValueError('Crosswalk/cluster length differs')
                family, genes = pair
                digest = hashlib.sha256(('\n'.join(sorted(genes)) + '\n').encode()).hexdigest()
                if (family != row['new_family'] or len(genes) != int(row['proteins'])
                        or digest != row['membership_sha256']):
                    raise ValueError('Crosswalk membership differs')
                kind = row['source_type']
                if kind not in ('retained', 'discovery'):
                    raise ValueError('Unknown source kind: ' + kind)
                size = len(genes)
                taxa = len({gene.split('_')[0] for gene in genes})
                bucket = 'singleton' if size == 1 else 'pair' if size == 2 else 'tree'
                counts['families'] += 1
                counts['proteins'] += size
                counts[kind + '_' + bucket + '_families'] += 1
                counts[kind + '_' + bucket + '_proteins'] += size
                counts['single_taxon_families'] += taxa == 1
                if size < 3:
                    if row['original_tree_candidate']:
                        raise ValueError('Small family incorrectly marked as tree candidate')
                    continue
                if (kind == 'retained') != bool(row['original_tree_candidate']):
                    raise ValueError('Tree candidate/source mismatch')
                task = dict(membership_sha256=digest, source_type=kind,
                            source_clade=row['source_clade'], source_family=row['source_family'],
                            proteins=size, taxa=taxa, representative_guide=guide['guide'],
                            representative_family=family,
                            original_tree_candidate=row['original_tree_candidate'])
                if digest in tasks:
                    previous = tasks[digest]
                    for field in ['source_type', 'proteins', 'taxa', 'original_tree_candidate']:
                        if previous[field] != task[field]:
                            raise ValueError('Inconsistent identical-membership task')
                    previous['guide_occurrences'] += 1
                else:
                    task['guide_occurrences'] = 1
                    tasks[digest] = task
        if counts['families'] != guide['families'] or counts['proteins'] != guide['proteins']:
            raise ValueError('Guide totals differ')
        summaries.append(dict(guide=guide['guide'], **counts))
    args.output.mkdir(parents=True)
    table = args.output / 'unique_tree_tasks.tsv'
    ordered = sorted(tasks.values(), key=lambda r: (-r['proteins'], r['membership_sha256']))
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ordered[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(ordered)
    result = dict(status='complete_expanded_family_tree_workload_census', guides=summaries,
                  unique_tree_tasks=len(tasks),
                  retained_tree_candidates=sum(t['source_type'] == 'retained' for t in ordered),
                  unique_new_tree_tasks=sum(t['source_type'] == 'discovery' for t in ordered),
                  tasks_shared_across_guides=sum(t['guide_occurrences'] == 2 for t in ordered),
                  largest_new_family=max(t['proteins'] for t in ordered if t['source_type'] == 'discovery'),
                  merge_receipt_sha256=sha(args.merged / 'receipt.json'),
                  merge_readback_sha256=sha(args.merged / 'readback.json'),
                  script_sha256=sha(Path(__file__)), artifacts={table.name: sha(table)},
                  interpretation='Exact membership deduplication across guide partitions. All families retained in source partitions; singleton/pair families require no inferred bifurcating gene tree. Retained trees are candidates requiring validation before reuse. New alignments, trees and reconciliation remain pending. Counts are not runtime estimates or biological validation.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
