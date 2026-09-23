#!/usr/bin/env python3
"""Stratify validated within-partition boundary disagreement by added residues."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--audit', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    audit = json.loads(args.audit.read_text())
    if audit['status'] != 'passed_full_domain_boundary_cluster_output_readback':
        raise ValueError('Full boundary readback required')
    pins = audit['source_hashes']
    for path, digest in pins.items():
        if sha(path) != digest:
            raise ValueError('Changed validated source: ' + path)
    def source(suffix):
        matches = [Path(p) for p in pins if p.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError('Ambiguous source: ' + suffix)
        return matches[0]
    intervals = {}
    with source('/intervals.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['interval_id']
            if key in intervals:
                raise ValueError('Repeated interval')
            start, end = int(row['start']), int(row['end'])
            if end - start + 1 != int(row['residues']):
                raise ValueError('Interval length differs')
            intervals[key] = (row['model_key'], start, end)
    labels = ['0', '1-4', '5-9', '10-19', '20-49', '50+']
    counts = Counter(); dispositions = Counter(); seen = set()
    with gzip.open(source('/boundary_cluster_dispositions.tsv.gz'), 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['model_key'], row['hit_id']
            if key in seen:
                raise ValueError('Repeated boundary pair')
            seen.add(key)
            am, a, b = intervals[row['alignment_interval']]
            em, c, d = intervals[row['envelope_interval']]
            if am != em or am != key[0] or not c <= a <= b <= d:
                raise ValueError('Envelope does not contain alignment')
            added = a - c + d - b
            label = labels[sum(added >= cut for cut in [1, 5, 10, 20, 50])]
            same = row['alignment_cluster'] == row['envelope_cluster']
            expected = ('identical_interval' if added == 0 else
                        'distinct_intervals_same_cluster' if same else
                        'distinct_intervals_different_clusters')
            if row['disposition'] != expected or (added == 0 and not same):
                raise ValueError('Disposition differs')
            counts[label, 'pairs'] += 1
            counts[label, 'different_clusters'] += int(not same)
            dispositions[expected] += 1
    if len(intervals) != audit['intervals'] or len(seen) != audit['candidate_pairs'] or dict(dispositions) != audit['dispositions']:
        raise ValueError('Full input scope differs')
    args.output.mkdir(parents=True, exist_ok=False)
    table = args.output / 'extension_cluster_disagreement.tsv'
    with table.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['added_residues', 'pairs', 'different_clusters', 'fraction_different'])
        for label in labels:
            n, k = counts[label, 'pairs'], counts[label, 'different_clusters']
            writer.writerow([label, n, k, format(k / n, '.17g') if n else 'NA'])
    for path, digest in pins.items():
        if sha(path) != digest:
            raise ValueError('Source changed during summary')
    result = dict(status='complete_boundary_extension_cluster_summary',
                  candidate_pairs=len(seen), intervals=len(intervals), bins=len(labels),
                  audit_sha256=sha(args.audit), script_sha256=sha(__file__),
                  artifacts={table.name: sha(table)},
                  scope='Descriptive full paired-boundary census in one partition. Pairs can share proteins, families and ancestry; no independence, causal effect, separate-run stability or evolutionary-event interpretation.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
