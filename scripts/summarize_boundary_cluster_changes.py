#!/usr/bin/env python3
"""Describe cluster changes by boundary extension, without treating hits as independent."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path

from screen_duplication_alignment_reuse import sha


def extension_bin(added):
    if added == 0:
        return '0'
    for upper, label in [(4, '1-4'), (9, '5-9'), (19, '10-19'), (49, '20-49')]:
        if added <= upper:
            return label
    return '50+'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--lengths', type=Path, required=True)
    parser.add_argument('--clusters', type=Path, required=True)
    parser.add_argument('--cluster-check', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    length_receipt = args.lengths / 'receipt.json'
    length_check = args.lengths / 'readback.json'
    cluster_receipt = args.clusters / 'receipt.json'
    lr, lc, cr, cc = [json.loads(p.read_text()) for p in
                      [length_receipt, length_check, cluster_receipt, args.cluster_check]]
    assert lr['status'] == 'complete_full_candidate_domain_boundary_sensitivity'
    assert lc['status'] == 'passed_all_boundary_pair_rows_and_threshold_counts'
    assert lc['producer_receipt_sha256'] == sha(length_receipt)
    assert cr['status'] == 'complete_all_domain_boundary_cluster_dispositions'
    assert cc['status'] == 'passed_full_domain_boundary_cluster_output_readback'
    assert cc['source_hashes'][str(cluster_receipt)] == sha(cluster_receipt)
    pins = {str(p): sha(p) for p in [length_receipt, length_check, cluster_receipt, args.cluster_check]}
    pins[str(Path(__file__))] = sha(__file__)
    for root, receipt in [(args.lengths, lr), (args.clusters, cr)]:
        pins.update({str(root / name): h for name, h in receipt['artifacts'].items()})
    def verify():
        for p, h in pins.items():
            assert sha(p) == h, p
    verify()
    assert extension_bin(0) == '0' and extension_bin(4) == '1-4'
    assert extension_bin(5) == '5-9' and extension_bin(49) == '20-49' and extension_bin(50) == '50+'
    lengths = {}
    with gzip.open(args.lengths / 'boundary_sensitivity.tsv.gz', 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            key = r['model_key'], r['hit_id']
            assert key not in lengths
            a, e = int(r['alignment_residues']), int(r['envelope_residues'])
            added = int(r['added_residues'])
            assert added == e-a >= 0 and a > 0
            lengths[key] = (added, a)
    seen, counts, bins, models = set(), Counter(), defaultdict(Counter), defaultdict(set)
    with gzip.open(args.clusters / 'boundary_cluster_dispositions.tsv.gz', 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            key = r['model_key'], r['hit_id']
            assert key not in seen and key in lengths
            seen.add(key)
            added, alignment_length = lengths[key]
            status = r['disposition']
            assert status in ['identical_interval', 'distinct_intervals_same_cluster', 'distinct_intervals_different_clusters']
            assert (added == 0) == (status == 'identical_interval')
            counts[status] += 1
            label = extension_bin(added)
            bins[label]['hits'] += 1
            bins[label]['changed_cluster'] += status == 'distinct_intervals_different_clusters'
            models[label].add(key[0])
    assert seen == set(lengths)
    assert len(seen) == lr['candidate_model_hit_pairs'] == cr['candidate_pairs'] == cc['candidate_pairs']
    assert dict(counts) == cr['dispositions'] == cc['dispositions']
    rows = [dict(added_residue_bin=label, **bins[label], models=len(models[label]),
                 changed_fraction=bins[label]['changed_cluster']/bins[label]['hits'])
            for label in ['0', '1-4', '5-9', '10-19', '20-49', '50+'] if bins[label]['hits']]
    assert sum(r['hits'] for r in rows) == len(seen)
    assert sum(r['changed_cluster'] for r in rows) == counts['distinct_intervals_different_clusters']
    verify()
    args.output.mkdir(parents=True, exist_ok=False)
    table = args.output / 'cluster_change_by_extension.tsv'
    with table.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    result = dict(status='complete_descriptive_boundary_extension_cluster_change_join',
                  candidate_hits=len(seen), dispositions=dict(counts), source_hashes=pins,
                  artifacts={table.name: sha(table)}, summaries=rows,
                  scope='Every checked model/hit matched once between audited boundary lengths and cluster assignments. Descriptive reporting bins, not exclusion thresholds. Models can occur in multiple bins; hits and model links are not independent replicates. Within-partition assignment sensitivity, not direct structural change, causal boundary effect, orthology or evolutionary inference.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'source_hashes'}))


if __name__ == '__main__':
    main()
