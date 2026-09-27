#!/usr/bin/env python3
"""Prepare fixed-tree sequence covariates for provisionally referenced duplicates."""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def derive(row):
    """Recover tip edges from shared reference paths; preserve numerical limits."""
    da, db, pair, shared = (float(row[k]) for k in (
        'distance_a_to_reference', 'distance_b_to_reference',
        'duplicate_pair_sequence_distance', 'reference_distance_from_duplicate_node'))
    if not all(math.isfinite(x) and x >= 0 for x in (da, db, pair, shared)):
        raise ValueError('Invalid source sequence distance')
    # Exported paths are sums of floating-point branch lengths. Cancellation
    # cannot support interpreting a tiny residual as an evolutionary difference.
    tolerance = 1e-12 * max(1., da, db, pair, shared)
    a, b = da - shared, db - shared
    if min(a, b) < -tolerance or abs(a + b - pair) > 4 * tolerance:
        raise ValueError('Incompatible shared-reference and duplicate paths')
    clamped = int(a < 0) + int(b < 0)
    a, b = max(0., a), max(0., b)
    signed = a - b
    if pair <= 4 * tolerance:
        state, normalized, direction = 'unresolved_near_zero_pair_distance', '', 'unresolved'
    else:
        normalized = signed / pair
        if abs(normalized) > 1 + 4 * tolerance / pair:
            raise ValueError('Impossible normalized contrast')
        state = 'resolved_pair_distance'
        direction = 'a_longer' if signed > 2 * tolerance else 'b_longer' if signed < -2 * tolerance else 'unresolved'
    return dict(sequence_tip_a=a, sequence_tip_b=b, sequence_signed_difference=signed,
                sequence_normalized_contrast=normalized, sequence_direction=direction,
                sequence_covariate_status=state, distance_roundoff_tolerance=tolerance,
                negative_roundoff_clamps=clamped)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, required=True)
    parser.add_argument('--readback', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    producer_path = args.inventory / 'receipt.json'
    audit_path = args.readback / 'receipt.json'
    producer = json.loads(producer_path.read_text())
    audit = json.loads(audit_path.read_text())
    if producer['status'] != 'complete_duplication_sister_reference_inventory' or audit['status'] != 'passed_full_duplication_sister_reference_readback' or audit['producer_receipt_sha256'] != sha(producer_path):
        raise ValueError('Missing or mismatched complete source readback')
    pins = {str(producer_path): sha(producer_path), str(audit_path): sha(audit_path)}
    for root, receipt in [(args.inventory, producer), (args.readback, audit)]:
        for name, digest in receipt['artifacts'].items():
            path = root / name
            if sha(path) != digest:
                raise ValueError(f'Changed source: {path}')
            pins[str(path)] = digest
    args.output.mkdir(parents=True, exist_ok=False)
    totals = {}
    rows_by_guide = {}
    for guide in ['profile', 'mafft']:
        counts, records = Counter(), {}
        with (args.inventory / f'{guide}_sister_references.tsv').open() as handle, (args.output / f'{guide}_sequence_covariates.tsv').open('w') as out:
            reader = csv.DictReader(handle, delimiter='\t')
            writer = None
            for row in reader:
                if row['status'] != 'provisional_reference_available':
                    continue
                derived = derive(row)
                key = tuple(sorted([row['gene_a'], row['gene_b']]))
                if key in records:
                    raise ValueError('Repeated gene pair')
                if row['gene_a'] != key[0]:
                    raise ValueError('Noncanonical source gene order')
                record = dict(row, **derived)
                records[key] = record
                counts[derived['sequence_covariate_status']] += 1
                counts[derived['sequence_direction']] += 1
                if writer is None:
                    writer = csv.DictWriter(out, fieldnames=list(record), delimiter='\t')
                    writer.writeheader()
                writer.writerow(record)
        expected = next(x for x in producer['guides'] if x['guide'] == guide)['counts']['provisional_reference_available']
        if len(records) != expected:
            raise ValueError('Candidate scope changed')
        totals[guide] = dict(rows=len(records), counts=dict(counts))
        rows_by_guide[guide] = records
    comparisons = Counter()
    path = args.output / 'guide_sequence_covariate_comparison.tsv'
    with path.open('w') as handle:
        fields = ['gene_a', 'gene_b', 'profile_family', 'mafft_family', 'reference_sets_equal',
                  'profile_status', 'mafft_status', 'profile_direction', 'mafft_direction',
                  'profile_normalized_contrast', 'mafft_normalized_contrast', 'comparison']
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
        for key in sorted(rows_by_guide['profile'].keys() & rows_by_guide['mafft'].keys()):
            a, b = (rows_by_guide[g][key] for g in ['profile', 'mafft'])
            same = set(json.loads(a['nearest_reference_genes'])) == set(json.loads(b['nearest_reference_genes']))
            directions = [x['sequence_direction'] for x in (a, b)]
            outcome = 'unresolved_in_one_or_both_guides' if 'unresolved' in directions else 'same_direction' if directions[0] == directions[1] else 'opposite_direction'
            comparisons[outcome] += 1
            writer.writerow(dict(gene_a=key[0], gene_b=key[1], profile_family=a['family'], mafft_family=b['family'],
                reference_sets_equal=int(same), profile_status=a['sequence_covariate_status'], mafft_status=b['sequence_covariate_status'],
                profile_direction=directions[0], mafft_direction=directions[1], profile_normalized_contrast=a['sequence_normalized_contrast'],
                mafft_normalized_contrast=b['sequence_normalized_contrast'], comparison=outcome))
    for path, digest in pins.items():
        if sha(path) != digest:
            raise ValueError('Source changed during execution')
    receipt = dict(status='complete_duplication_sequence_covariates', source_sha256=pins,
        script_sha256=sha(__file__), guides=totals, shared_pair_comparison=dict(comparisons),
        artifacts={p.name:sha(p) for p in sorted(args.output.glob('*.tsv'))},
        scope='Fixed sequence-tree relative divergence covariates for all provisionally referenced pairs. Lexical gene orientation, no calendar-time rate, ancestral state, structural asymmetry, significance test or biological orthology validation. Near-zero values remain unresolved; guide outputs are sensitivity analyses, not independent observations.')
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k:receipt[k] for k in ['status','guides','shared_pair_comparison']}, indent=2))


if __name__ == '__main__':
    main()
