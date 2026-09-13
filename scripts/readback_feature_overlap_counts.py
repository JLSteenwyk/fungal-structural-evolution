#!/usr/bin/env python3
"""Verify full feature-overlap alignment grids and totals against paired FASTAs."""
import argparse
import json
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['overlap', 'inputs', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    r = checked_receipt(a.overlap); checked_receipt(a.inputs)
    if r['source_receipts']['inputs']['sha256'] != sha(a.inputs / 'receipt.json'):
        raise ValueError('Input receipt mismatch')
    expected = {}
    for path in a.inputs.glob('*/3di.faa'):
        for record in SeqIO.parse(path, 'fasta'):
            key = path.parent.name, record.id
            if key in expected:
                raise ValueError('Duplicate FASTA identity')
            expected[key] = len(record.seq), sum(c != '?' for c in record.seq)
    rows = read_table(a.overlap / 'feature_overlap.tsv')
    seen = set()
    for row in rows:
        key = row['marker'], row['taxon_id']
        if key in seen or key not in expected:
            raise ValueError('Unexpected alignment identity')
        seen.add(key); length, observed = expected[key]
        if (length, observed) != (int(row['alignment_columns']), int(row['observed_features'])):
            raise ValueError('FASTA dimensions or state counts differ')
        total = int(row['shared_coordinate_feature_pairs'])
        b10 = int(row['overlap_pairs_beyond_circular_block_10']); b30 = int(row['overlap_pairs_beyond_circular_block_30'])
        if not 0 <= b30 <= b10 <= total or not 1 <= int(row['largest_overlap_component']) <= observed or not 1 <= int(row['overlap_components']) <= observed:
            raise ValueError('Invalid graph count bounds')
        for block, count in [(10, b10), (30, b30)]:
            reported = row[f'fraction_overlap_pairs_beyond_circular_block_{block}']
            if total and abs(float(reported) - count / total) > 1e-12:
                raise ValueError('Fraction mismatch')
            if not total and reported != '':
                raise ValueError('Nonempty fraction without denominator')
    if seen != set(expected) or len(rows) != r['taxon_marker_alignments'] or len({k[0] for k in seen}) != r['markers']:
        raise ValueError('Incomplete alignment universe')
    for field in ['observed_features', 'shared_coordinate_feature_pairs']:
        if sum(int(row[field]) for row in rows) != r[field]:
            raise ValueError('Receipt total mismatch')
    for block in [10, 30]:
        if sum(int(row[f'overlap_pairs_beyond_circular_block_{block}']) for row in rows) != r['overlap_pairs_beyond_circular_blocks'][str(block)]:
            raise ValueError('Block count total mismatch')
    a.output.mkdir(parents=True)
    result = {'status': 'passed_complete_feature_overlap_counts_readback', 'source_receipt_sha256': sha(a.overlap / 'receipt.json'),
              'script_sha256': sha(Path(__file__)), 'alignments': len(rows), 'observed_features': r['observed_features'],
              'block_fractions': {b: n / r['shared_coordinate_feature_pairs'] for b, n in r['overlap_pairs_beyond_circular_blocks'].items()},
              'interpretation': 'Every identity, FASTA dimension/state count, summary total, fraction and graph count bound checked. Not an independent full graph reconstruction or covariance estimate.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
