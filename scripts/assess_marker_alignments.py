#!/usr/bin/env python3
"""Audit verified alignments and record transparent column-occupancy masks."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')


def assess(sequences, thresholds=(.25, .5, .75)):
    if not sequences or len({len(s) for s in sequences}) != 1:
        raise ValueError('Expected nonempty rectangular alignment')
    length = len(sequences[0])
    if length == 0:
        raise ValueError('Empty alignment')
    masks = {str(t): [] for t in thresholds}
    informative, constant, all_unknown = 0, 0, 0
    for position, column in enumerate(zip(*sequences), 1):
        observed = Counter(aa for aa in column if aa in AMINO_ACIDS)
        occupied = sum(observed.values()) / len(sequences)
        for threshold in thresholds:
            if occupied >= threshold:
                masks[str(threshold)].append(position)
        if not observed:
            all_unknown += 1
        elif len(observed) == 1:
            constant += 1
        if sum(count >= 2 for count in observed.values()) >= 2:
            informative += 1
    return {'taxa': len(sequences), 'columns': length, 'constant_columns': constant,
            'parsimony_informative_columns': informative, 'all_gap_or_ambiguous_columns': all_unknown,
            'unambiguous_residue_fraction': sum(aa in AMINO_ACIDS for seq in sequences for aa in seq) / (len(sequences) * length)}, masks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alignments', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--allow-incomplete', action='store_true')
    args = parser.parse_args()
    complete_path = args.alignments / 'receipt.json'
    if not complete_path.exists() and not args.allow_incomplete:
        raise ValueError('Alignment batch is incomplete; use staging mode only for inspection')
    if args.output.exists():
        raise FileExistsError('Use a new audit snapshot output directory')
    with (ROOT / 'metadata/analysis_manifest.tsv').open() as handle:
        taxa = {r['taxon_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    rows, coverage, masks, receipts = [], [], {}, []
    for path in sorted(args.alignments.glob('*.receipt.json')):
        receipt = json.loads(path.read_text())
        if receipt['status'] != 'aligned':
            continue
        marker = receipt['marker']
        alignment = args.alignments / f'{marker}.faa'
        if hashlib.sha256(alignment.read_bytes()).hexdigest() != receipt['output_sha256']:
            raise ValueError(f'Changed alignment: {marker}')
        with alignment.open() as handle:
            records = list(SeqIO.parse(handle, 'fasta'))
        if len({r.id for r in records}) != len(records) or any(r.id not in taxa for r in records):
            raise ValueError('Duplicate or unexpected taxon')
        sequences = [str(r.seq).upper() for r in records]
        result, retained = assess(sequences)
        result['marker'] = marker
        for threshold, positions in retained.items():
            result['columns_occupancy_' + threshold] = len(positions)
        rows.append(result)
        masks[marker] = retained
        for record, seq in zip(records, sequences):
            coverage.append({'marker': marker, 'taxon_id': record.id,
                             'unambiguous_residues': sum(aa in AMINO_ACIDS for aa in seq),
                             'residues_after_50pct_mask': sum(seq[p - 1] in AMINO_ACIDS for p in retained['0.5'])})
        receipts.append({'marker': marker, 'alignment_sha256': receipt['output_sha256']})
    if not rows:
        raise ValueError('No verified alignments available')
    if complete_path.exists():
        expected = {r['marker'] for r in json.loads(complete_path.read_text())['alignments'] if r['status'] == 'aligned'}
        if expected != set(masks):
            raise ValueError('Alignment batch receipt and audit inputs differ')
    args.output.mkdir(parents=True)
    for name, records in [('marker_statistics.tsv', rows), ('taxon_marker_coverage.tsv', coverage)]:
        with (args.output / name).open('w') as handle:
            writer = csv.DictWriter(handle, list(records[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(records)
    (args.output / 'column_masks.json').write_text(json.dumps(masks) + '\n')
    result = {'status': 'complete_audit' if complete_path.exists() else 'incomplete_staging',
              'markers': len(rows), 'planned_taxa': len(taxa), 'thresholds': [.25, .5, .75],
              'alignment_receipts': receipts,
              'mask_coordinate_system': '1-based positions in untrimmed alignment',
              'interpretation': 'Occupancy counts unambiguous amino acids among taxa present in each marker. This is not alignment-confidence estimation or a final filtering decision.'}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Audited', len(rows), 'alignments;', result['status'])


if __name__ == '__main__':
    main()
