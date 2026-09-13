#!/usr/bin/env python3
"""Concatenate audited markers with explicit missing-data and coordinate provenance."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
AA = set('ACDEFGHIKLMNPQRSTVWY')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def concatenate(taxa, markers):
    chunks = {name: [] for name in taxa}
    partitions, sites, coverage = [], [], {name: 0 for name in taxa}
    offset = 0
    for marker, sequences, columns in markers:
        if set(sequences) - set(taxa):
            raise ValueError('Unexpected taxon in marker')
        if len({len(seq) for seq in sequences.values()}) != 1:
            raise ValueError('Nonrectangular marker alignment')
        length = len(next(iter(sequences.values())))
        if columns != sorted(set(columns)) or any(p < 1 or p > length for p in columns):
            raise ValueError('Invalid 1-based column mask')
        if not columns:
            continue
        for name in taxa:
            seq = sequences.get(name)
            if seq is None:
                value = '-' * len(columns)
            else:
                value = ''.join(seq[p - 1] if seq[p - 1] in AA or seq[p - 1] == '-' else 'X' for p in columns)
            chunks[name].append(value)
            coverage[name] += sum(aa in AA for aa in value)
        partitions.append((marker, offset + 1, offset + len(columns)))
        sites.extend((offset + i, marker, p) for i, p in enumerate(columns, 1))
        offset += len(columns)
    if not offset or any(n == 0 for n in coverage.values()):
        raise ValueError('Empty matrix or taxon with no unambiguous residues')
    matrix = {name: ''.join(parts) for name, parts in chunks.items()}
    if {len(s) for s in matrix.values()} != {offset}:
        raise ValueError('Concatenation length mismatch')
    return matrix, partitions, sites, coverage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alignments', required=True, type=Path)
    parser.add_argument('--audit', required=True, type=Path)
    parser.add_argument('--occupancy', choices=['0.25', '0.5', '0.75'], default='0.5')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    audit = json.loads((args.audit / 'receipt.json').read_text())
    if audit['status'] != 'complete_audit':
        raise ValueError('Full alignment audit required')
    with (ROOT / 'metadata/analysis_manifest.tsv').open() as handle:
        taxa = {row['taxon_id']: row for row in csv.DictReader(handle, delimiter='\t')}
    masks = json.loads((args.audit / 'column_masks.json').read_text())
    expected = {row['marker']: row['alignment_sha256'] for row in audit['alignment_receipts']}
    if set(masks) != set(expected):
        raise ValueError('Mask/receipt marker mismatch')
    markers = []
    for marker in sorted(expected):
        path = args.alignments / f'{marker}.faa'
        if sha(path) != expected[marker]:
            raise ValueError('Changed audited alignment')
        with path.open() as handle:
            records = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        sequences = {name: str(record.seq).upper() for name, record in records.items()}
        selected = masks[marker][args.occupancy]
        # Recompute the declared rule before accepting a serialized mask.
        recomputed = [i for i, column in enumerate(zip(*sequences.values()), 1)
                      if sum(aa in AA for aa in column) / len(sequences) >= float(args.occupancy)]
        if selected != recomputed:
            raise ValueError('Mask differs from its declared occupancy rule')
        markers.append((marker, sequences, selected))
    matrix, partitions, sites, coverage = concatenate(sorted(taxa), markers)
    if args.output.exists():
        raise FileExistsError('Use an immutable new matrix directory')
    args.output.mkdir(parents=True)
    (args.output / 'matrix.faa').write_text(''.join(f'>{name}\n{seq}\n' for name, seq in matrix.items()))
    (args.output / 'partitions.nex').write_text('#nexus\nbegin sets;\n' + ''.join(
        f'  charset marker_{marker} = {start}-{end};\n' for marker, start, end in partitions) + 'end;\n')
    with (args.output / 'site_mapping.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['matrix_column_1based', 'marker', 'alignment_column_1based'])
        writer.writerows(sites)
    with (args.output / 'taxon_coverage.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['taxon_id', 'species_name', 'study_role', 'unambiguous_residues', 'unambiguous_fraction'])
        for name in sorted(taxa):
            writer.writerow([name, taxa[name]['species_name'], taxa[name]['study_role'], coverage[name], coverage[name] / len(matrix[name])])
    result = {'taxa': len(taxa), 'markers': len(partitions), 'columns': len(sites),
              'occupancy_threshold': args.occupancy, 'ambiguous_residues': 'normalized to X',
              'manifest_sha256': sha(ROOT / 'metadata/analysis_manifest.tsv'),
              'audit_receipt_sha256': sha(args.audit / 'receipt.json'),
              'mask_sha256': sha(args.audit / 'column_masks.json'),
              'artifacts': {p.name: sha(p) for p in sorted(args.output.iterdir())}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
