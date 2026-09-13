#!/usr/bin/env python3
"""Build paired AA/3Di marker alignments from audited coordinate-derived states."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from compare_marker_structures import ROOT, AA, sha
from assess_pae_sensitivity import checked_receipt


def write_table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def paired_row(sequence, columns, positions, encoding):
    """Keep physical correspondence and reject invalid or uncertain feature contexts."""
    amino, structural, reasons = [], [], Counter()
    for column in columns:
        aa = sequence[column - 1]
        reason = ''
        if aa not in AA:
            reason = 'noncanonical_or_missing_sequence'
        elif column not in positions:
            reason = 'no_structural_mapping'
        else:
            position, confidence = positions[column]
            i = position - 1
            if encoding is None or not 0 <= i < len(encoding['sequence']):
                raise ValueError('Missing encoding or out-of-range mapped residue')
            if encoding['sequence'][i] != aa or not np.isclose(encoding['ca_plddt'][i], confidence):
                raise ValueError('Sequence or confidence correspondence differs')
            if not encoding['valid'][i]:
                reason = 'invalid_native_feature'
            elif not np.isfinite(encoding['feature_min_plddt'][i]) or encoding['feature_min_plddt'][i] < 70:
                reason = 'low_feature_plddt'
            elif not np.isfinite(encoding['feature_max_pae'][i]) or encoding['feature_max_pae'][i] > 10:
                reason = 'high_feature_pae'
            elif encoding['states'][i] not in AA:
                raise ValueError('Unknown structural state')
        if reason:
            amino.append('?'); structural.append('?'); reasons[reason] += 1
        else:
            amino.append(aa); structural.append(encoding['states'][i]); reasons['observed'] += 1
    return ''.join(amino), ''.join(structural), reasons


def site_counts(sequences):
    variable = informative = 0
    for column in zip(*sequences):
        counts = Counter(c for c in column if c != '?')
        variable += len(counts) > 1
        informative += sum(v >= 2 for v in counts.values()) >= 2
    return variable, informative


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['encodings', 'snapshot', 'matrix', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    receipts = {name: checked_receipt(getattr(args, name)) for name in ['encodings', 'snapshot', 'matrix']}
    if receipts['encodings']['mapping_receipt_sha256'] != sha(args.snapshot / 'receipt.json'):
        raise ValueError('Encoding and mapping snapshots differ')
    if receipts['snapshot']['matrix_receipt_sha256'] != sha(args.matrix / 'receipt.json'):
        raise ValueError('Mapping and alignment snapshots differ')
    if args.output.exists():
        raise FileExistsError('Use a new immutable output directory')
    matrix_records = list(SeqIO.parse(args.matrix / 'matrix.faa', 'fasta'))
    matrix = {r.id: str(r.seq) for r in matrix_records}
    if len(matrix) != len(matrix_records) or len({len(s) for s in matrix.values()}) != 1:
        raise ValueError('Duplicate taxa or nonrectangular matrix')
    markers = defaultdict(list)
    source_sites = {}
    for row in csv.DictReader((args.matrix / 'site_mapping.tsv').open(), delimiter='\t'):
        column = int(row['matrix_column_1based'])
        if column in source_sites:
            raise ValueError('Duplicate matrix column')
        markers[row['marker']].append(column)
        source_sites[column] = row
    if sorted(source_sites) != list(range(1, len(next(iter(matrix.values()))) + 1)):
        raise ValueError('Incomplete site map')
    encoded = {}
    for row in csv.DictReader((args.encodings / 'model_summary.tsv').open(), delimiter='\t'):
        path = ROOT / row['encoding_path']
        if sha(path) != row['encoding_sha256'] or row['model_name'] in encoded:
            raise ValueError('Changed or duplicate encoding')
        with np.load(path, allow_pickle=False) as source:
            data = {k: source[k].copy() for k in source.files}
        data['sequence'], data['states'] = str(data['sequence']), str(data['states'])
        if hashlib.sha256(data['sequence'].encode()).hexdigest() != row['sequence_sha256']:
            raise ValueError('Encoding sequence checksum differs')
        if len(data['states']) != len(data['sequence']):
            raise ValueError('State length mismatch')
        encoded[row['model_name']] = data
    links = {}
    for row in csv.DictReader((args.snapshot / 'marker_structure_links.tsv').open(), delimiter='\t'):
        key = row['marker'], row['taxon_id']
        name = Path(row['model_path']).stem
        if key in links or key[0] not in markers or key[1] not in matrix:
            raise ValueError('Duplicate or unknown marker/taxon link')
        if sha(ROOT / row['model_path']) != row['model_sha256']:
            raise ValueError('Changed model coordinates')
        if hashlib.sha256(encoded[name]['sequence'].encode()).hexdigest() != row['sequence_sha256']:
            raise ValueError('Linked model sequence differs')
        links[key] = name
    mapped = defaultdict(dict)
    with gzip.open(args.snapshot / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['marker'], row['taxon_id']
            column = int(row['matrix_column_1based'])
            if key not in links or column not in markers[key[0]] or column in mapped[key]:
                raise ValueError('Invalid or duplicate residue map')
            mapped[key][column] = int(row['protein_residue_1based']), float(row['ca_plddt'])
    args.output.mkdir(parents=True)
    coverage, summary = [], []
    all_reasons = ['observed', 'noncanonical_or_missing_sequence', 'no_structural_mapping',
                   'invalid_native_feature', 'low_feature_plddt', 'high_feature_pae']
    for marker, columns in markers.items():
        minimum = max(50, math.ceil(.3 * len(columns)))
        eligible = {}
        for taxon, sequence in matrix.items():
            key = marker, taxon
            name = links.get(key)
            a, s, reasons = paired_row(sequence, columns, mapped.get(key, {}), encoded.get(name))
            if sum(reasons.values()) != len(columns):
                raise ValueError('Site accounting differs')
            include = reasons['observed'] >= minimum
            coverage.append({'marker': marker, 'taxon_id': taxon, 'model_name': name or '',
                'original_marker_columns': len(columns), 'required_observed_columns': minimum,
                **{reason: reasons[reason] for reason in all_reasons}, 'taxon_eligible': include})
            if include:
                eligible[taxon] = (a, s)
        ready = len(eligible) >= 4
        kept = [i for i in range(len(columns)) if any(pair[0][i] != '?' for pair in eligible.values())]
        stats = {}
        if eligible:
            for channel, index in [('aa', 0), ('3di', 1)]:
                sequences = [''.join(pair[index][i] for i in kept) for pair in eligible.values()]
                variable, informative = site_counts(sequences)
                stats[channel + '_variable_columns'] = variable
                stats[channel + '_parsimony_informative_columns'] = informative
        else:
            stats = {channel + suffix: 0 for channel in ['aa', '3di']
                     for suffix in ['_variable_columns', '_parsimony_informative_columns']}
        if ready:
            folder = args.output / marker
            folder.mkdir()
            for channel, index in [('aa', 0), ('3di', 1)]:
                with (folder / (channel + '.faa')).open('w') as handle:
                    for taxon, pair in eligible.items():
                        handle.write('>' + taxon + '\n' + ''.join(pair[index][i] for i in kept) + '\n')
            write_table(folder / 'columns.tsv', [dict(paired_column_1based=j + 1, **source_sites[columns[i]])
                                               for j, i in enumerate(kept)])
            # Independent FASTA readback checks the emitted files, not only in-memory strings.
            a = {r.id: str(r.seq) for r in SeqIO.parse(folder / 'aa.faa', 'fasta')}
            s = {r.id: str(r.seq) for r in SeqIO.parse(folder / '3di.faa', 'fasta')}
            if set(a) != set(s) or set(a) != set(eligible):
                raise ValueError('Paired file taxon mismatch')
            for taxon in a:
                if len(a[taxon]) != len(kept) or len(s[taxon]) != len(kept):
                    raise ValueError('Paired file length mismatch')
                if [c == '?' for c in a[taxon]] != [c == '?' for c in s[taxon]]:
                    raise ValueError('Paired file observation masks differ')
        summary.append({'marker': marker, 'total_taxa_audited': len(matrix),
            'original_marker_columns': len(columns), 'eligible_taxa': len(eligible),
            'retained_columns': len(kept), 'status': 'ready_for_inference' if ready else 'insufficient_structural_coverage', **stats})
    write_table(args.output / 'taxon_coverage.tsv', coverage)
    write_table(args.output / 'marker_summary.tsv', summary)
    receipt = {'status': 'complete_paired_phylogenetic_input_preparation',
        'taxa_audited': len(matrix), 'markers_audited': len(markers), 'taxon_marker_rows': len(coverage),
        'ready_markers': sum(r['status'] == 'ready_for_inference' for r in summary),
        'source_receipts': {name: {'path': str(getattr(args, name)), 'sha256': sha(getattr(args, name) / 'receipt.json')} for name in receipts},
        'script_sha256': sha(Path(__file__)),
        'mask': 'Canonical AA with exact structure correspondence; native valid state, all six feature residues pLDDT>=70, maximum directional feature-context PAE<=10. The AA and 3Di masks are identical.',
        'eligibility': 'At least max(50,ceil(0.3*original marker columns)) observations per taxon; at least four taxa per marker. Keep invariant sites; remove only all-missing columns among eligible taxa.',
        'interpretation': 'Coverage-limited acquisition checkpoint within full sampling; no branch inference or evolutionary claims yet. ? means unobserved, never a structural state. Alphabet models and paired branch uncertainty remain required.',
        'artifacts': {str(p.relative_to(args.output)): sha(p) for p in sorted(args.output.rglob('*')) if p.is_file()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'artifacts'}, indent=2))


if __name__ == '__main__':
    main()
