#!/usr/bin/env python3
"""Independently read candidate Stockholm blocks and verify frozen-site projections."""
import argparse
import csv
import hashlib
import json
import statistics
from collections import defaultdict
from pathlib import Path

AA = set('ACDEFGHIKLMNPQRSTVWY')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def fasta(path):
    result = {}
    ident = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            ident = line[1:].split()[0]
            if ident in result:
                raise ValueError('Duplicate FASTA identifier')
            result[ident] = ''
        elif line.strip():
            if ident is None:
                raise ValueError('Sequence before identifier')
            result[ident] += line.strip()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--candidates', type=Path, required=True)
    parser.add_argument('--matrix', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    receipt = json.loads((args.source / 'receipt.json').read_text())
    for name, digest in receipt['artifacts'].items():
        assert sha(args.source / name) == digest
    assert sha(args.candidates / 'receipt.json') == receipt['candidate_source_receipt_sha256']
    assert sha(args.matrix / 'receipt.json') == receipt['baseline_matrix_receipt_sha256']
    originals = fasta(args.candidates / 'candidate_proteins.faa')
    baseline = fasta(args.matrix / 'matrix.faa')
    sites = defaultdict(list)
    for row in rows(args.matrix / 'site_mapping.tsv'):
        sites[row['marker']].append(row)
    mapping = defaultdict(list)
    for row in rows(args.source / 'candidate_residue_mapping.tsv'):
        mapping[row['review_sequence_id']].append(row)
    summary_rows = rows(args.source / 'candidate_alignment_review.tsv')
    summaries = {r['review_sequence_id']: r for r in summary_rows}
    assert len(summaries) == len(summary_rows) == len(originals)
    checked = observed = baseline_characters = 0
    dispositions = []
    seen = set()
    for proof in receipt['marker_proofs']:
        marker = proof['marker']
        folder = args.source / marker
        for name, digest in proof['artifacts'].items():
            assert sha(folder / name) == digest
        # Parse interleaved Stockholm blocks directly, without Bio.AlignIO.
        sequences = defaultdict(str)
        rf = ''
        lines = (folder / 'candidates.sto').read_text().splitlines()
        assert lines[0] == '# STOCKHOLM 1.0' and lines[-1] == '//'
        for line in lines[1:-1]:
            if line.startswith('#=GC RF '):
                rf += line.split()[2]
            elif line and not line.startswith('#'):
                ident, fragment = line.split()
                sequences[ident] += fragment
        assert all(len(s) == len(rf) for s in sequences.values())
        projected = fasta(folder / 'candidate_projected.faa')
        augmented = fasta(folder / 'baseline_plus_candidates.faa')
        assert set(augmented) == set(baseline) | set(sequences)
        assert set(projected) == set(sequences)
        for taxon, sequence in baseline.items():
            expected = ''.join(sequence[int(r['matrix_column_1based']) - 1] for r in sites[marker])
            assert augmented[taxon] == expected
            baseline_characters += len(expected)
            dispositions.append({'marker': marker, 'sequence_id': taxon, 'kind': 'baseline',
                                 'observed_residues': sum(c in AA for c in expected),
                                 'disposition': 'retain_for_review' if any(c in AA for c in expected) else 'omit_all_missing_before_tree_inference'})
        for ident, aligned in sequences.items():
            assert ident not in seen
            seen.add(ident)
            assert aligned.upper().replace('.', '').replace('-', '') == originals[ident]
            assert len(mapping[ident]) == len(sites[marker])
            states = {}
            residue = state = 0
            for character, reference in zip(aligned.upper(), rf):
                present = character not in '.-'
                residue += present
                if reference not in '.- ':
                    state += 1
                    states[state] = (character if present else '-', residue if present else '')
            expected = ''
            for index, (site, row) in enumerate(zip(sites[marker], mapping[ident]), 1):
                aa, position = states[int(site['alignment_column_1based'])]
                analysis = aa if aa in AA else '?'
                assert row == {'marker': marker, 'review_sequence_id': ident,
                               'marker_matrix_column_1based': str(index),
                               'matrix_column_1based': site['matrix_column_1based'],
                               'profile_match_state_1based': site['alignment_column_1based'],
                               'protein_residue_1based': str(position), 'amino_acid': aa,
                               'analysis_state': analysis}
                if position:
                    assert originals[ident][position - 1] == aa
                expected += analysis
                checked += 1
                observed += analysis != '?'
            assert expected == projected[ident] == augmented[ident]
            count = sum(c in AA for c in expected)
            summary = summaries[ident]
            assert int(summary['observed_canonical_residues']) == count
            assert int(summary['matrix_marker_columns']) == len(expected)
            assert abs(float(summary['observed_fraction']) - count / len(expected)) < 1e-12
            dispositions.append({'marker': marker, 'sequence_id': ident, 'kind': 'alternative_gene_copy',
                                 'observed_residues': count,
                                 'disposition': 'retain_for_review' if count else 'omit_all_missing_before_tree_inference'})
    assert seen == set(originals) == set(mapping) == set(summaries)
    assert checked == receipt['projected_candidate_site_rows']
    assert observed == receipt['observed_canonical_candidate_residues']
    fractions = [float(r['observed_fraction']) for r in summary_rows]
    args.output.mkdir(parents=True)
    path = args.output / 'sequence_disposition.tsv'
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(dispositions[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(dispositions)
    result = {'status': 'complete_independent_stockholm_and_projection_readback',
              'source_receipt_sha256': sha(args.source / 'receipt.json'),
              'script_sha256': sha(Path(__file__)), 'candidates': len(seen),
              'markers': len(receipt['marker_proofs']), 'candidate_site_rows_checked': checked,
              'canonical_candidate_residues': observed, 'baseline_characters_checked': baseline_characters,
              'candidate_observed_fraction_min_median_max': [min(fractions), statistics.median(fractions), max(fractions)],
              'all_missing_baseline_rows': sum(r['kind'] == 'baseline' and r['observed_residues'] == 0 for r in dispositions),
              'all_missing_candidates': sum(r['kind'] == 'alternative_gene_copy' and r['observed_residues'] == 0 for r in dispositions),
              'artifacts': {path.name: sha(path)},
              'interpretation': 'All sequences and site mappings checked. Disposition identifies zero-information rows for later tree preparation; retained candidates remain unaccepted gene-copy hypotheses.'}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
