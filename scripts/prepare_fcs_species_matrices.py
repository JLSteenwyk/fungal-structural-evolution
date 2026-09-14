#!/usr/bin/env python3
"""Mask FCS EXCLUDE/FIX/TRIM marker observations in unchanged full species matrices."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table

AA = set('ACDEFGHIKLMNPQRSTVWY')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['matrix', 'mapping', 'audit', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    baseline = checked_receipt(a.matrix)
    checked_receipt(a.mapping)
    audit = json.loads((a.audit / 'receipt.json').read_text())
    if audit['status'] != 'passed_full_region_cds_intersection_and_protein_marker_join_audit' or audit['mapping_receipt_sha256'] != sha(a.mapping / 'receipt.json'):
        raise ValueError('Full FCS mapping audit required')
    columns = defaultdict(list)
    for row in read_table(a.matrix / 'site_mapping.tsv'):
        columns[row['marker']].append(int(row['matrix_column_1based']) - 1)
    original_records = list(SeqIO.parse(a.matrix / 'matrix.faa', 'fasta'))
    original = {r.id: str(r.seq) for r in original_records}
    if len(original) != len(original_records) or len(original) != baseline['taxa']:
        raise ValueError('Taxon grid differs')
    if any(len(s) != baseline['columns'] for s in original.values()):
        raise ValueError('Matrix is not rectangular')
    transformed = {t: list(s) for t, s in original.items()}
    flagged = []
    mask_positions = defaultdict(set)
    observed_removed = 0
    for row in read_table(a.mapping / 'marker_overlap_review.tsv'):
        actions = set(row['fcs_actions'].split(';'))
        if not actions & {'EXCLUDE', 'FIX', 'TRIM'}:
            continue
        marker, taxon = row['marker'], row['taxon_id']
        if marker not in columns or taxon not in original:
            raise ValueError('Flagged observation absent from full matrix universe')
        removed = sum(original[taxon][i] in AA for i in columns[marker])
        observed_removed += removed
        for i in columns[marker]:
            if i in mask_positions[taxon]:
                raise ValueError('Duplicate flagged observation')
            mask_positions[taxon].add(i)
            transformed[taxon][i] = '?'
        flagged.append(dict(row, matrix_columns_masked=len(columns[marker]), canonical_residues_removed=removed))
    a.output.mkdir(parents=True)
    path = a.output / 'matrix.faa'
    path.write_text(''.join('>' + t + '\n' + ''.join(s) + '\n' for t, s in transformed.items()))
    # Independent full readback of all changed and unchanged positions.
    actual_records = list(SeqIO.parse(path, 'fasta'))
    actual = {r.id: str(r.seq) for r in actual_records}
    if len(actual_records) != len(original) or set(actual) != set(original):
        raise ValueError('Output taxon grid changed')
    changed = checked = 0
    for taxon, before in original.items():
        after = actual[taxon]
        if len(after) != len(before):
            raise ValueError('Column count changed')
        for i, (x, y) in enumerate(zip(before, after)):
            expected = '?' if i in mask_positions[taxon] else x
            if y != expected:
                raise ValueError('Unexpected matrix edit')
            changed += x != y
            checked += 1
        if not any(c in AA for c in after):
            raise ValueError('Masking would leave an entirely unobserved taxon')
    # Preserve the original retained columns exactly; no occupancy retrimming.
    (a.output / 'site_mapping.tsv').write_bytes((a.matrix / 'site_mapping.tsv').read_bytes())
    write_table(a.output / 'masked_observations.tsv', flagged)
    summaries = []
    for taxon in original:
        before, after = original[taxon], actual[taxon]
        summaries.append({'taxon_id': taxon, 'columns': len(after),
                          'canonical_before': sum(c in AA for c in before),
                          'canonical_after': sum(c in AA for c in after),
                          'masked_marker_observations': sum(r['taxon_id'] == taxon for r in flagged)})
    write_table(a.output / 'taxon_coverage.tsv', summaries)
    result = {'status': 'complete_full_species_fcs_mask_sensitivity',
              'taxa': len(original), 'columns': baseline['columns'], 'markers': len(columns),
              'masked_marker_taxon_observations': len(flagged),
              'affected_markers': len({r['marker'] for r in flagged}),
              'affected_taxa': sorted(t for t, positions in mask_positions.items() if positions), 'canonical_residues_removed': observed_removed,
              'characters_changed': changed, 'all_matrix_characters_checked': checked,
              'baseline_matrix_receipt_sha256': sha(a.matrix / 'receipt.json'),
              'mapping_receipt_sha256': sha(a.mapping / 'receipt.json'),
              'mapping_audit_receipt_sha256': sha(a.audit / 'receipt.json'),
              'manifest_sha256': baseline['manifest_sha256'], 'script_sha256': sha(Path(__file__)),
              'artifacts': {p.name: sha(p) for p in a.output.iterdir()},
              'interpretation': 'Whole marker/taxon observations overlapping FCS EXCLUDE/FIX/TRIM coding regions masked as unknown. REVIEW-only observations retained. Original taxa and alignment columns preserved, no occupancy retrimming, alternative-copy replacement or confirmed-contamination claim. All resulting taxa retain canonical observations; coupled taxon-identity and heterogeneous-model sensitivities remain necessary.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
