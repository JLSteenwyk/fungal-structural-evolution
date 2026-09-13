#!/usr/bin/env python3
"""Add observed Pfam functional correspondences to the matched site-evolution frame."""
import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import read_table, sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['frame', 'functional', 'projection', 'inputs', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable output')
    receipts = {name: checked_receipt(getattr(a, name))
                for name in ['frame', 'functional', 'projection', 'inputs']}
    if (receipts['frame']['status'] != 'complete_matched_site_evolution_frame'
            or receipts['functional']['status'] != 'complete_functional_site_structure_alignment_join'
            or receipts['frame']['source_receipts']['projection'] != sha(a.projection / 'receipt.json')
            or any(h != sha(a.inputs / 'receipt.json') for h in [
                receipts['frame']['source_receipts']['inputs'],
                receipts['functional']['source_receipts']['paired']['sha256'],
                receipts['projection']['source_receipts']['inputs']])):
        raise ValueError('Mismatched input lineage')
    frame = read_table(a.frame / 'site_evolution_frame.tsv')
    by_site = {(r['marker'], r['paired_column_1based']): r for r in frame}
    if len(by_site) != len(frame):
        raise ValueError('Duplicate frame site')
    annotations = read_table(a.functional / 'site_structure_links.tsv')
    observed = [r for r in annotations if r['observed_in_paired_alignment'] == 'True']
    if len(observed) != receipts['functional']['paired_observed_rows']:
        raise ValueError('Functional observation count differs')
    wanted = {(r['marker'], r['taxon_id'], r['paired_column_1based']) for r in observed}
    coordinates = {}
    with gzip.open(a.projection / 'paired_site_accessibility.tsv.gz', 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            key = r['marker'], r['taxon_id'], r['paired_column_1based']
            if key in wanted:
                if key in coordinates:
                    raise ValueError('Duplicate observed coordinate')
                coordinates[key] = r
    if coordinates.keys() != wanted:
        raise ValueError('Missing functional coordinate observation')
    grouped = defaultdict(list)
    ledger = []
    for r in observed:
        key = r['marker'], r['paired_column_1based']
        c = coordinates[r['marker'], r['taxon_id'], r['paired_column_1based']]
        if key not in by_site or by_site[key]['matrix_column_1based'] != r['matrix_column_1based']:
            raise ValueError('Functional/frame column differs')
        fields = {'protein_id': 'protein_id', 'model_id': 'model_id',
                  'protein_residue_1based': 'protein_residue_1based',
                  'matrix_column_1based': 'matrix_column_1based',
                  'observed_amino_acid': 'amino_acid', 'native_state': '3di_state'}
        if any(r[x] != c[y] for x, y in fields.items()) or r['conserved_candidate'] not in ['True', 'False']:
            raise ValueError('Functional coordinate identity differs')
        grouped[key].append(r)
        ledger.append({**r, 'coordinate_sasa_angstrom_squared': c['sasa_angstrom_squared'],
                       'coordinate_ca_plddt': c['ca_plddt'], 'coordinate_context': c['context']})
    for key, row in by_site.items():
        selected = grouped.get(key, [])
        taxa = {r['taxon_id'] for r in selected}
        candidates = {r['taxon_id'] for r in selected if r['conserved_candidate'] == 'True'}
        other = {r['taxon_id'] for r in selected if r['conserved_candidate'] == 'False'}
        n = int(row['observed_taxa'])
        if len(taxa) > n:
            raise ValueError('More annotated taxa than observed taxa')
        row.update(functional_correspondence_rows=len(selected),
                   functional_correspondence_taxa=len(taxa),
                   conserved_functional_candidate_taxa=len(candidates),
                   other_functional_correspondence_taxa=len(other),
                   functional_correspondence_fraction_observed_taxa=len(taxa) / n,
                   conserved_functional_candidate_fraction_observed_taxa=len(candidates) / n,
                   functional_pfam_accessions=';'.join(sorted({r['pfam_accession'] for r in selected})),
                   functional_annotation_status='observed_correspondence' if selected else 'no_observed_correspondence')
    summary = []
    for marker in sorted({r['marker'] for r in frame}):
        rows = [r for r in frame if r['marker'] == marker]
        summary.append({'marker': marker, 'sites': len(rows),
                        'sites_with_correspondence': sum(r['functional_correspondence_rows'] > 0 for r in rows),
                        'sites_with_conserved_candidate': sum(r['conserved_functional_candidate_taxa'] > 0 for r in rows),
                        'correspondence_rows': sum(r['functional_correspondence_rows'] for r in rows)})
    a.output.mkdir(parents=True)
    for name, rows in [('site_evolution_frame.tsv', frame), ('observed_functional_links.tsv', ledger),
                       ('marker_summary.tsv', summary)]:
        write_table(a.output / name, rows)
    result = {'status': 'complete_functionally_annotated_site_evolution_frame',
              'markers': len(summary), 'sites': len(frame), 'observed_functional_rows': len(ledger),
              'sites_with_correspondence': len(grouped),
              'sites_with_conserved_candidate': sum(r['conserved_functional_candidate_taxa'] > 0 for r in frame),
              'source_site_status_counts': dict(Counter(r['site_status'] for r in annotations)),
              'source_receipts': {name: sha(getattr(a, name) / 'receipt.json') for name in receipts},
              'script_sha256': sha(Path(__file__)),
              'interpretation': 'All original site-estimate rows retained. Counts describe observed projected Pfam correspondences, with conserved candidates distinguished from other correspondences; neither is experimental functional validation. No observed correspondence is unknown annotation status, not evidence of nonfunction. Taxon counts use sets; candidate and other correspondence sets may overlap. Fractions use all observed taxa, not annotation sensitivity. No functional enrichment, independent transition, branch effect or selection test.',
              'artifacts': {f.name: sha(f) for f in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
