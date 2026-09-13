#!/usr/bin/env python3
"""Measure full-design matrix coverage and explicit taxon-filter sensitivities."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import ROOT, sha, read_table
from assess_pae_sensitivity import checked_receipt


def coverage(path, expected):
    receipt = checked_receipt(path)
    records = list(SeqIO.parse(path / 'matrix.faa', 'fasta'))
    if len(records) != len(expected) or {r.id for r in records} != expected:
        raise ValueError('Matrix taxon identity/count differs')
    counts = {}
    alphabet = set('ACDEFGHIKLMNPQRSTVWY')
    for record in records:
        sequence = str(record.seq)
        if len(sequence) != receipt['columns'] or not set(sequence) <= alphabet | set('X-?'):
            raise ValueError('Unexpected matrix length/alphabet')
        counts[record.id] = sum(sequence.count(aa) for aa in alphabet)
    rows = read_table(path / 'taxon_coverage.tsv')
    if len(rows) != len(expected) or {r['taxon_id'] for r in rows} != expected:
        raise ValueError('Coverage table scope differs')
    for row in rows:
        observed = counts[row['taxon_id']]
        if observed != int(row['unambiguous_residues']) or abs(observed/receipt['columns']-float(row['unambiguous_fraction'])) > 1e-12:
            raise ValueError('Archived matrix coverage differs from direct count')
    return receipt, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=Path, required=True)
    parser.add_argument('--mafft', type=Path, required=True)
    parser.add_argument('--copies', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError('Use a new immutable sensitivity output')
    manifest_path = ROOT / 'metadata/analysis_manifest.tsv'
    manifest = read_table(manifest_path)
    expected = {r['taxon_id'] for r in manifest}
    if len(expected) != len(manifest): raise ValueError('Repeated manifest taxa')
    profile, pc = coverage(args.profile, expected)
    mafft, mc = coverage(args.mafft, expected)
    copy = checked_receipt(args.copies)
    if any(r['manifest_sha256'] != sha(manifest_path) for r in [profile, mafft, copy]):
        raise ValueError('Different sampling provenance')
    cr = read_table(args.copies / 'taxon_gene_copy_summary.tsv')
    copies = {r['taxon_id']: r for r in cr}
    if set(copies) != expected or len(copies) != len(cr): raise ValueError('Different copy audit scope')
    rows, membership, groups = [], [], defaultdict(list)
    for taxon in manifest:
        name = taxon['taxon_id']
        p, m = pc[name]/profile['columns'], mc[name]/mafft['columns']
        row = {'taxon_id': name, 'species_name': taxon['species_name'], 'study_role': taxon['study_role'],
            'manifest_lineage_group': taxon['lineage'].split(';')[0], 'profile_unambiguous_residues': pc[name],
            'profile_fraction': p, 'mafft_unambiguous_residues': mc[name], 'mafft_fraction': m,
            'minimum_fraction_across_methods': min(p,m),
            'raw_single_copy_markers': int(copies[name]['raw_complete_single']),
            'raw_duplicated_markers': int(copies[name]['raw_duplicated']),
            'raw_fragmented_markers': int(copies[name]['raw_fragmented']), 'raw_missing_markers': int(copies[name]['raw_missing'])}
        rows.append(row); groups[(row['study_role'], row['manifest_lineage_group'])].append(row)
        for threshold in [0.1, 0.3, 0.5, 0.7]:
            membership.append({'taxon_id': name, 'minimum_fraction_threshold': threshold,
                'passes_profile': p >= threshold, 'passes_mafft': m >= threshold, 'passes_both': min(p,m) >= threshold})
    sensitivity = []
    for (role, group), subset in sorted(groups.items()):
        for threshold in [0.1, 0.3, 0.5, 0.7]:
            kept = sum(r['minimum_fraction_across_methods'] >= threshold for r in subset)
            sensitivity.append({'study_role': role, 'manifest_lineage_group': group, 'minimum_fraction_threshold': threshold,
                'original_taxa': len(subset), 'taxa_passing_both': kept, 'taxa_excluded': len(subset)-kept,
                'entire_manifest_group_lost': kept == 0})
    args.output.mkdir(parents=True)
    for name, table in [('taxon_matrix_coverage.tsv', rows), ('coverage_sensitivity_membership.tsv', membership), ('lineage_sensitivity.tsv', sensitivity)]:
        with (args.output / name).open('w') as handle:
            writer = csv.DictWriter(handle, list(table[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(table)
    totals = []
    for threshold in [0.1, 0.3, 0.5, 0.7]:
        kept = [r for r in rows if r['minimum_fraction_across_methods'] >= threshold]
        totals.append({'minimum_fraction_threshold': threshold, 'taxa_retained': len(kept),
            'roles_retained': dict(Counter(r['study_role'] for r in kept)),
            'entire_groups_lost': [r['manifest_lineage_group'] for r in sensitivity if r['minimum_fraction_threshold']==threshold and r['entire_manifest_group_lost']]})
    result = {'status': 'complete_full_taxon_matrix_coverage_sensitivity', 'taxa': len(rows),
        'profile_columns': profile['columns'], 'mafft_columns': mafft['columns'], 'threshold_summaries': totals,
        'profile_receipt_sha256': sha(args.profile/'receipt.json'), 'mafft_receipt_sha256': sha(args.mafft/'receipt.json'),
        'copy_audit_receipt_sha256': sha(args.copies/'receipt.json'), 'script_sha256': sha(Path(__file__)),
        'interpretation': 'Direct canonical-residue counts verified in both complete matrices. Thresholds 10/30/50/70 percent are descriptive taxon-sensitivity definitions, not calibrated reliability thresholds or final exclusion decisions. Different alignment scopes have different denominators. No production taxon is removed, no filtered tree inferred, and passing coverage does not establish correct placement. Lost groups are descriptive manifest bins, not equal taxonomic ranks.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__': main()
