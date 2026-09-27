#!/usr/bin/env python3
"""Freeze the full local-alignment tree diagnostic grid and resource plan."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from Bio import SeqIO


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    with Path(p).open() as h:
        rows = list(csv.DictReader(h, delimiter='\t'))
    assert len({r['case_id'] for r in rows}) == len(rows)
    return {r['case_id']: r for r in rows}


def main():
    root = Path('results/cds/full-group-alignment-readiness-20260927-v1')
    receipt = json.loads((root / 'receipt.json').read_text())
    proof_path = Path('metadata/full_codon_alignment_integration_readback_20260927.json')
    proof = json.loads(proof_path.read_text())
    assert proof['status'] == 'passed_full_codon_alignment_integration_readback'
    assert proof['source_receipt_sha256'] == sha(root / 'receipt.json')
    assert receipt['artifacts']['cases.tsv'] == sha(root / 'cases.tsv')
    rows = read(root / 'cases.tsv')
    fcs_path = Path('metadata/fcs_codon_omission_disposition.tsv')
    fcs = read(fcs_path)
    historical = json.loads(Path('metadata/codon_analysis_readiness_receipt.json').read_text())
    assert sha(fcs_path) == historical['sources'][str(fcs_path)]
    assert set(fcs) <= set(rows) and len(rows) == 1712
    mapping_root = Path('results/qc/fcs-cds-overlap-v1')
    mapping_path = mapping_root / 'marker_overlap_review.tsv'
    mapping_receipt = json.loads((mapping_root / 'receipt.json').read_text())
    mapping_proof_path = Path('metadata/fcs_cds_overlap_audit_receipt.json')
    mapping_proof = json.loads(mapping_proof_path.read_text())
    assert mapping_proof['status'] == 'passed_full_region_cds_intersection_and_protein_marker_join_audit'
    assert mapping_proof['mapping_receipt_sha256'] == sha(mapping_root / 'receipt.json')
    assert sha(mapping_path) == mapping_receipt['artifacts'][mapping_path.name]
    with mapping_path.open() as h:
        flagged = {(r['marker'], r['taxon_id']) for r in csv.DictReader(h, delimiter='\t')
                   if {'EXCLUDE', 'FIX', 'TRIM'} & set(r['fcs_actions'].split(';'))}
    inputs = Path('results/cds/full-group-codon-realignment-projection-20260927-v1')
    input_receipt = json.loads((inputs / 'receipt.json').read_text())
    output = []
    for case, r in rows.items():
        path = inputs / case / 'codons.fna'
        assert sha(path) == input_receipt['artifacts'][f'{case}/codons.fna']
        taxa = {record.id for record in SeqIO.parse(path, 'fasta')}
        assert len(taxa) == int(r['local_taxa'])
        marker = case.split('__')[1]
        omitted = {taxon for taxon in taxa if (marker, taxon) in flagged}
        f = dict(omitted_taxa=';'.join(sorted(omitted)), remaining_taxa=len(taxa - omitted),
                 status='below_existing_four_taxon_minimum' if omitted else 'unchanged_by_declared_fcs_omission')
        if case in fcs:
            assert int(fcs[case]['original_taxa']) == len(taxa)
            assert fcs[case]['omitted_taxa'] == f['omitted_taxa']
            assert int(fcs[case]['remaining_taxa']) == f['remaining_taxa']
            assert fcs[case]['status'] == f['status'] == r['fcs_omission_status']
        information_passes = r['local_information_screen_passes'] == 'True'
        # No source-flagged taxon is silently retained or replaced. All affected
        # groups currently have three taxa after omission and stay in the ledger.
        if f['omitted_taxa']:
            assert int(f['remaining_taxa']) < 4
            status = 'fcs_omission_below_four_taxa'
        elif information_passes:
            status = 'ready_for_supported_tree_diagnostic'
        else:
            status = 'insufficient_sequence_information_for_this_tree_workflow'
        output.append(dict(case_id=case, taxa=r['local_taxa'],
                           nucleotide_columns=r['local_nucleotide_columns'],
                           distinct_aligned_sequences=r['local_distinct_aligned_sequences'],
                           parsimony_informative_nucleotide_columns=r['local_parsimony_informative_nucleotide_columns'],
                           marker_copy_caveat=r['copy_caveat'],
                           historical_review_flags=r['case_specific_review_flags'],
                           local_information_screen_passes=r['local_information_screen_passes'],
                           fcs_omitted_taxa=f['omitted_taxa'],
                           fcs_remaining_taxa=f['remaining_taxa'],
                           fcs_omission_status=f['status'], status=status))
    info = Path('metadata/local_codon_tree_information_20260927.tsv')
    with info.open('x') as h:
        writer = csv.DictWriter(h, fieldnames=list(output[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(output)
    counts = dict(Counter(r['status'] for r in output))
    pins = [root / 'receipt.json', root / 'cases.tsv', proof_path, fcs_path,
            mapping_path, mapping_root / 'receipt.json', mapping_proof_path,
            Path('metadata/codon_analysis_readiness_receipt.json'),
            Path('scripts/run_local_codon_trees.py'),
            Path('scripts/run_paired_marker_fits.py'),
            Path('scripts/assess_pae_sensitivity.py'), Path('scripts/audit_busco_gene_copies.py'),
            Path('scripts/audit_genus_codon_trees.py'), Path(__file__)]
    plan = dict(status='prelaunch_full_local_codon_tree_diagnostics', cases_screened=len(rows),
                ready_cases=counts['ready_for_supported_tree_diagnostic'], dispositions=counts,
                input_receipt_sha256=sha(inputs / 'receipt.json'), information_table_sha256=sha(info),
                workers=4, threads_per_job=1, memory_per_job_gb=2, service_memory_gib=16,
                output_allowance_gb=10, planning_hours=[1, 96],
                model='GTR+F+G4', sh_alrt=1000, ultrafast_bootstrap=1000, bootstrap_nni=True,
                pins={str(p): sha(p) for p in pins},
                scope='Full local-alignment diagnostic set after existing information screen and source-contamination omission policy. Copy caveats retained for diagnostic sensitivity only; no selection eligibility. Same seeds/model/support design as original fits. All 1712 cases remain in the disposition table.',
                cost='Existing authorized local CPU host; no GPU work or new charges')
    target = Path('metadata/local_codon_tree_resource_plan_20260927.json')
    with target.open('x') as h:
        json.dump(plan, h, indent=2); h.write('\n')
    print(json.dumps(counts))


if __name__ == '__main__':
    main()
