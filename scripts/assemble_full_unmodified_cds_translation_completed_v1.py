#!/usr/bin/env python3
"""Merge closed full producer/reader proof and publish every taxon diagnostic."""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['plan', 'producer', 'producer-transport', 'reader', 'reader-transport', 'table', 'receipt']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.table.exists() and not args.receipt.exists()
    plan, producer, reader, pt, rt = [json.loads(p.read_text()) for p in
        [args.plan, args.producer, args.reader, args.producer_transport, args.reader_transport]]
    assert producer['status'] == 'complete_full_unmodified_cds_fixed_code_translation_pending_independent_readback'
    assert reader['status'] == 'passed_full_unmodified_cds_fixed_code_translation_independent_readback'
    assert pt['original_tool_terminal_exit_code'] == rt['original_tool_terminal_exit_code'] == 0
    assert pt['whole_wrapper_initial_and_terminal_payloads_matched'] and rt['whole_wrapper_initial_and_terminal_payloads_matched']
    assert pt['validation_sha256'] == sha(args.producer) and rt['validation_sha256'] == sha(args.reader)
    assert sha(args.plan) == producer['source_hashes'][str(args.plan)] == reader['source_hashes'][str(args.plan)]
    keys = ['source_products', 'selected_representatives', 'target_records', 'original_dna_bases', 'changed_codon_rows',
            'target_role_status_counts', 'inherited_status_cross_counts', 'changed_codon_role_counts',
            'product_role_status_counts', 'coding_model_translation_cross_counts']
    assert all(producer[key] == reader[key] for key in keys)
    assert producer['taxa'] == reader['taxa'] == plan['expected_taxa'] == 526
    assert producer['source_products'] == 5927745 and producer['selected_representatives'] == 5815847
    assert producer['target_records'] == 5923039
    assert sum(producer['target_role_status_counts'].values()) == 3 * 5923039
    assert sum(producer['product_role_status_counts'].values()) == 3 * 5927745
    assert sum(producer['coding_model_translation_cross_counts'].values()) == 3 * 5927745
    assert sum(producer['changed_codon_role_counts'].values()) == producer['changed_codon_rows']
    pins = {}
    for mapping in [pt['source_hashes'], rt['source_hashes'], producer['source_hashes'], reader['source_hashes'], plan['pins']]:
        for path, value in mapping.items():
            bind(pins, path, value)
    for path in [args.plan, args.producer, args.reader, args.producer_transport, args.reader_transport, Path(__file__)]:
        bind(pins, path)
    verify(pins)
    entries = {e['taxon_id']: e for e in plan['entries']}
    reports = {r['taxon_id']: r for r in producer['taxa_reports']}
    assert len(reports) == 526 and reports.keys() == entries.keys()
    fields = ['taxon_id', 'study_role', 'manifest_species_name', 'canonical_species_taxid', 'identity_review_flags',
        'mapping_mode', 'original_cds_source_kind', 'source_products', 'selected_representatives', 'target_records',
        'original_dna_bases', 'snapshot_nuclear_code', 'snapshot_mitochondrial_code', 'changed_codon_rows',
        'target_role_status_counts_json', 'inherited_status_cross_counts_json', 'changed_codon_role_counts_json',
        'product_role_status_counts_json', 'coding_model_translation_cross_counts_json',
        'scientific_eligibility', 'genetic_code_admission', 'biological_codon_eligibility']
    rows = []
    for taxon in sorted(entries):
        e, report = entries[taxon], reports[taxon]
        assert report['source_products'] == e['source_products']
        assert report['selected_representatives'] == e['selected_representatives']
        assert report['target_records'] == e['target_records']
        row = dict(taxon_id=taxon, study_role=e['study_role'], manifest_species_name=e['species_name'],
            canonical_species_taxid=e['context']['canonical_species_taxid'],
            identity_review_flags=e['context']['identity_review_flags'], mapping_mode=e['mapping_mode'],
            original_cds_source_kind=e['original_cds_source_kind'],
            snapshot_nuclear_code=e['context']['nuclear']['code'], snapshot_mitochondrial_code=e['context']['mitochondrial']['code'],
            scientific_eligibility=False, genetic_code_admission=False, biological_codon_eligibility=False)
        for field in ['source_products', 'selected_representatives', 'target_records', 'original_dna_bases', 'changed_codon_rows']:
            row[field] = report[field]
        for field in keys[5:]:
            row[field + '_json'] = json.dumps(report[field], sort_keys=True, separators=(',', ':'))
        assert row.keys() == set(fields)
        rows.append(row)
    with args.table.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    with args.table.open() as handle:
        saved = list(csv.DictReader(handle, delimiter='\t'))
    assert len(saved) == 526 and len({r['taxon_id'] for r in saved}) == 526
    assert saved == [{key: str(row[key]) for key in fields} for row in rows]
    assert sum(row['study_role'] == 'ingroup' for row in saved) == 501
    assert sum(row['study_role'] == 'outgroup' for row in saved) == 25
    bind(pins, args.table)
    result = dict(status='complete_verified_all526_unmodified_cds_fixed_code_translation',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526,
        **{key: producer[key] for key in keys}, source_bindings=len(pins), source_hashes=pins,
        public_taxon_table=str(args.table), table_readback_complete=True,
        original_transports=[dict(unit=t['unit'], invocation_id=t['invocation_id'],
            original_tool_session_id=t['original_tool_session_id'], actual_terminal_exit_code=t['original_tool_terminal_exit_code'],
            whole_original_payloads_verified=t['whole_wrapper_initial_and_terminal_payloads_matched']) for t in [pt, rt]],
        scientific_eligibility=False, genetic_code_admission=False, biological_codon_eligibility=False,
        best_code_selected=False, dna_modified=False, protein_modified=False,
        gpu=False, new_predictions=0, all_eight_aims_incomplete=True,
        scope='Full original producer and independent reader source/output/execution closure and complete526taxon '
              'diagnostic table. Reuses already independently reconstructed record aggregates, with all declared '
              'hashes rechecked; no third translation parser. Every original/alternative/separate-derived '
              'disposition preserved. Fixed-code differences are not evolutionary events, accepted biological '
              'codes/compartments, codon alignment/selection eligibility or completed biological aims.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
