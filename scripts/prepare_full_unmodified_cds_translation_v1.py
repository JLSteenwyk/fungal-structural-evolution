#!/usr/bin/env python3
"""Freeze every original CDS/protein and qualified fixed-code diagnostic dependency."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import Bio.SeqIO.FastaIO

from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    assert not args.plan.exists()
    pins = {}
    source_completion = Path('metadata/full_coding_structure_source_completed_20261005_v1.json')
    code_completion = Path('metadata/full_genetic_code_context_completed_20261006_v1.json')
    source_closed, code_closed = [json.loads(p.read_text()) for p in [source_completion, code_completion]]
    assert source_closed['status'] == 'complete_verified_full_coding_structure_source_coupling'
    assert code_closed['status'] == 'complete_verified_all526_coding_and_taxonomic_genetic_code_context'
    for closed in [source_closed, code_closed]:
        assert all(t['actual_terminal_exit_code'] == 0 for t in closed['original_transports'])
    context_path = Path('metadata/full_taxonomic_genetic_code_context_20261006_v1.json')
    bind(pins, context_path, code_closed['source_hashes'][str(context_path)])
    context = json.loads(context_path.read_text())
    assert context['status'] == 'complete_pinned_all526_taxonomic_code_context'
    genomic_plan = Path('metadata/full_genome_annotation_cds_plan_20261005_v1.json')
    genomic = json.loads(genomic_plan.read_text())
    sources = {e['taxon_id']: e for e in genomic['entries']}
    coding_producer = Path('metadata/full_coding_structure_source_coupling_20261005_v1.json')
    bind(pins, coding_producer, source_closed['source_hashes'][str(coding_producer)])
    coding = json.loads(coding_producer.read_text())
    entries = []
    for report in coding['taxa_reports']:
        taxon = report['taxon_id']
        original = sources[taxon]
        targets, products = report['artifact_paths']
        for path in [targets, products]:
            bind(pins, path, source_closed['source_hashes'][path])
        proteome = original['proteome_path']
        bind(pins, proteome, original['pins'][proteome])
        original_cds = original['target_cds']
        if original_cds is not None:
            bind(pins, original_cds['path'], original_cds['sha256'])
            assert original_cds['expected_records'] == report['target_records']
        else:
            assert report['target_records'] == 0 and original['mapping_mode'] == 'creolimax_gtf'
        assert report['source_products'] == original['source_products']
        assert report['selected_representatives'] == original['selected_representatives']
        entries.append(dict(taxon_id=taxon, mapping_mode=original['mapping_mode'], study_role=report['study_role'],
            species_name=report['manifest_species_name'], proteome=proteome,
            original_cds=original_cds['path'] if original_cds is not None else None,
            original_cds_source_kind=original_cds.get('source_kind', original['mapping_mode'])
                if original_cds is not None else 'no_qualified_original_target',
            targets=targets, products=products, source_products=report['source_products'],
            selected_representatives=report['selected_representatives'], target_records=report['target_records'],
            context=context['contexts'][taxon]))
    assert len(entries) == len({e['taxon_id'] for e in entries}) == 526
    assert {e['taxon_id'] for e in entries} == set(sources) == set(context['contexts'])
    assert Counter(e['study_role'] for e in entries) == dict(ingroup=501, outgroup=25)
    assert sum(e['source_products'] for e in entries) == 5927745
    assert sum(e['selected_representatives'] for e in entries) == 5815847
    assert sum(e['target_records'] for e in entries) == 5923039
    qualifications = [
        ('metadata/unmodified_cds_translation_software_validation_20261006_v1.json',
         'metadata/unmodified_cds_translation_software_transport_20261006_v1.json',
         'passed_exhaustive_fixed_code_translation_and_unmodified_source_controls'),
        ('metadata/full_unmodified_cds_integration_validation_20261006_v1.json',
         'metadata/full_unmodified_cds_integration_transport_20261006_v1.json',
         'passed_literal_unmodified_cds_writer_independent_reader_integration')]
    for validation_path, transport_path, status in qualifications:
        validation, transport = [json.loads(Path(p).read_text()) for p in [validation_path, transport_path]]
        assert validation['status'] == status and transport['original_tool_terminal_exit_code'] == 0
        verify(transport['source_hashes'])
        for path, sha in transport['source_hashes'].items():
            bind(pins, path, sha)
        for path in [validation_path, transport_path]:
            bind(pins, path)
    for path in [source_completion, code_completion, genomic_plan, Path(__file__),
        Path('scripts/full_unmodified_cds_translation_v1.py'), Path('scripts/build_full_unmodified_cds_translation_v1.py'),
        Path('scripts/readback_full_unmodified_cds_translation_v1.py'), Path('scripts/unmodified_cds_translation_v1.py'),
        Path('scripts/independent_codon_translation_v1.py'), Path(Bio.SeqIO.FastaIO.__file__),
        Path('metadata/full_unmodified_cds_translation_resources_20261006_v1.json'),
        Path('metadata/full_unmodified_cds_translation_readback_resources_20261006_v1.json')]:
        bind(pins, path)
    verify(pins)
    result = dict(status='prepared_full526_unmodified_original_cds_fixed_code_translation',
        prepared_utc=datetime.now(timezone.utc).isoformat(), cpu=4, expected_taxa=526,
        source_products=5927745, selected_representatives=5815847, target_records=5923039,
        entries=entries, context_receipt=str(context_path), pins=pins,
        source_input_bytes=sum(Path(p).stat().st_size for p in pins),
        output='results/cds/full-unmodified-cds-translation-20261006-v1',
        reader_output='results/cds/full-unmodified-cds-translation-readback-20261006-v1',
        resources='metadata/full_unmodified_cds_translation_resources_20261006_v1.json',
        reader_resources='metadata/full_unmodified_cds_translation_readback_resources_20261006_v1.json',
        planned_wall_hours_uncalibrated=[0.5, 24], planned_output_gib_uncalibrated=[0.5, 32], new_cost_usd=0,
        fixed_possible_positive_codes=[1, 3, 4, 5, 6, 12, 16, 26], unspecified_code_zero_tested=False,
        scientific_eligibility=False, genetic_code_admission=False, biological_codon_eligibility=False,
        reader_launch_gate='Actual original producer terminal zero, complete receipt, immutable hashes and whole original wrapper/invocation transport closure; no reader is launched or queued by preparation.',
        scope='Full original source FASTA/target/product/protein-order/hash and fixed code diagnostics. '
              '519 NCBI, four original publisher CDS and two dependent original ORF target sources; '
              'Creolimax has no qualified original target and all8694products retained. Preserve all '
              '111898alternatives, no-target and separate-derived cases. No original DNA/protein edits, '
              'frame/phase/initiation/boundary correction, best-code search, code/compartment adoption, '
              'family codon alignment or selection. Taxonomic code snapshot and source library dependencies '
              'are explicit; uncertain biological species identities and all eight aims remain unresolved.')
    with args.plan.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['pins', 'entries']}, indent=2))


if __name__ == '__main__':
    main()
