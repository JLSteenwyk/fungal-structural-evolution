#!/usr/bin/env python3
"""Freeze every source product/target and inherited coding/model evidence before joining."""
from collections import Counter
from datetime import datetime, timezone
import csv
import json
from pathlib import Path

import psutil

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def main():
    plan_path = Path('metadata/full_coding_structure_source_coupling_plan_20261005_v1.json')
    resources_path = Path('metadata/full_coding_structure_source_coupling_resources_20261005_v1.json')
    reader_resources_path = Path('metadata/full_coding_structure_source_coupling_readback_resources_20261005_v1.json')
    assert not any(p.exists() for p in [plan_path, resources_path, reader_resources_path])
    pins = {}
    for stem, status in [('coding_structure_source_coupling', 'passed_literal_coding_structure_source_coupling_contracts'),
                         ('coding_structure_source_reader', 'passed_separate_coding_structure_full_reader_literal_source_integration')]:
        receipt_path = Path('metadata', stem + '_software_validation_20261005_v1.json')
        transport_path = Path('metadata', stem + '_software_transport_20261005_v1.json')
        receipt, transport = read(receipt_path), read(transport_path)
        assert receipt['status'] == status and transport['original_tool_terminal_exit_code'] == 0
        assert transport['validation_sha256'] == sha(receipt_path)
        assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
        verify(transport['source_hashes'])
        for path, digest in transport['source_hashes'].items():
            bind(pins, path, digest)
        for path in [receipt_path, transport_path]:
            bind(pins, path)
    genome_path = Path('metadata/full_genomic_cds_completed_20261005_v1.json')
    atlas_path = Path('metadata/full_prediction_atlas_union_completed_20261005_v1.json')
    genome, atlas = read(genome_path), read(atlas_path)
    assert genome['status'] == 'complete_verified_full_genome_annotation_cds_comparison'
    assert atlas['status'] == 'complete_verified_full_prediction_atlas_availability_union'
    assert genome['taxa'] == atlas['taxa'] == 526
    assert all(t['actual_terminal_exit_code'] == 0 for d in [genome, atlas] for t in d['original_transports'])
    genome_receipt_path = Path('metadata/full_genome_annotation_cds_20261005_v1.json')
    assert genome['source_hashes'][str(genome_receipt_path)] == sha(genome_receipt_path)
    reports = read(genome_receipt_path)['taxa_reports']
    manifest_path = Path('metadata/analysis_manifest.tsv')
    with manifest_path.open() as handle:
        manifest = {r['taxon_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    assert len(manifest) == len(reports) == 526
    assert Counter(r['study_role'] for r in manifest.values()) == dict(ingroup=501, outgroup=25)
    ncbi_root = Path('results/cds/ncbi-strict-translation-v1')
    ncbi = read(ncbi_root / 'receipt.json')
    assert ncbi == read('metadata/ncbi_cds_translation_receipt.json')
    assert ncbi['status'] == 'complete_full_ncbi_strict_cds_translation_audit' and ncbi['taxa_count'] == 519
    assert sha(ncbi_root / 'config.json') == ncbi['config_sha256']
    ncbi_reader_path = Path('metadata/ncbi_cds_translation_readback.json')
    ncbi_reader = read(ncbi_reader_path)
    assert ncbi_reader['status'] == 'complete_full_ncbi_cds_audit_readback'
    assert ncbi_reader['source_receipt_sha256'] == sha(ncbi_root / 'receipt.json')
    by_taxon = {r['taxon_id']: r for r in ncbi['taxa']}
    assert len(by_taxon) == 519
    external_roots = {
        'strict': Path('results/cds/published-outgroup-translation-v1'),
        'projection': Path('results/cds/published-outgroup-codon-projection-v2'),
        'creolimax': Path('results/cds/creolimax-extraction-v1')}
    external = {k: read(root / 'receipt.json') for k, root in external_roots.items()}
    assert external['strict']['status'] == 'complete_published_outgroup_cds_translation_audit'
    assert external['projection']['status'] == 'complete_published_outgroup_genome_codon_projection_audit'
    assert external['creolimax']['status'] == 'complete_creolimax_cds_extraction_audit'
    assert external['projection']['strict_receipt_sha256'] == sha(external_roots['strict'] / 'receipt.json')
    projection_reader_path = Path('metadata/published_outgroup_codon_readback.json')
    assert read(projection_reader_path)['receipt_sha256'] == sha(external_roots['projection'] / 'receipt.json')
    entries = []
    for report in sorted(reports, key=lambda r: r['taxon_id']):
        taxon = report['taxon_id']
        label = manifest[taxon]
        entry = dict(taxon_id=taxon, mapping_mode=report['mapping_mode'], study_role=label['study_role'],
            species_name=label['species_name'], lineage=label['lineage'], source_products=report['source_products'],
            selected_representatives=report['selected_representatives'], target_records=report['target_records'],
            genomic_targets=report['artifact_paths'][2], genomic_products=report['artifact_paths'][3])
        for path in [entry['genomic_targets'], entry['genomic_products']]:
            bind(pins, path, report['source_hashes'][path])
        if entry['mapping_mode'] == 'ncbi_protein_gff':
            old = by_taxon[taxon]
            old_receipt = ncbi_root / (taxon + '.receipt.json')
            assert read(old_receipt) == old
            assert old['cds_records'] == entry['target_records']
            assert old['source_hashes']['cds'] == report['target_cds_source']['sha256']
            entry['translation_table'] = str(ncbi_root / (taxon + '.audit.tsv'))
            bind(pins, entry['translation_table'], old['artifacts'][taxon + '.audit.tsv'])
            bind(pins, old_receipt)
        elif entry['mapping_mode'] == 'transcript_gff':
            entry['translation_table'] = str(external_roots['strict'] / 'protein_cds_audit.tsv')
            entry['supplementary_table'] = str(external_roots['projection'] / 'codon_projection_audit.tsv')
        elif entry['mapping_mode'] == 'verified_orf_coordinates':
            entry['translation_table'] = str(Path('results/sanchytrid_coordinates') / (taxon + '.tsv'))
            old = next(r for r in read('metadata/sanchytrid_coordinate_audit.json')['taxa'] if r['taxon_id'] == taxon)
            bind(pins, entry['translation_table'], old['mapping_sha256'])
            assert report['target_cds_source']['sha256'] == old['cds_sha256']
        else:
            assert entry['mapping_mode'] == 'creolimax_gtf'
            entry['supplementary_table'] = str(external_roots['creolimax'] / 'protein_cds_audit.tsv')
            assert report['target_cds_source'] is None
        entries.append(entry)
    for key, root in external_roots.items():
        bind(pins, root / 'receipt.json')
        for name in ['protein_cds_audit.tsv', 'codon_projection_audit.tsv']:
            if name in external[key]['artifacts']:
                bind(pins, root / name, external[key]['artifacts'][name])
    atlas_root = Path('results/structures/full-prediction-atlas-union-20261005-v1')
    database, availability_table = atlas_root / 'prediction_atlas.sqlite', atlas_root / 'protein_source_links.tsv'
    for path in [database, availability_table]:
        bind(pins, path, atlas['source_hashes'][str(path)])
    for path in [genome_path, atlas_path, genome_receipt_path, manifest_path, ncbi_root / 'receipt.json',
                 ncbi_root / 'config.json', ncbi_reader_path, projection_reader_path,
                 Path('metadata/sanchytrid_coordinate_audit.json'), Path(__file__),
                 *[Path('scripts', name + '.py') for name in ['coding_structure_source_coupling_v1',
                     'build_full_coding_structure_source_coupling_v1', 'readback_full_coding_structure_source_coupling_v1']]]:
        bind(pins, path)
    assert sum(e['source_products'] for e in entries) == 5927745
    assert sum(e['selected_representatives'] for e in entries) == 5815847
    assert sum(e['target_records'] for e in entries) == 5923039
    resources = dict(prepared_utc=datetime.now(timezone.utc).isoformat(), cpus=4, memory_gib=64,
        swap_gib=0, blas_threads=1, address_space_gib=56, cpu_seconds_per_stage=604800,
        wall_seconds_per_stage=604800, per_file_limit_mib=16384, minimum_available_ram_gib=64,
        minimum_free_disk_gib=256, output_allowance_gib=32, estimated_output_gib=[0.5,16],
        estimated_wall_hours=[0.5,12], runtime_is_uncalibrated=True, gpu=False, new_predictions=0, new_cost_usd=0,
        qualified_join_input_bytes=sum(Path(path).stat().st_size for path in pins),
        available_ram_gib=psutil.virtual_memory().available / 2**30, free_disk_gib=psutil.disk_usage('.').free / 2**30,
        scope='Complete526entry original coding/product/model source joins and gzip record exports under4CPU64GiB0swap56GiBAS. '
              'Four per-taxon dictionaries; model SQLite views, inherited translation TSVs and validated genomic source outputs. '
              'Planning0.5-12h/0.5-16GiB is uncalibrated;7day safety caps are notETA. No raw genome replay, new translation, '
              'GPU, model inference, paid resource or biological codon admission.')
    assert resources['available_ram_gib'] >= 64 and resources['free_disk_gib'] >= 256
    reader_resources = dict(resources, output_allowance_gib=8, estimated_output_gib=[0.01,2],
        scope='Full526entry/5927745product/5923039target independent source-record join replay; complete model TSV '
              'byte indexing and per-taxon decoding independent of producer SQLite. Inherited translation is not '
              'independently recomputed.4CPU64GiB0swap56GiBAS; no new biological/model inference.')
    save(resources_path, resources)
    save(reader_resources_path, reader_resources)
    for path in [resources_path, reader_resources_path]:
        bind(pins, path)
    verify(pins)
    plan = dict(status='prepared_full526_coding_structure_source_coupling_after_genomic_and_model_union_closure',
        prepared_utc=resources['prepared_utc'], cpu=4, output='results/cds/full-coding-structure-source-coupling-20261005-v1',
        reader_output='results/cds/full-coding-structure-source-readback-20261005-v1', expected_taxa=526,
        source_products=5927745, selected_representatives=5815847, target_records=5923039, entries=entries,
        availability_database=str(database), availability_table=str(availability_table),
        resources=str(resources_path), reader_resources=str(reader_resources_path), pins=pins,
        producer_controls=17, reader_semantic_controls=3,
        control_scope='Reader control invokes semantic replay directly, without an artifact-hash gate; '
                      'legacy control labels mention self-consistent hashes but no hash-rebinding exercise was performed.',
        inherited_translation_scope='Existing audited statuses and recorded code policy retained; no new independent '
                                    'nonmarker translation, coding-boundary correction or scientific eligibility claim.',
        scope='Full original targets/products, all501fungalentries+25outgroups, all alternatives and source missingness. '
              'Separate annotation-derived projections/extractions remain additional dependent evidence, never substitutes '
              'for original CDS/genomic target outcomes. Representative structure availability is before confidence and '
              'predictor accuracy qualification. Taxonomic labels, gene/isoform uncertainty, alignment/divergence, '
              'homology/reconciliation/phylogeny and all8evolutionary aims remain unqualified or incomplete.')
    assert not Path(plan['output']).exists() and not Path(plan['reader_output']).exists()
    save(plan_path, plan)
    print(json.dumps({k: v for k, v in plan.items() if k not in ['pins', 'entries']}, indent=2))


if __name__ == '__main__':
    main()
