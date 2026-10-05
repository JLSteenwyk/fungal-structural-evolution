#!/usr/bin/env python3
"""Freeze all 526 annotation/product sources and estimate full registry resources."""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import psutil
from reference_measurement_union_sources import bind, verify


def save(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def main():
    fixture_path = Path('metadata/full_annotation_coordinate_fixture_20261005_v2.json')
    transport_path = Path('metadata/full_annotation_coordinate_fixture_transport_20261005_v2.json')
    fixture, transport = [json.loads(p.read_text()) for p in (fixture_path, transport_path)]
    assert fixture['status'] == 'passed_offline_literal_full_annotation_registry_controls' and fixture['count'] == 16
    assert transport['status'] == 'verified_original_single_call_software_result_and_whole_wrapper_payloads'
    assert transport['original_tool_terminal_exit_code'] == 0
    verify(transport['source_hashes'])
    files = ['metadata/analysis_manifest.tsv', 'metadata/annotation_download_receipts.json',
             'metadata/external_genome_receipts.json', 'metadata/gene_mapping_snapshot.json',
             'metadata/qc_input_receipts.json', 'metadata/gene_representatives_receipt.json',
             'metadata/sanchytrid_coordinate_audit.json', 'metadata/full_assembly_dna_plan_20261005_v1.json']
    with open(files[0]) as handle:
        taxa = list(csv.DictReader(handle, delimiter='\t'))
    def load(name):
        return json.loads(Path('metadata', name + '.json').read_text())
    annotations = {r['taxon_id']: r for r in load('annotation_download_receipts')}
    external = {r['path']: r for r in load('external_genome_receipts')}
    maps = {r['taxon_id']: r for r in load('gene_mapping_snapshot')['taxa']}
    proteins = {r['taxon_id']: r for r in load('qc_input_receipts')}
    representatives = {r['taxon_id']: r for r in load('gene_representatives_receipt')['taxa']}
    orfs = {r['taxon_id']: r for r in load('sanchytrid_coordinate_audit')['taxa']}
    published = {'OFS5426470': 'data/external/5426470/Clim_long.annot.gff',
                 'OFS5426494': 'data/external/5426494/Nk52_long.annot.gff',
                 'OFS5426506': 'data/external/5426506/Pgem_long.annot.gff',
                 'OFS5426458': 'data/external/5426458/Awhi_long.annot.gff',
                 'OFS1403592': 'data/external/1403592/Creolimax_fragrantissima.gtf.gz'}
    identifiers = {r['taxon_id'] for r in taxa}
    assert len(taxa) == len(identifiers) == 526
    assert identifiers == set(maps) == set(proteins) == set(representatives)
    assert set(annotations) | set(published) | set(orfs) == identifiers
    pins = {}
    for path in files + [str(fixture_path), str(transport_path), str(Path(__file__)),
                         'scripts/build_full_annotation_coordinate_registry_v2.py',
                         'scripts/readback_full_annotation_coordinate_registry_v1.py']:
        bind(pins, path)
    entries = []
    for taxon in taxa:
        name = taxon['taxon_id']
        mapping, protein, representative = maps[name], proteins[name], representatives[name]
        if name in annotations:
            annotation = annotations[name]
            assert annotation['status'] == 'validated'
            assert annotation['assembly_accession'] == taxon['assembly_accession']
            assert mapping['mapping_mode'] == 'ncbi_protein_gff'
            annotation_path, annotation_sha = annotation['path'], annotation['sha256']
        elif name in published:
            annotation = external[published[name]]
            annotation_path, annotation_sha = annotation['path'], annotation['sha256']
        else:
            annotation = orfs[name]
            assert mapping['mapping_mode'] == 'verified_orf_coordinates'
            annotation_path, annotation_sha = annotation['mapping_path'], annotation['mapping_sha256']
        assert mapping['annotation_sha256'] == annotation_sha
        assert mapping['proteome_sha256'] == representative['source_sha256'] == protein['sha256']
        assert mapping['mapping_sha256'] == representative['mapping_sha256']
        assert mapping['proteins'] == representative['source_proteins'] == protein['proteins']
        decisions = str(Path(representative['path']).with_suffix('.decisions.tsv'))
        entry_pins = {}
        for path, digest in [(annotation_path, annotation_sha), (protein['input_path'], protein['sha256']),
                             (mapping['mapping_path'], mapping['mapping_sha256']),
                             (representative['path'], representative['sha256']),
                             (decisions, representative['decisions_sha256'])]:
            bind(entry_pins, path, digest)
            bind(pins, path, digest)
        entry = dict(taxon_id=name, taxon=taxon, mapping_mode=mapping['mapping_mode'],
            annotation_path=annotation_path, annotation_receipt=annotation, proteome_path=protein['input_path'],
            mapping_path=mapping['mapping_path'], representative_path=representative['path'], decisions_path=decisions,
            source_products=protein['proteins'], selected_representatives=representative['selected_proteins'], pins=entry_pins)
        if name in annotations:
            entry['expected_feature_counts'] = annotation['feature_counts']
        entries.append(entry)
    verify(pins)
    source_products = sum(r['source_products'] for r in entries)
    selected = sum(r['selected_representatives'] for r in entries)
    assert source_products == 5927745 and selected == 5815847
    now = datetime.now(timezone.utc).isoformat()
    ncbi_rows = sum(sum(r['feature_counts'].values()) for r in annotations.values())
    ncbi_cds = sum(r['feature_counts'].get('CDS', 0) for r in annotations.values())
    compressed = sum(r['compressed_bytes'] for r in annotations.values())
    resource = dict(prepared_utc=now, cpus=4, memory_gib=32, swap_gib=0, blas_threads=1,
        address_space_gib=28, cpu_seconds_per_stage=604800, wall_seconds_per_stage=604800,
        per_file_limit_mib=16384, minimum_available_ram_gib=32, minimum_free_disk_gib=2048,
        output_allowance_gib=256, estimated_output_gib=[20,128], estimated_wall_hours=[1,48],
        runtime_is_uncalibrated=True, ncbi_feature_rows=ncbi_rows, ncbi_cds_segments=ncbi_cds,
        ncbi_annotation_compressed_bytes=compressed,
        largest_ncbi_annotation_feature_rows=max(sum(r['feature_counts'].values()) for r in annotations.values()),
        available_ram_gib=psutil.virtual_memory().available/2**30, free_disk_gib=shutil.disk_usage('.').free/2**30,
        gpu=False, new_predictions=0, new_cost_usd=0,
        basis='All519NCBI annotation receipts contain60,434,823features/23,306,261CDSparts in1,183,724,940compressedbytes. '
              'Includes five external annotations and two full provisional ORF tables,5,927,745sourceproducts and '
              '5,815,847representatives. Four CPU workers,32GiB,no swap, one BLAS thread. SQLite raw lines,decoded attributes '
              'and indexes estimated20-128GiB;256GiB allowance includes journals/temp/index scratch. '
              '1-48h uncalibrated planning range,7day CPU/wall safety caps are not an ETA; no new paid resources.')
    plan_path = 'metadata/full_annotation_coordinate_plan_20261005_v2.json'
    resources_path = 'metadata/full_annotation_coordinate_resources_20261005_v2.json'
    root = 'results/annotation-coordinate-registry-20261005-v2'
    assert not Path(root).exists()
    plan = dict(prepared_utc=now, output=root, cpu=4, expected_taxa=526,
        expected_source_products=source_products, expected_selected_representatives=selected,
        expected_ncbi_feature_rows=ncbi_rows, expected_ncbi_cds_segments=ncbi_cds,
        entries=entries, pins=pins, resources=resources_path,
        budget=dict(emergency_free_disk_gib=1024, output_allowance_gib=192),
        storage_guard_scope='Periodic guard on stable SQLite payload paths at192GiB leaves64GiB declared journal/temp/in-flight scratch. '
                            'Native16GiB per-file cap and1TiB emergency disk reserve also apply.',
        scope='Every selected original annotation line and every normalized source product, unchanged representative baseline; '
              'exact CDS coordinates/strand/phase, attributes/exceptions, Parent and direct candidate product links indexed. '
              'No phase trimming, coordinate normalization, orphan pruning or promotion of provisional ORFs to genes. '
              'Full independent index replay and independently closed genome-DNA/CDS reconstruction remain required; '
              'gene correctness, expression, contamination, duplication and evolutionary events are not established.')
    save(resources_path, resource)
    save(plan_path, plan)
    reader_pins = {}
    for path in (plan_path, 'scripts/readback_full_annotation_coordinate_registry_v1.py', str(fixture_path), str(transport_path)):
        bind(reader_pins, path)
    save('metadata/full_annotation_coordinate_readback_resources_20261005_v1.json', dict(resource,
        output_allowance_gib=8, estimated_output_gib=[0.01,1],
        basis='Independently replay every original annotation line and all5,927,745source/5,815,847selected products against immutable SQLite; '
              'separate percent decoding/FASTA parser; no write to producer databases. FourCPU32GiB0swap28GiBAS, '
              '1-48h uncalibrated planning range and7day hardcaps; no new prediction/paid resources.'))
    save('metadata/full_annotation_coordinate_readback_plan_20261005_v1.json', dict(prepared_utc=now,
        producer_plan=plan_path, producer_receipt='metadata/full_annotation_coordinate_20261005_v2.json',
        producer_transport='metadata/full_annotation_coordinate_transport_20261005_v2.json',
        output='results/annotation-coordinate-independent-readback-20261005-v1', cpu=4,
        expected_taxa=526, launch_state='not_launched_or_queued', pins=reader_pins,
        gate='Actual original full producer API/native zero, whole journal and all declared source/output bindings; '
             'unavailable receipt/transport keeps this reader unlaunched. Genome-CDS work additionally requires qualified assembly DNA.'))
    print(json.dumps(dict(taxa=526, source_products=source_products, selected_representatives=selected,
                         ncbi_feature_rows=ncbi_rows, ncbi_cds_segments=ncbi_cds, complete_source_bindings=len(pins)), indent=2))


if __name__ == '__main__':
    main()
