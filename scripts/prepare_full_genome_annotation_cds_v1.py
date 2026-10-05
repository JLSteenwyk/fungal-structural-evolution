#!/usr/bin/env python3
"""Freeze the complete genome/CDS comparison only after both full source readbacks close."""
import argparse
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import shutil

from Bio import SeqIO
import psutil
from reference_measurement_union_sources import bind, verify


def save(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2); handle.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--resources', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    args = parser.parse_args()
    assert not args.plan.exists() and not args.resources.exists() and not args.output_root.exists()
    references = {
        'assembly_readback': 'metadata/full_assembly_dna_readback_20261005_v1.json',
        'assembly_readback_transport': 'metadata/full_assembly_dna_readback_transport_20261005_v1.json',
        'annotation_readback': 'metadata/full_annotation_coordinate_readback_20261005_v1.json',
        'annotation_readback_transport': 'metadata/full_annotation_coordinate_readback_transport_20261005_v1.json'}
    pins = {}
    for name, expected in [('assembly_readback', 'passed_full_selected_genome_disposition_and_original_fasta_reconstruction'),
                           ('annotation_readback', 'passed_full_independent_annotation_coordinate_registry_replay')]:
        path = Path(references[name])
        receipt = json.loads(path.read_text())
        assert receipt['status'] == expected and receipt['taxa'] == 526
        transport_path = Path(references[name+'_transport'])
        transport = json.loads(transport_path.read_text())
        assert transport['status'] == 'verified_original_software_wait_and_whole_wrapper_payloads'
        assert transport['original_tool_terminal_exit_code'] == 0
        for mapping in (receipt['source_hashes'], transport['source_hashes']):
            verify(mapping)
            for file, digest in mapping.items(): bind(pins, file, digest)
        bind(pins,path);bind(pins,transport_path)
    fixture_path = Path('metadata/full_genome_annotation_cds_fixture_20261005_v1.json')
    fixture_transport = Path('metadata/full_genome_annotation_cds_fixture_transport_20261005_v1.json')
    fixture, closure = [json.loads(p.read_text()) for p in (fixture_path, fixture_transport)]
    assert fixture['status'] == 'passed_offline_literal_genome_annotation_cds_controls' and fixture['count'] == 19
    assert closure['original_tool_terminal_exit_code'] == 0
    verify(closure['source_hashes'])
    annotation_plan = Path('metadata/full_annotation_coordinate_plan_20261005_v2.json')
    annotation_receipt = Path('metadata/full_annotation_coordinate_20261005_v2.json')
    assembly_plan = Path('metadata/full_assembly_dna_plan_20261005_v1.json')
    source = json.loads(annotation_plan.read_text())
    reports = {r['taxon_id']: r for r in json.loads(annotation_receipt.read_text())['taxa_reports']}
    genome_plan = json.loads(assembly_plan.read_text())
    genome_table = Path(genome_plan['output']) / 'taxon_dispositions.jsonl'
    genomes = {r['taxon']['taxon_id']: r for r in map(json.loads,genome_table.read_text().splitlines())}
    ncbi_path = Path('metadata/cds_download_receipts.json')
    ncbi_receipt = json.loads(ncbi_path.read_text())
    assert ncbi_receipt['status'] == 'complete_ncbi_cds_acquisition'
    ncbi = {r['taxon_id']: r for r in ncbi_receipt['taxa']}
    external_path = Path('metadata/external_genome_receipts.json')
    external = {r['path']: r for r in json.loads(external_path.read_text())}
    orf_path = Path('metadata/sanchytrid_coordinate_audit.json')
    orfs = {r['taxon_id']:r for r in json.loads(orf_path.read_text())['taxa']}
    published = {'OFS5426470':'data/external/5426470/Clim_long.cds.fasta',
                 'OFS5426494':'data/external/5426494/Nk52_long.cds.fasta',
                 'OFS5426506':'data/external/5426506/Pgem_long.cds.fasta',
                 'OFS5426458':'data/external/5426458/Awhi_long.cds.fasta'}
    for path in (fixture_path,fixture_transport,annotation_plan,annotation_receipt,assembly_plan,
                 genome_table,ncbi_path,external_path,orf_path,Path(__file__),
                 Path('scripts/audit_full_genome_annotation_cds_v1.py'),Path('scripts/genomic_cds_join_v1.py'),
                 Path('scripts/build_full_annotation_coordinate_registry_v2.py'),
                 Path('scripts/check_full_annotation_coordinate_registry_cases_v2.py')):
        bind(pins,path)
    entries=[];target_records=target_bases=0
    for original in source['entries']:
        name=original['taxon_id'];registry=reports[name];genome=genomes[name]
        assert original['taxon']==genome['taxon']
        entry_pins=dict(original['pins'])
        for path,digest in registry['source_hashes'].items():bind(entry_pins,path,digest)
        if genome['status']=='verified_publisher_bound_genomic_dna':
            for key,hash_key in [('path','sha256'),('contig_metadata','contig_metadata_sha256'),
                                  ('receipt_path','receipt_sha256')]:bind(entry_pins,genome[key],genome[hash_key])
        else:assert genome['status']=='genome_acquisition_or_format_error'
        target=None;reason=None
        if name in ncbi:
            receipt=ncbi[name]
            assert receipt['assembly_accession']==original['taxon']['assembly_accession']
            assert receipt['status']=='validated_cds_fasta'
            target=dict(path=receipt['path'],sha256=receipt['sha256'],identifier_mode='ncbi_protein_id_header',
                        expected_records=receipt['cds_records'],source_receipt=receipt)
            target_records+=receipt['cds_records'];target_bases+=receipt['nucleotide_bases']
            bind(entry_pins,receipt['path'],receipt['sha256'])
            bind(entry_pins,receipt['checksum_listing_path'],receipt['checksum_listing_sha256'])
        elif name in published or name in orfs:
            if name in published:
                receipt=external[published[name]];path=receipt['path'];digest=receipt['sha256']
                origin='original_publisher_cds'
            else:
                receipt=orfs[name];path=receipt['cds_path'];digest=receipt['cds_sha256']
                origin='previously_extracted_provisional_orf_cds_not_independent_publisher_target'
            # Inventory every external CDS record; this is not a benchmark/pilot.
            opener=gzip.open if path.endswith('.gz') else open
            count=bases=0
            with opener(path,'rt') as handle:
                for record in SeqIO.parse(handle,'fasta'):count+=1;bases+=len(record.seq)
            target=dict(path=path,sha256=digest,identifier_mode='exact_fasta_id',expected_records=count,
                        source_receipt=receipt,source_kind=origin)
            target_records+=count;target_bases+=bases;bind(entry_pins,path,digest)
        else:
            assert name=='OFS1403592'
            reason='Published mRNA is not an independently qualified CDS target; retain all reconstructed candidates and source products for coding-boundary review'
        for path,digest in entry_pins.items():bind(pins,path,digest)
        entries.append(dict(original,registry_report=registry,genome_report=genome,target_cds=target,
                            missing_target_reason=reason,pins=entry_pins))
    assert len(entries)==len(reports)==len(genomes)==526
    assert sum(e['source_products'] for e in entries)==5927745
    assert sum(e['selected_representatives'] for e in entries)==5815847
    verify(pins)
    now=datetime.now(timezone.utc).isoformat()
    total_genome_bases=sum(r['fasta']['total_length'] for r in genomes.values() if 'fasta' in r)
    largest=max(r['fasta']['total_length'] for r in genomes.values() if 'fasta' in r)
    resources=dict(prepared_utc=now,cpus=4,memory_gib=64,swap_gib=0,blas_threads=1,address_space_gib=56,
        cpu_seconds_per_stage=604800,wall_seconds_per_stage=604800,per_file_limit_mib=16384,
        minimum_available_ram_gib=64,minimum_free_disk_gib=2048,output_allowance_gib=128,
        estimated_output_gib=[2,32],estimated_wall_hours=[2,72],runtime_is_uncalibrated=True,
        source_genome_bases=total_genome_bases,largest_source_genome_bases=largest,
        source_cds_records=target_records,source_cds_bases=target_bases,
        available_ram_gib=psutil.virtual_memory().available/2**30,free_disk_gib=shutil.disk_usage('.').free/2**30,
        gpu=False,new_predictions=0,new_cost_usd=0,
        basis='All526source entries,5,927,745products and60,917,860annotation features. Whole genomes and per-taxon CDS candidates/targets '
              'are loaded by fourCPUworkers;64GiB and56GiBAS cover the largest deposited genome plus annotation/target caches. '
              'Every original coordinate/candidate/target/product disposition retained in gzip tables;2-32GiB estimated output,128GiB allowance. '
              '2-72h uncalibrated planning range;7day CPU/wall caps are notETA. No phase trimming, corpus pilot, GPU or paid resources.')
    plan=dict(prepared_utc=now,**references,output=str(args.output_root),cpu=4,expected_taxa=526,
        expected_source_products=5927745,expected_selected_representatives=5815847,expected_cds_targets=target_records,
        budget=dict(emergency_free_disk_gib=1024),entries=entries,pins=pins,resources=str(args.resources),
        launch_state='not_launched_or_queued',scientific_eligibility=False,
        independent_readback_state='not_prepared_or_launched',
        scope='Full selected assembly/annotation/CDS sources and all products. All exceptional/ambiguous/unsupported joins remain explicit. '
              'Source labels and representatives unchanged. Full independent comparison replay is required before interpretation; '
              'no translation/selection, duplication, contamination, gene-copy or evolutionary acceptance.')
    save(args.resources,resources);save(args.plan,plan)
    print(json.dumps(dict(taxa=526,cds_targets=target_records,cds_bases=target_bases,
                         genome_bases=total_genome_bases,source_bindings=len(pins),launch_state=plan['launch_state']),indent=2))


if __name__=='__main__':main()
