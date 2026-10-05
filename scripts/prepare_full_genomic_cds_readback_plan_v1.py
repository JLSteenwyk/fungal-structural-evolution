#!/usr/bin/env python3
"""Prepare the complete independent comparison reader after full source freezing."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from reference_measurement_union_sources import bind, verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['producer-plan','producer-receipt','producer-transport','plan','resources','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    assert not args.plan.exists() and not args.resources.exists() and not args.output.exists()
    original=json.loads(args.producer_plan.read_text());verify(original['pins'])
    assert original['expected_taxa']==len(original['entries'])==526
    assert original['expected_source_products']==5927745 and original['expected_selected_representatives']==5815847
    fixture=Path('metadata/full_genomic_cds_readback_fixture_20261005_v1.json')
    closure=Path('metadata/full_genomic_cds_readback_fixture_transport_20261005_v1.json')
    checked=json.loads(fixture.read_text());transport=json.loads(closure.read_text())
    assert checked['status']=='passed_offline_independent_genomic_CDS_full_source_replay_controls' and checked['count']==18
    assert transport['original_tool_terminal_exit_code']==0;verify(transport['source_hashes'])
    pins={}
    for path in [args.producer_plan,fixture,closure,Path(__file__),Path('scripts/readback_full_genome_annotation_cds_v1.py')]:bind(pins,path)
    now=datetime.now(timezone.utc).isoformat()
    resources=dict(prepared_utc=now,cpus=4,memory_gib=64,swap_gib=0,blas_threads=1,address_space_gib=56,
        cpu_seconds_per_stage=604800,wall_seconds_per_stage=604800,per_file_limit_mib=16384,
        minimum_available_ram_gib=64,minimum_free_disk_gib=2048,output_allowance_gib=8,
        estimated_output_gib=[0.01,1],estimated_wall_hours=[2,72],runtime_is_uncalibrated=True,
        expected_taxa=526,expected_source_products=5927745,expected_representatives=5815847,
        expected_cds_targets=original['expected_cds_targets'],gpu=False,new_predictions=0,new_cost_usd=0,
        basis='Complete original genomes, indexed annotations/CDS references and all available CDS targets; '
              'every coordinate/candidate/target/source product independently replayed. Own FASTA decoder and '
              'byte-table complement, fourCPU64GiB0swap56GiBAS, one BLAS thread. Source-derived planning2-72h '
              'is uncalibrated;7day safety caps are notETA; no full producer repeat or output mutation.')
    plan=dict(prepared_utc=now,producer_plan=str(args.producer_plan),producer_receipt=str(args.producer_receipt),
        producer_transport=str(args.producer_transport),output=str(args.output),cpu=4,expected_taxa=526,
        pins=pins,launch_state='not_launched_or_queued',scientific_eligibility=False,
        gate='Actual original full comparison producer API/native zero, complete original journal and all '
             'source/output bindings; unavailable receipt/transport keeps the full reader unlaunched.')
    for path,value in [(args.resources,resources),(args.plan,plan)]:
        with path.open('x') as handle:json.dump(value,handle,indent=2);handle.write('\n')
    print(json.dumps(dict(taxa=526,source_products=5927745,representatives=5815847,
                         cds_targets=original['expected_cds_targets'],launch_state='not_launched_or_queued'),indent=2))


if __name__=='__main__':main()
