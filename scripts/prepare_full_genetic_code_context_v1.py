#!/usr/bin/env python3
"""Freeze all526 source-coded products and resolve pinned taxonomy inheritance."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

from genetic_code_context_v1 import ancestry, canonical, inherited_code, read_taxonomy
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--context-output', type=Path, required=True)
    args = p.parse_args()
    assert not args.plan.exists() and not args.context_output.exists()
    pins = {}
    schema_path = Path('metadata/ncbi_genetic_code_schema_source_20261006_v1.json')
    schema = json.loads(schema_path.read_text())
    bind(pins, schema['path'], schema['sha256'])
    identity_path = Path('metadata/selected_taxon_identity_snapshot_completed_20261003_v2.json')
    identity = json.loads(identity_path.read_text())
    assert identity['status'] == 'complete_selected_taxon_frozen_snapshot_identity_audit' and identity['entries'] == 526
    table_path = Path('results/taxonomy/selected-identity-snapshot-20261003-v2/taxon_identity.tsv')
    archive = Path('data/raw/new_taxdump.tar.gz')
    bind(pins, table_path, identity['artifacts'][str(table_path)])
    bind(pins, archive, identity['source_hashes'][str(archive)])
    with table_path.open() as f:
        labels = list(csv.DictReader(f, delimiter='\t'))
    assert len(labels) == len({r['taxon_id'] for r in labels}) == 526
    assert Counter(r['study_role'] for r in labels) == dict(ingroup=501, outgroup=25)
    nodes, merged, codes = read_taxonomy(archive)
    contexts = {}
    for row in labels:
        species = canonical(row['effective_species_taxid'], merged)
        assert species == int(row['canonical_effective_taxid'])
        path = ancestry(species, nodes)
        assert path == [int(x) for x in row['node_path_taxids'].split(';')]
        n, m = [inherited_code(species, nodes, c) for c in ['nuclear', 'mitochondrial']]
        assert n['code'] in codes and m['code'] in codes
        catalog_id = row['catalogue_taxid']
        assembly = None
        if catalog_id:
            canonical_catalog = canonical(catalog_id, merged)
            assembly = dict(original_catalogue_taxid=catalog_id, canonical_taxid=canonical_catalog,
                nuclear=inherited_code(canonical_catalog, nodes, 'nuclear'),
                mitochondrial=inherited_code(canonical_catalog, nodes, 'mitochondrial'))
        contexts[row['taxon_id']] = dict(taxon_id=row['taxon_id'], study_role=row['study_role'],
            frozen_species_name=row['frozen_species_name'], snapshot_species_name=row['snapshot_scientific_name'],
            canonical_species_taxid=species, source_parent_path=path, nuclear=n, mitochondrial=m,
            assembly_catalogue_code_context=assembly, identity_review_flags=row['identity_review_flags'],
            biological_species_delimitation=row['biological_species_delimitation'],
            compartment_assignment='unverified', genetic_code_admission=False)
    completion_path = Path('metadata/full_coding_structure_source_completed_20261005_v1.json')
    completion = json.loads(completion_path.read_text())
    assert completion['status'] == 'complete_verified_full_coding_structure_source_coupling'
    assert all(r['actual_terminal_exit_code'] == 0 for r in completion['original_transports'])
    producer_path = Path('metadata/full_coding_structure_source_coupling_20261005_v1.json')
    bind(pins, producer_path, completion['source_hashes'][str(producer_path)])
    producer = json.loads(producer_path.read_text())
    entries = []
    for report in producer['taxa_reports']:
        target_path, product_path = report['artifact_paths']
        for path in [target_path, product_path]:
            bind(pins, path, completion['source_hashes'][path])
        entries.append(dict(taxon_id=report['taxon_id'], targets=target_path, products=product_path,
            source_products=report['source_products'], selected_representatives=report['selected_representatives'],
            target_records=report['target_records'], context=contexts[report['taxon_id']]))
    assert set(contexts) == {e['taxon_id'] for e in entries}
    assert sum(e['source_products'] for e in entries) == 5927745
    assert sum(e['selected_representatives'] for e in entries) == 5815847
    assert sum(e['target_records'] for e in entries) == 5923039
    software_path = Path('metadata/genetic_code_context_software_validation_20261006_v1.json')
    transport_path = Path('metadata/genetic_code_context_software_transport_20261006_v1.json')
    software, transport = [json.loads(q.read_text()) for q in [software_path, transport_path]]
    assert software['status'] == 'passed_literal_genetic_code_context_inheritance_controls'
    assert transport['original_tool_terminal_exit_code'] == 0
    verify(transport['source_hashes'])
    for path, digest in transport['source_hashes'].items():
        bind(pins, path, digest)
    reader_validation = Path('metadata/genetic_code_context_reader_software_validation_20261006_v2.json')
    reader_transport = Path('metadata/genetic_code_context_reader_software_transport_20261006_v2.json')
    rv, rt = [json.loads(q.read_text()) for q in [reader_validation, reader_transport]]
    assert rv['status'] == 'passed_independent_genetic_code_context_reader_literal_integration'
    assert rt['original_tool_terminal_exit_code'] == 0
    verify(rt['source_hashes'])
    for path, digest in rt['source_hashes'].items():
        bind(pins, path, digest)
    for path in [reader_validation, reader_transport]:
        bind(pins, path)
    for path in [schema_path, identity_path, completion_path, software_path, transport_path,
                 Path('metadata/full_genetic_code_context_resources_20261006_v1.json'),
                 Path('metadata/full_genetic_code_context_readback_resources_20261006_v1.json'),
                 Path(__file__), Path('scripts/genetic_code_context_v1.py'),
                 Path('scripts/build_full_genetic_code_context_v1.py'),
                 Path('scripts/readback_full_genetic_code_context_v1.py')]:
        bind(pins, path)
    verify(pins)
    context_result = dict(status='complete_pinned_all526_taxonomic_code_context',
        checked_utc=datetime.now(timezone.utc).isoformat(), contexts=contexts, codebook=codes,
        selected_nuclear_code_counts=dict(Counter(c['nuclear']['code'] for c in contexts.values())),
        selected_mitochondrial_code_counts=dict(Counter(c['mitochondrial']['code'] for c in contexts.values())),
        source_hashes=dict(pins), scientific_eligibility=False, genetic_code_admission=False,
        scope='Exact historical taxonomy ancestor inheritance for all526canonical selected species records '
              'and521assembly catalogue taxids; all earlier identity/hybrid flags retained. Current primary '
              'readme documents first13field schema only. No current taxonomy, biological code/compartment '
              'qualification, translation, source relabeling or source-target replacement.')
    with args.context_output.open('x') as f:
        json.dump(context_result, f, indent=2)
        f.write('\n')
    bind(pins, args.context_output)
    plan = dict(status='prepared_full526_coding_structure_genetic_code_context_audit',
        prepared_utc=datetime.now(timezone.utc).isoformat(), cpu=4, expected_taxa=526,
        output='results/cds/full-genetic-code-context-20261006-v1',
        reader_output='results/cds/full-genetic-code-context-readback-20261006-v1',
        source_products=5927745, selected_representatives=5815847, target_records=5923039,
        entries=entries, context_receipt=str(args.context_output), taxonomy_archive=str(archive),
        identity_table=str(table_path), pins=pins, source_input_bytes=sum(Path(s).stat().st_size for s in pins),
        resources='metadata/full_genetic_code_context_resources_20261006_v1.json',
        reader_resources='metadata/full_genetic_code_context_readback_resources_20261006_v1.json',
        scientific_eligibility=False, genetic_code_admission=False, biological_codon_eligibility=False,
        scope='All original products/targets and separate derived code tests compared to pinned taxonomy '
              'context without replacing codes or translations. Nuclear versus mitochondrial relationship '
              'does not determine sequence compartment or actual genetic code. Additional new-taxdump '
              'plastid/hydrogenosome fields are not interpreted. Shared source evidence and parser libraries '
              'remain explicit; no codon alignment, selection or evolutionary admission.')
    with args.plan.open('x') as f:
        json.dump(plan, f, indent=2)
        f.write('\n')
    print(json.dumps({k: v for k, v in plan.items() if k not in ['entries', 'pins']}, indent=2))
    print(json.dumps({k: context_result[k] for k in ['selected_nuclear_code_counts', 'selected_mitochondrial_code_counts']}, indent=2))


if __name__ == '__main__':
    main()
