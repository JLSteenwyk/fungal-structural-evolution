#!/usr/bin/env python3
"""Replay every original code-context record with independent taxonomy traversal."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import re
import tarfile
import time

from reference_measurement_union_sources import bind, verify


def pairs(items):
    result = {}
    for key, value in items:
        assert key not in result, 'duplicate JSON key'
        result[key] = value
    return result


def reject_constant(value):
    raise ValueError('nonstandard JSON constant: ' + value)


def loads(text):
    return json.loads(text, object_pairs_hook=pairs, parse_constant=reject_constant)


def source_contexts(archive, identity_table):
    nodes, aliases, codes = {}, {}, {}
    with tarfile.open(archive) as tar:
        for member in ['nodes.dmp', 'merged.dmp', 'gencode.dmp']:
            with tar.extractfile(member) as stream:
                for blob in stream:
                    row = re.split(r'\s*\|\s*', blob.decode().rstrip('\n'))
                    assert row.pop() == ''
                    if member == 'nodes.dmp':
                        assert len(row) >= 13 and row[7] in ('0', '1') and row[9] in ('0', '1')
                        taxid = int(row[0])
                        assert taxid not in nodes
                        nodes[taxid] = dict(parent=int(row[1]), nuclear=int(row[6]), nuclear_inherited=int(row[7]),
                                           mitochondrial=int(row[8]), mitochondrial_inherited=int(row[9]))
                    elif member == 'merged.dmp':
                        assert len(row) == 2 and int(row[0]) not in aliases
                        aliases[int(row[0])] = int(row[1])
                    else:
                        assert len(row) == 5 and row[0] not in codes
                        codes[row[0]] = dict(abbreviation=row[1], name=row[2], amino_acids=row[3], starts=row[4])
    def resolve_id(original):
        node, seen = int(original), set()
        while node in aliases:
            assert node not in seen
            seen.add(node)
            node = aliases[node]
        return node
    def path(start):
        result = [start]
        while result[-1] != 1:
            parent = nodes[result[-1]]['parent']
            assert parent not in result and parent in nodes
            result.append(parent)
        assert nodes[1]['parent'] == 1
        return result
    def code(start, name):
        trace, seen = [], set()
        while True:
            assert start in nodes and start not in seen
            seen.add(start)
            value = nodes[start]
            flag, number = value[name + '_inherited'], value[name]
            trace.append(dict(taxid=start, recorded_code=number, inherited=bool(flag)))
            if flag == 0:
                return dict(code=number, source_taxid=start, inheritance_path=trace,
                            recorded_values_agree=all(q['recorded_code'] == number for q in trace))
            start = value['parent']
    contexts = {}
    with Path(identity_table).open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            taxon = row['taxon_id']
            assert taxon not in contexts
            species = resolve_id(row['effective_species_taxid'])
            assert species == int(row['canonical_effective_taxid'])
            taxon_path = path(species)
            assert taxon_path == list(map(int, row['node_path_taxids'].split(';')))
            assembly = None
            if row['catalogue_taxid']:
                catalog = resolve_id(row['catalogue_taxid'])
                assembly = dict(original_catalogue_taxid=row['catalogue_taxid'], canonical_taxid=catalog,
                                nuclear=code(catalog, 'nuclear'), mitochondrial=code(catalog, 'mitochondrial'))
            contexts[taxon] = dict(taxon_id=taxon, study_role=row['study_role'],
                frozen_species_name=row['frozen_species_name'], snapshot_species_name=row['snapshot_scientific_name'],
                canonical_species_taxid=species, source_parent_path=taxon_path,
                nuclear=code(species, 'nuclear'), mitochondrial=code(species, 'mitochondrial'),
                assembly_catalogue_code_context=assembly, identity_review_flags=row['identity_review_flags'],
                biological_species_delimitation=row['biological_species_delimitation'],
                compartment_assignment='unverified', genetic_code_admission=False)
            assert str(contexts[taxon]['nuclear']['code']) in codes
            assert str(contexts[taxon]['mitochondrial']['code']) in codes
    return contexts, codes


def relation(number, ctx):
    if number is None or number == '':
        return 'inherited_test_code_not_recorded'
    number = int(number)
    if number == 0:
        return 'inherited_test_code_unspecified'
    same_nuclear = number == ctx['nuclear']['code']
    same_mitochondrial = number == ctx['mitochondrial']['code']
    if same_nuclear and same_mitochondrial:
        return 'matches_snapshot_nuclear_and_mitochondrial_codes_compartment_unverified'
    if same_nuclear:
        return 'matches_snapshot_nuclear_code_compartment_unverified'
    if same_mitochondrial:
        return 'matches_snapshot_mitochondrial_code_requires_compartment_verification'
    return 'differs_from_snapshot_nuclear_and_mitochondrial_codes_requires_review'


def original_target(source, ctx):
    old = source['translation_evidence']
    return dict(ordinal=source['ordinal'], cds_id=source['cds_id'], protein_id=source['protein_id'],
        original_target_dna_sha256=source['original_target_dna_sha256'], original_joined_status=source['joined_status'],
        inherited_translation_status=old['inherited_status'], inherited_test_code=old['translation_table'],
        inherited_code_provenance=old['code_provenance'], snapshot_code_relationship=relation(old['translation_table'], ctx),
        taxon_context_id=ctx['taxon_id'], genetic_code_admission=False)


def sources(path):
    with gzip.open(path, 'rt') as f:
        for line in f:
            yield loads(line)


def check_taxon(entry, report, ctx):
    expected_targets, expected_products, expected_cross = Counter(), Counter(), Counter()
    targets = iter(sources(report['artifact_paths'][0]))
    for raw in sources(entry['targets']):
        expected = original_target(raw, ctx)
        assert next(targets) == expected
        expected_targets[expected['snapshot_code_relationship']] += 1
    assert next(targets, None) is None
    outputs = iter(sources(report['artifact_paths'][1]))
    selected, supplements = 0, 0
    for raw in sources(entry['products']):
        original = [original_target(t, ctx) for t in raw['original_target_evidence']]
        derived = None
        if raw['separate_derived_evidence'] is not None:
            source = raw['separate_derived_evidence']
            assert source['kind'] in ['separate_annotation_boundary_codon_projection', 'separate_strict_annotation_derived_cds']
            derived = dict(kind=source['kind'], test_code=1, source='table_1_test_in_original_separate_derived_audit',
                exact_derived_translation=source['exact_derived_translation'], snapshot_code_relationship=relation(1, ctx),
                shared_original_genome_dependency=source['shared_original_genome_dependency'])
        labels = {t['snapshot_code_relationship'] for t in original}
        if len(labels) == 1:
            label = list(labels)[0]
        elif len(labels) > 1:
            label = 'multiple_original_code_context_relationships_require_review'
        elif derived is not None:
            label = 'separate_derived_test:' + derived['snapshot_code_relationship']
        else:
            label = 'no_original_or_separate_derived_test_code'
        availability = raw['representative_availability']
        availability = availability['availability'] if availability is not None else 'alternative_not_assigned'
        expected = dict(protein_id=raw['protein_id'], sequence_sha256=raw['sequence_sha256'],
            selected_representative=raw['selected_representative'], original_joined_status=raw['joined_status'],
            original_target_code_contexts=original, separate_derived_code_context=derived, availability=availability,
            snapshot_code_classification=label, taxon_context_id=ctx['taxon_id'], original_target_count=len(original),
            genetic_code_admission=False, biological_codon_eligibility=False)
        assert next(outputs) == expected
        expected_products[label] += 1
        expected_cross[raw['joined_status'] + ':' + availability + ':' + label] += 1
        selected += int(raw['selected_representative'])
        supplements += int(derived is not None)
    assert next(outputs, None) is None
    proof = dict(taxon_id=entry['taxon_id'], target_records=sum(expected_targets.values()),
                 source_products=sum(expected_products.values()), selected_representatives=selected,
                 separate_derived_evidence_products=supplements, target_relationship_counts=dict(expected_targets),
                 product_classification_counts=dict(expected_products), coding_model_code_cross_counts=dict(expected_cross))
    for key, value in proof.items():
        assert report[key] == value
    return proof


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'producer-receipt', 'producer-transport', 'output', 'receipt'):
        p.add_argument('--' + name, type=Path, required=True)
    args = p.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    plan, producer, transport = [loads(q.read_text()) for q in [args.plan, args.producer_receipt, args.producer_transport]]
    assert producer['status'] == 'complete_full_genetic_code_context_pending_independent_readback'
    assert transport['original_tool_terminal_exit_code'] == 0
    verify(transport['source_hashes'])
    verify(plan['pins'])
    start = time.monotonic()
    contexts, codes = source_contexts(plan['taxonomy_archive'], plan['identity_table'])
    context_receipt = loads(Path(plan['context_receipt']).read_text())
    assert contexts == context_receipt['contexts'] and codes == context_receipt['codebook']
    assert len(contexts) == 526
    assert Counter(c['study_role'] for c in contexts.values()) == dict(ingroup=501, outgroup=25)
    reports = {r['taxon_id']: r for r in producer['taxa_reports']}
    assert len(reports) == len(producer['taxa_reports']) == 526
    assert set(contexts) == set(reports) == {e['taxon_id'] for e in plan['entries']}
    args.output.mkdir(exist_ok=False)
    context_path = args.output / 'independent_taxonomic_code_context.json'
    context_path.write_text(json.dumps(contexts, indent=2) + '\n')
    proofs = []
    with ProcessPoolExecutor(max_workers=plan['cpu']) as pool:
        futures = []
        for entry in plan['entries']:
            assert entry['context'] == contexts[entry['taxon_id']]
            futures.append(pool.submit(check_taxon, entry, reports[entry['taxon_id']], contexts[entry['taxon_id']]))
        for future in as_completed(futures):
            proofs.append(future.result())
            print('independent_genetic_code_context_replay', len(proofs), '/526', flush=True)
    assert sum(r['source_products'] for r in proofs) == 5927745
    assert sum(r['selected_representatives'] for r in proofs) == 5815847
    assert sum(r['target_records'] for r in proofs) == 5923039
    pins = dict(producer['source_hashes'])
    for path, digest in transport['source_hashes'].items():
        bind(pins, path, digest)
    for path in [args.plan, args.producer_receipt, args.producer_transport, context_path, Path(__file__)]:
        bind(pins, path)
    verify(pins)
    result = dict(status='passed_all526_taxonomy_inheritance_and_full_coding_code_context_replay',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, source_products=5927745,
        selected_representatives=5815847, target_records=5923039, taxa_checks=sorted(proofs, key=lambda r: r['taxon_id']),
        elapsed_seconds=time.monotonic() - start, source_hashes=pins, scientific_eligibility=False,
        genetic_code_admission=False, biological_codon_eligibility=False, translations_recomputed=False,
        gpu=False, new_predictions=0,
        scope='Independent raw taxonomy/codebook parser and complete ancestor traversal plus all original '
              'target/product code-context records replayed without producer/helper imports. Original '
              'coding/identity evidence and Python libraries remain shared; no independent translation, '
              'compartment determination, biological taxonomy/genetic-code adjudication, codon selection '
              'or evolutionary admission.')
    with args.receipt.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'taxa_checks']}, indent=2))


if __name__ == '__main__':
    main()
