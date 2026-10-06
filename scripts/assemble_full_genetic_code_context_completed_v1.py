#!/usr/bin/env python3
"""Summarize closed full genetic-code context and every taxon disposition."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--taxon-table', type=Path, required=True)
    args = p.parse_args()
    assert not args.output.exists() and not args.taxon_table.exists()
    paths = [Path('metadata', stem + '_20261006_v1.json') for stem in
             ['full_genetic_code_context', 'full_genetic_code_context_transport',
              'full_genetic_code_context_readback', 'full_genetic_code_context_readback_transport']]
    producer, pt, reader, rt = [json.loads(q.read_text()) for q in paths]
    assert producer['status'] == 'complete_full_genetic_code_context_pending_independent_readback'
    assert reader['status'] == 'passed_all526_taxonomy_inheritance_and_full_coding_code_context_replay'
    context_path = Path('metadata/full_taxonomic_genetic_code_context_20261006_v1.json')
    contexts = json.loads(context_path.read_text())['contexts']
    originals = {r['taxon_id']: r for r in producer['taxa_reports']}
    proofs = {r['taxon_id']: r for r in reader['taxa_checks']}
    assert len(originals) == len(producer['taxa_reports']) == len(proofs) == len(reader['taxa_checks']) == 526
    assert set(originals) == set(proofs) == set(contexts)
    fields = ['source_products', 'selected_representatives', 'target_records',
              'separate_derived_evidence_products', 'target_relationship_counts',
              'product_classification_counts', 'coding_model_code_cross_counts']
    totals = Counter()
    aggregates = {f: Counter() for f in fields if f.endswith('_counts')}
    for taxon, proof in proofs.items():
        for field in fields:
            assert originals[taxon][field] == proof[field], (taxon, field)
            if field in aggregates:
                aggregates[field].update(proof[field])
            else:
                totals[field] += proof[field]
        assert sum(proof['target_relationship_counts'].values()) == proof['target_records']
        assert sum(proof['product_classification_counts'].values()) == proof['source_products']
        assert sum(proof['coding_model_code_cross_counts'].values()) == proof['source_products']
    assert {f: totals[f] for f in fields[:3]} == dict(source_products=5927745,
        selected_representatives=5815847, target_records=5923039)
    assert totals['separate_derived_evidence_products'] == 67810
    for field in aggregates:
        assert dict(aggregates[field]) == producer[field]
    for field in fields[:3]:
        assert producer[field] == reader[field] == totals[field]
    pins, transports = {}, []
    for path, receipt, transport in [(paths[0], producer, pt), (paths[2], reader, rt)]:
        assert transport['status'] == 'verified_original_software_wait_and_whole_wrapper_payloads'
        assert transport['original_tool_terminal_exit_code'] == 0
        assert transport['validation_sha256'] == sha(path)
        assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
        assert transport['manager_start_records'] == transport['manager_completion_records'] == 1
        for mapping in [receipt['source_hashes'], transport['source_hashes']]:
            for name, digest in mapping.items():
                assert name not in pins or pins[name] == digest, name
                pins[name] = digest
        transports.append(dict(unit=transport['unit'], invocation_id=transport['invocation_id'],
                               original_tool_session_id=transport['original_tool_session_id'], actual_terminal_exit_code=0))
    exact = 'original_target_genome_agrees_and_inherited_strict_translation_exact'
    exact_modeled, model_counts = Counter(), Counter()
    per_taxon_exact = {}
    for taxon, proof in proofs.items():
        counts = Counter()
        for key, count in proof['coding_model_code_cross_counts'].items():
            status, availability, classification = key.split(':', 2)
            model_counts[availability] += count
            if status == exact and availability in ('afdb_only', 'esmfold_only', 'both'):
                counts[classification] += count
        exact_modeled.update(counts)
        per_taxon_exact[taxon] = counts
    assert dict(model_counts) == dict(afdb_only=2994146, esmfold_only=24801, both=722,
                                     neither=2796178, alternative_not_assigned=111898)
    assert sum(exact_modeled.values()) == 2952711
    assert Counter(c['study_role'] for c in contexts.values()) == dict(ingroup=501, outgroup=25)
    classifications = sorted(aggregates['product_classification_counts'])
    exact_classes = sorted(exact_modeled)
    columns = ['taxon_id', 'study_role', 'frozen_species_name', 'snapshot_species_name',
               'canonical_species_taxid', 'nuclear_code', 'nuclear_code_source_taxid',
               'mitochondrial_code', 'mitochondrial_code_source_taxid', 'assembly_nuclear_code',
               'identity_review_flags', *fields[:4], *['products:' + c for c in classifications],
               *['original_exact_modeled:' + c for c in exact_classes]]
    with args.taxon_table.open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for taxon in sorted(proofs):
            c, proof = contexts[taxon], proofs[taxon]
            assembly = c['assembly_catalogue_code_context']
            w.writerow(dict(taxon_id=taxon, **{k: c[k] for k in columns[1:5]},
                nuclear_code=c['nuclear']['code'], nuclear_code_source_taxid=c['nuclear']['source_taxid'],
                mitochondrial_code=c['mitochondrial']['code'], mitochondrial_code_source_taxid=c['mitochondrial']['source_taxid'],
                assembly_nuclear_code=assembly['nuclear']['code'] if assembly else '',
                identity_review_flags=c['identity_review_flags'], **{k: proof[k] for k in fields[:4]},
                **{'products:' + label: proof['product_classification_counts'].get(label, 0) for label in classifications},
                **{'original_exact_modeled:' + label: per_taxon_exact[taxon].get(label, 0) for label in exact_classes}))
    for path in [*paths, context_path, args.taxon_table, Path(__file__)]:
        name, digest = str(path), sha(path)
        assert name not in pins or pins[name] == digest, name
        pins[name] = digest
    result = dict(status='complete_verified_all526_coding_and_taxonomic_genetic_code_context',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, **dict(totals),
        selected_nuclear_code_counts=dict(Counter(c['nuclear']['code'] for c in contexts.values())),
        selected_mitochondrial_code_counts=dict(Counter(c['mitochondrial']['code'] for c in contexts.values())),
        **{f: dict(sorted(v.items())) for f, v in aggregates.items()},
        original_exact_modeled_code_classification_counts=dict(exact_modeled),
        representative_availability_and_alternatives=dict(model_counts), taxon_table=str(args.taxon_table),
        original_transports=transports, complete_bound_files=len(pins), source_hashes=pins,
        scientific_eligibility=False, genetic_code_admission=False, biological_codon_eligibility=False,
        translations_recomputed=False, gpu=False, new_predictions=0, all_eight_aims_incomplete=True,
        scope='Merger of closed full526source producer/independent reader; all taxon code/status/model '
              'dispositions retained without a third corpus scan. Original exact agreement and every '
              'review/alternative/derived source remain explicit. Code matches are pinned taxonomic '
              'assignments, not biological code/compartment verification; no original translation '
              'replacement, codon selection or evolutionary admission.')
    with args.output.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
