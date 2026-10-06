#!/usr/bin/env python3
"""Summarize closed full coding/model source joins without biological admission."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


EXACT = 'original_target_genome_agrees_and_inherited_strict_translation_exact'
REVIEW = 'original_target_review_or_no_strict_translation_agreement'
DERIVED = 'no_original_target_with_separate_derived_translation_evidence'
MISSING = 'no_original_target_or_exact_derived_translation_evidence'
PRODUCT_STATUSES = [EXACT, REVIEW, DERIVED, MISSING]
AVAILABILITY = ['afdb_only', 'esmfold_only', 'both', 'neither', 'alternative_not_assigned']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--taxon-table', type=Path, required=True)
    args = p.parse_args()
    assert not args.output.exists() and not args.taxon_table.exists()
    paths = [Path('metadata', name + '_20261005_v1.json') for name in
             ['full_coding_structure_source_coupling', 'full_coding_structure_source_coupling_transport',
              'full_coding_structure_source_coupling_readback',
              'full_coding_structure_source_coupling_readback_transport']]
    producer, pt, reader, rt = [json.loads(q.read_text()) for q in paths]
    assert producer['status'] == 'complete_full_coding_structure_source_coupling_pending_independent_readback'
    assert reader['status'] == 'passed_all526_coding_structure_source_join_and_independent_model_TSV_replay'
    plan_path = Path('metadata/full_coding_structure_source_coupling_plan_20261005_v1.json')
    plan = json.loads(plan_path.read_text())
    entries = {r['taxon_id']: r for r in plan['entries']}
    originals = {r['taxon_id']: r for r in producer['taxa_reports']}
    proofs = {r['taxon_id']: r for r in reader['taxa_checks']}
    assert len(entries) == len(plan['entries']) == len(originals) == len(producer['taxa_reports']) == 526
    assert len(proofs) == len(reader['taxa_checks']) == 526
    assert entries.keys() == originals.keys() == proofs.keys()
    assert Counter(e['study_role'] for e in entries.values()) == dict(ingroup=501, outgroup=25)
    fields = ['source_products', 'selected_representatives', 'target_records', 'target_status_counts',
              'product_status_counts', 'coding_structure_cross_counts', 'derived_evidence_products']
    totals = Counter()
    aggregates = {f: Counter() for f in fields if f.endswith('_counts')}
    for taxon, proof in proofs.items():
        original, entry = originals[taxon], entries[taxon]
        for field in fields:
            assert proof[field] == original[field], (taxon, field)
        for field in ['source_products', 'selected_representatives', 'target_records']:
            assert proof[field] == entry[field]
        for report_field, plan_field in [('mapping_mode', 'mapping_mode'), ('study_role', 'study_role'),
                                        ('manifest_species_name', 'species_name'), ('manifest_lineage', 'lineage')]:
            assert original[report_field] == entry[plan_field]
        assert sum(proof['target_status_counts'].values()) == proof['target_records']
        assert sum(proof['product_status_counts'].values()) == proof['source_products']
        assert set(proof['target_status_counts']) <= {EXACT, REVIEW}
        assert set(proof['product_status_counts']) <= set(PRODUCT_STATUSES)
        by_status, by_model = Counter(), Counter()
        for key, value in proof['coding_structure_cross_counts'].items():
            status, model = key.rsplit(':', 1)
            assert status in PRODUCT_STATUSES and model in AVAILABILITY
            by_status[status] += value
            by_model[model] += value
        assert dict(by_status) == proof['product_status_counts']
        assert sum(v for k, v in by_model.items() if k != 'alternative_not_assigned') == proof['selected_representatives']
        assert by_model['alternative_not_assigned'] == proof['source_products'] - proof['selected_representatives']
        for field in fields:
            if field in aggregates:
                aggregates[field].update(proof[field])
            else:
                totals[field] += proof[field]
    assert {k: totals[k] for k in fields[:3]} == dict(source_products=5927745,
        selected_representatives=5815847, target_records=5923039)
    for field in fields[:3]:
        assert producer[field] == reader[field] == totals[field]
    for field in aggregates:
        assert dict(aggregates[field]) == producer[field]
    model_totals = Counter()
    for key, value in aggregates['coding_structure_cross_counts'].items():
        model_totals[key.rsplit(':', 1)[1]] += value
    assert dict(model_totals) == producer['representative_availability_and_alternative_counts']
    assert dict(model_totals) == dict(afdb_only=2994146, esmfold_only=24801, both=722,
                                     neither=2796178, alternative_not_assigned=111898)
    pins, waits = {}, []
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
        waits.append(dict(unit=transport['unit'], invocation_id=transport['invocation_id'],
                          original_tool_session_id=transport['original_tool_session_id'],
                          actual_terminal_exit_code=0))
    cross_fields = [s + ':' + a for s in PRODUCT_STATUSES for a in AVAILABILITY]
    columns = ['taxon_id', 'study_role', 'manifest_species_name', 'manifest_lineage', 'mapping_mode',
               'source_products', 'selected_representatives', 'target_records', 'derived_evidence_products',
               'original_exact_with_any_model', *cross_fields]
    with args.taxon_table.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for taxon in sorted(proofs):
            original, proof = originals[taxon], proofs[taxon]
            cross = proof['coding_structure_cross_counts']
            exact_model = sum(cross.get(EXACT + ':' + a, 0) for a in AVAILABILITY[:3])
            writer.writerow(dict(taxon_id=taxon,
                **{f: original[f] for f in columns[1:5]},
                **{f: proof[f] for f in columns[5:9]}, original_exact_with_any_model=exact_model,
                **{f: cross.get(f, 0) for f in cross_fields}))
    for path in [*paths, plan_path, args.taxon_table, Path(__file__)]:
        name, digest = str(path), sha(path)
        assert name not in pins or pins[name] == digest, name
        pins[name] = digest
    exact_models = sum(aggregates['coding_structure_cross_counts'][EXACT + ':' + a] for a in AVAILABILITY[:3])
    result = dict(status='complete_verified_full_coding_structure_source_coupling',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, fungal_entries=501, outgroup_entries=25,
        **dict(totals), **{f: dict(sorted(v.items())) for f, v in aggregates.items()},
        representative_availability_and_alternative_counts=dict(model_totals),
        original_exact_with_any_model=exact_models, taxon_table=str(args.taxon_table),
        original_transports=waits, complete_bound_files=len(pins), source_hashes=pins,
        scientific_eligibility=False, biological_codon_eligibility=False,
        inherited_translation_independently_recomputed=False, gpu=False, new_predictions=0,
        all_eight_aims_incomplete=True,
        scope='Summary of separately closed full526entry original producer/reader source joins and every '
              'taxon disposition. Full reader uses model TSV, producer uses SQLite; inherited translation '
              'evidence is shared and not independently retranslated. No third original corpus scan by '
              'this merger. Original and derived CDS evidence remain distinct; gene/isoform ambiguity and '
              'all111898alternative products retained. Counts are source agreement before genetic-code '
              'review, alignment/divergence, homology/copy/taxonomy/contamination checks or confidence '
              'calibration. No codon-selection, physical structural-rate or evolutionary admission.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
