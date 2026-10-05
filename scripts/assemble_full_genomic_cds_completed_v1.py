#!/usr/bin/env python3
"""Join closed full genomic-CDS evidence and publish all taxon dispositions."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--taxon-table', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.taxon_table.exists()
    paths = [Path('metadata/' + name + '_20261005_v1.json') for name in
             ['full_genome_annotation_cds', 'full_genome_annotation_cds_transport',
              'full_genomic_cds_readback', 'full_genomic_cds_readback_transport']]
    producer, pt, reader, rt = [json.loads(p.read_text()) for p in paths]
    assert producer['status'] == 'complete_full_genome_annotation_cds_comparison_pending_independent_readback'
    assert reader['status'] == 'passed_full_independent_genome_annotation_cds_and_source_product_replay'
    assert producer['taxa'] == reader['taxa'] == 526
    originals = {r['taxon_id']: r for r in producer['taxa_reports']}
    proofs = {r['taxon_id']: r for r in reader['taxa_checks']}
    assert len(originals) == len(producer['taxa_reports']) == len(proofs) == len(reader['taxa_checks']) == 526
    assert originals.keys() == proofs.keys()
    fields = ['source_products', 'selected_representatives', 'target_records', 'feature_counts',
              'coordinate_counts', 'candidate_counts', 'target_status_counts']
    aggregates = {field: Counter() for field in fields if field.endswith('_counts')}
    for taxon, proof in proofs.items():
        for field in fields:
            assert proof[field] == originals[taxon][field], (taxon, field)
        for field in aggregates:
            aggregates[field].update(proof[field])
        assert sum(proof['target_status_counts'].values()) == proof['target_records']
        assert sum(proof['feature_counts'].values()) == sum(proof['coordinate_counts'].values())
    totals = {field: sum(r[field] for r in proofs.values()) for field in fields[:3]}
    assert totals == dict(source_products=5927745, selected_representatives=5815847, target_records=5923039)
    assert sum(aggregates['feature_counts'].values()) == 60917860 + 16588
    assert aggregates['feature_counts']['provisional_orf_row'] == 16588
    assert reader['target_records'] == totals['target_records']
    pins = {}
    waits = []
    for path, receipt, transport in [(paths[0], producer, pt), (paths[2], reader, rt)]:
        assert transport['status'] == 'verified_original_software_wait_and_whole_wrapper_payloads'
        assert transport['original_tool_terminal_exit_code'] == 0
        assert transport['validation_sha256'] == digest(path)
        assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
        assert transport['manager_start_records'] == transport['manager_completion_records'] == 1
        for mapping in [receipt['source_hashes'], transport['source_hashes']]:
            for name, sha in mapping.items():
                assert name not in pins or pins[name] == sha, name
                pins[name] = sha
        waits.append(dict(unit=transport['unit'], invocation_id=transport['invocation_id'],
                          original_tool_session_id=transport['original_tool_session_id'],
                          actual_terminal_exit_code=transport['original_tool_terminal_exit_code']))
    modes = Counter(r['mapping_mode'] for r in originals.values())
    statuses = sorted(aggregates['target_status_counts'])
    with args.taxon_table.open('x', newline='') as handle:
        columns = ['taxon_id', 'mapping_mode', *fields[:3], 'target_source_kind',
                   'missing_target_reason', *statuses]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        for taxon in sorted(proofs):
            original, proof = originals[taxon], proofs[taxon]
            target = original['target_cds_source']
            kind = 'none' if target is None else (
                'genome_derived_dependent_orf' if original['mapping_mode'] == 'verified_orf_coordinates' else
                'publisher_ncbi' if target['identifier_mode'] == 'ncbi_protein_id_header' else 'publisher_external')
            writer.writerow(dict(taxon_id=taxon, mapping_mode=original['mapping_mode'],
                **{f: proof[f] for f in fields[:3]}, target_source_kind=kind,
                missing_target_reason=original['missing_target_reason'] or '',
                **{status: proof['target_status_counts'].get(status, 0) for status in statuses}))
    for path in [*paths, args.taxon_table, Path(__file__)]:
        name, sha = str(path), digest(path)
        assert name not in pins or pins[name] == sha, name
        pins[name] = sha
    result = dict(status='complete_verified_full_genome_annotation_cds_comparison',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, **totals,
        feature_coordinate_rows=60917860 + 16588, original_annotation_features=60917860,
        provisional_orf_rows=16588, mapping_mode_counts=dict(modes),
        **{field: dict(sorted(counts.items())) for field, counts in aggregates.items()},
        taxon_table=str(args.taxon_table), original_transports=waits,
        complete_bound_files=len(pins), source_hashes=pins, scientific_eligibility=False,
        gpu=False, new_predictions=0, all_eight_aims_incomplete=True,
        scope='Merger of independently closed full producer/reader evidence with per-taxon and aggregate '
              'disposition agreement. No third genome/CDS corpus scan. Exact DNA agreement does not qualify '
              'translation, expression, canonical isoforms, gene copies, contamination, taxonomy, haplotigs, '
              'selection or evolutionary results. Two ORF-derived targets share their genome dependency; '
              'missing targets and every review record remain in scope.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
