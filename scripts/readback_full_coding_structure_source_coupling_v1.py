#!/usr/bin/env python3
"""Replay every coding/model join from source tables and an independent TSV model view."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import io
import json
from pathlib import Path
import time

from reference_measurement_union_sources import bind, verify


EXACT = 'original_target_genome_agrees_and_inherited_strict_translation_exact'
REVIEW = 'original_target_review_or_no_strict_translation_agreement'
DERIVED = 'no_original_target_with_separate_derived_translation_evidence'
MISSING = 'no_original_target_or_exact_derived_translation_evidence'


def load(text):
    def pairs(items):
        out = {}
        for key, value in items:
            assert key not in out
            out[key] = value
        return out
    return json.loads(text, object_pairs_hook=pairs,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def rows(path):
    with open(path) as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def records(path):
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            yield load(line)


def availability_index(path):
    ranges = {}
    count = 0
    with open(path, 'rb') as handle:
        header = handle.readline()
        previous = None
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line:
                if previous is not None:
                    ranges[previous]['end'] = offset
                break
            taxon = line.split(b'\t', 1)[0].decode()
            if taxon != previous:
                if previous is not None:
                    ranges[previous]['end'] = offset
                assert taxon not in ranges
                ranges[taxon] = dict(start=offset, records=0)
                previous = taxon
            ranges[taxon]['records'] += 1
            count += 1
    assert len(ranges) == 526 and count == 5815847
    return header.decode(), ranges


def availability_block(path, header, span, taxon):
    with open(path, 'rb') as handle:
        handle.seek(span['start'])
        block = handle.read(span['end'] - span['start']).decode()
    result = {}
    for row in csv.DictReader(io.StringIO(header + block), delimiter='\t'):
        assert row['taxon_id'] == taxon and row['protein_id'] not in result
        afdb, esm = bool(row['afdb_model_id']), bool(row['esmfold_model_id'])
        state = 'both' if afdb and esm else 'afdb_only' if afdb else 'esmfold_only' if esm else 'neither'
        assert state == row['availability']
        result[row['protein_id']] = dict(protein_id=row['protein_id'], sequence_sha256=row['sequence_sha256'],
                                        length=int(row['length']), availability=state)
    assert len(result) == span['records']
    return result


def translation_evidence(entry):
    primary, supplemental = {}, {}
    mode = entry['mapping_mode']
    if mode == 'ncbi_protein_gff':
        for number, row in enumerate(rows(entry['translation_table']), 1):
            assert row['taxon_id'] == entry['taxon_id']
            if row['status'] == 'protein_without_cds_record':
                continue
            assert row['cds_id'] not in primary
            primary[row['cds_id']] = dict(inherited_status=row['status'], protein_id=row['protein_id'] or None,
                translation_table=row['translation_table'] or None, code_provenance=row['code_source'] or None,
                terminal_stop=row['terminal_stop'] or None, annotation_flags=row['gff_flags'], original_header=row['cds_header'],
                source_kind='ncbi_unmodified_publisher_cds', evidence_path=entry['translation_table'],
                evidence_ordinal=number, translation_independently_recomputed_here=False)
    elif mode == 'transcript_gff':
        for number, row in enumerate(rows(entry['translation_table']), 1):
            if row['taxon_id'] == entry['taxon_id']:
                assert row['protein_id'] not in primary
                primary[row['protein_id']] = dict(inherited_status=row['status'], protein_id=row['protein_id'],
                    translation_table='1', code_provenance='table_1_test_in_original_publisher_audit',
                    terminal_stop=row['terminal_stop_in_cds'], annotation_flags=None,
                    source_kind='external_unmodified_publisher_cds', evidence_path=entry['translation_table'],
                    evidence_ordinal=number, original_dna_sha256=row['cds_sequence_sha256'],
                    original_protein_sha256=row['protein_sequence_sha256'], translation_independently_recomputed_here=False)
        for number, row in enumerate(rows(entry['supplementary_table']), 1):
            if row['taxon_id'] == entry['taxon_id']:
                assert row['protein_id'] not in supplemental
                supplemental[row['protein_id']] = dict(kind='separate_annotation_boundary_codon_projection',
                    evidence_path=entry['supplementary_table'], evidence_ordinal=number, original_row=row,
                    shared_original_genome_dependency=True,
                    exact_derived_translation=row['projection_status'] == 'exact_genome_linked_codon_translation')
    elif mode == 'verified_orf_coordinates':
        for number, row in enumerate(rows(entry['translation_table']), 1):
            assert row['protein_id'] not in primary
            primary[row['protein_id']] = dict(inherited_status=row['status'], protein_id=row['protein_id'],
                translation_table='1', code_provenance='table_1_test_in_original_dependent_orf_audit',
                terminal_stop=None, annotation_flags=None, source_kind='dependent_genome_derived_orf_cds',
                evidence_path=entry['translation_table'], evidence_ordinal=number,
                shared_original_genome_dependency=True, translation_independently_recomputed_here=False)
    else:
        assert mode == 'creolimax_gtf'
        for number, row in enumerate(rows(entry['supplementary_table']), 1):
            assert row['protein_id'] not in supplemental
            supplemental[row['protein_id']] = dict(kind='separate_strict_annotation_derived_cds',
                evidence_path=entry['supplementary_table'], evidence_ordinal=number, original_row=row,
                shared_original_genome_dependency=True, exact_derived_translation=row['status'] == 'exact_translation')
    return primary, supplemental


def check_taxon(entry, report, model_table, header, span):
    inherited, supplements = translation_evidence(entry)
    models = availability_block(model_table, header, span, entry['taxon_id'])
    linked = defaultdict(list)
    target_counts, product_counts, cross = Counter(), Counter(), Counter()
    targets = 0
    outputs = iter(records(report['artifact_paths'][0]))
    for raw in records(entry['genomic_targets']):
        proof = inherited[raw['cds_id']]
        assert raw['product_id'] == proof['protein_id']
        if 'original_header' in proof:
            assert raw['description'] == proof['original_header']
        if 'original_dna_sha256' in proof:
            assert raw['sequence_sha256'] == proof['original_dna_sha256']
        status = EXACT if raw['status'] == 'one_exact_unmodified_genomic_cds_candidate' and proof['inherited_status'] == 'exact_translation' else REVIEW
        expected = dict(ordinal=raw['ordinal'], cds_id=raw['cds_id'], protein_id=raw['product_id'],
            original_target_dna_sha256=raw['sequence_sha256'], original_target_length=raw['sequence_length'],
            genomic_status=raw['status'], genomic_candidate_count=raw['candidate_count'],
            matching_candidate_indices=raw['matching_candidate_indices'], translation_evidence=proof,
            joined_status=status, biological_codon_eligibility=False)
        assert next(outputs) == expected
        target_counts[status] += 1
        targets += 1
        if raw['product_id'] is not None:
            linked[raw['product_id']].append(expected)
    assert next(outputs, None) is None
    outputs = iter(records(report['artifact_paths'][1]))
    products = selected = 0
    used_models, used_supplements = set(), set()
    for raw in records(entry['genomic_products']):
        pid = raw['protein_id']
        evidence = linked.get(pid, [])
        assert raw['original_cds_targets'] == [dict(cds_id=e['cds_id'], ordinal=e['ordinal'], status=e['genomic_status']) for e in evidence]
        model = models.get(pid) if raw['selected_representative'] else None
        if raw['selected_representative']:
            assert model is not None and model['sequence_sha256'] == raw['sequence_sha256'] and model['length'] == raw['protein_length']
            used_models.add(pid)
        supplement = supplements.get(pid)
        if supplement is not None:
            used_supplements.add(pid)
        status = (EXACT if len(evidence) == 1 and evidence[0]['joined_status'] == EXACT else REVIEW if evidence else
                  DERIVED if supplement is not None and supplement['exact_derived_translation'] else MISSING)
        expected = dict(protein_id=pid, protein_length=raw['protein_length'], sequence_sha256=raw['sequence_sha256'],
            selected_representative=raw['selected_representative'], original_mapping=raw['original_mapping'],
            original_decision=raw['original_decision'], original_genomic_candidate_count=raw['genomic_candidate_count'],
            original_missing_target_reason=raw['missing_target_reason'], original_target_evidence=evidence,
            separate_derived_evidence=supplement, representative_availability=model,
            availability_scope='selected_representative_baseline' if raw['selected_representative'] else 'alternative_product_not_assigned_by_this_baseline',
            joined_status=status, biological_codon_eligibility=False, scientific_eligibility=False)
        for target in evidence:
            if 'original_protein_sha256' in target['translation_evidence']:
                assert target['translation_evidence']['original_protein_sha256'] == raw['sequence_sha256']
        assert next(outputs) == expected
        products += 1
        selected += int(raw['selected_representative'])
        product_counts[status] += 1
        cross[status + ':' + (model['availability'] if model else 'alternative_not_assigned')] += 1
    assert next(outputs, None) is None
    assert used_models == set(models) and used_supplements == set(supplements)
    proof = dict(taxon_id=entry['taxon_id'], source_products=products, selected_representatives=selected,
                 target_records=targets, target_status_counts=dict(target_counts), product_status_counts=dict(product_counts),
                 coding_structure_cross_counts=dict(cross), derived_evidence_products=len(used_supplements))
    for key, value in proof.items():
        assert report[key] == value
    return proof


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['plan', 'producer-receipt', 'producer-transport', 'output', 'receipt']:
        p.add_argument('--' + name, type=Path, required=True)
    args = p.parse_args()
    assert not args.receipt.exists()
    plan, producer, transport = [json.loads(q.read_text()) for q in [args.plan, args.producer_receipt, args.producer_transport]]
    assert producer['status'] == 'complete_full_coding_structure_source_coupling_pending_independent_readback'
    assert transport['original_tool_terminal_exit_code'] == 0
    verify(transport['source_hashes'])
    verify(plan['pins'])
    root = args.output
    root.mkdir(exist_ok=False)
    header, spans = availability_index(plan['availability_table'])
    assert set(spans) == {e['taxon_id'] for e in plan['entries']}
    index_path = root / 'independent_availability_byte_index.json'
    index_path.write_text(json.dumps(spans, indent=2) + '\n')
    by_taxon = {r['taxon_id']: r for r in producer['taxa_reports']}
    assert len(by_taxon) == 526
    start = time.monotonic()
    proofs = []
    with ProcessPoolExecutor(max_workers=plan['cpu']) as pool:
        futures = [pool.submit(check_taxon, e, by_taxon[e['taxon_id']], plan['availability_table'], header, spans[e['taxon_id']])
                   for e in plan['entries']]
        for future in as_completed(futures):
            proofs.append(future.result())
            print('independent_coding_structure_source_replay', len(proofs), '/526', flush=True)
    assert sum(r['source_products'] for r in proofs) == 5927745
    assert sum(r['selected_representatives'] for r in proofs) == 5815847
    assert sum(r['target_records'] for r in proofs) == 5923039
    pins = dict(producer['source_hashes'])
    for path, digest in transport['source_hashes'].items():
        bind(pins, path, digest)
    for path in [args.plan, args.producer_receipt, args.producer_transport, Path(__file__), index_path]:
        bind(pins, path)
    verify(pins)
    result = dict(status='passed_all526_coding_structure_source_join_and_independent_model_TSV_replay',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, source_products=5927745,
        selected_representatives=5815847, target_records=5923039, taxa_checks=sorted(proofs, key=lambda r: r['taxon_id']),
        elapsed_seconds=time.monotonic() - start, source_hashes=pins, scientific_eligibility=False,
        biological_codon_eligibility=False, inherited_translation_independently_recomputed=False,
        gpu=False, new_predictions=0,
        scope='All original source/protein/target classifications replayed without producer imports. Models reconstructed '
              'from the complete TSV rather than producer SQLite. Shared original evidence/JSON/CSV/gzip/hash dependencies '
              'remain explicit. Existing translation classifications are inherited, not independently retranslated here. '
              'No old-derived/original source substitution, canonical isoform claim, confidence/copy/taxonomy/homology, '
              'codon alignment/divergence, selection or evolutionary admission.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'taxa_checks']}, indent=2))


if __name__ == '__main__':
    main()
