#!/usr/bin/env python3
"""Literal coding/model-link controls; no biological source subset is run."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path

from coding_structure_source_coupling_v1 import target_record, product_record, EXACT, REVIEW, DERIVED, MISSING
from reference_measurement_union_sources import bind


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--receipt', type=Path, required=True)
    args = p.parse_args()
    assert not args.receipt.exists()
    target = dict(ordinal=1, cds_id='dna1', product_id='p1', description='dna1 original header',
        sequence_sha256='d' * 64, sequence_length=9, status='one_exact_unmodified_genomic_cds_candidate',
        candidate_count=1, matching_candidate_indices=[0])
    inherited = dict(inherited_status='exact_translation', protein_id='p1',
                     original_header=target['description'], translation_table='12',
                     code_provenance='explicit_gff_transl_table')
    joined = target_record(target, inherited)
    assert joined['joined_status'] == EXACT and not joined['biological_codon_eligibility']
    assert joined['translation_evidence']['translation_table'] == '12'
    controls = ['one_source_bound_exact_pair_preserves_code_and_no_admission']
    assumed = target_record(target, dict(inherited, translation_table='1',
        code_provenance='table_1_assumption_no_explicit_code'))
    assert assumed['translation_evidence']['code_provenance'] == 'table_1_assumption_no_explicit_code'
    assert not assumed['biological_codon_eligibility']
    controls.append('default_code_assumption_retained_without_admission')
    for status in ['no_exact_unmodified_genomic_cds_match',
                   'exact_genomic_candidate_with_annotation_exception_requires_review',
                   'multiple_exact_genomic_candidates_require_review']:
        changed = dict(target, status=status)
        assert target_record(changed, inherited)['joined_status'] == REVIEW
        controls.append('retained_' + status)
    assert target_record(target, dict(inherited, inherited_status='translation_mismatch'))['joined_status'] == REVIEW
    controls.append('exact_genomic_dna_does_not_override_translation_mismatch')
    for name, changed in [('product_identity', dict(inherited, protein_id='foreign')),
                          ('header', dict(inherited, original_header='foreign')),
                          ('dna_digest', dict(inherited, original_dna_sha256='x' * 64))]:
        try:
            target_record(target, changed)
        except AssertionError:
            controls.append('reject_' + name)
        else:
            raise AssertionError(name)
    product = dict(protein_id='p1', protein_length=3, sequence_sha256='a' * 64,
        selected_representative=True, original_mapping={'status': 'unresolved_gene'},
        original_decision={'decision': 'baseline_only'}, genomic_candidate_count=1,
        original_cds_targets=[dict(cds_id='dna1', ordinal=1, status=target['status'])], missing_target_reason=None)
    availability = dict(protein_id='p1', sequence_sha256='a' * 64, length=3, availability='both')
    row = product_record(product, [joined], None, availability)
    assert row['joined_status'] == EXACT and row['original_mapping'] == product['original_mapping']
    assert not row['biological_codon_eligibility'] and row['representative_availability']['availability'] == 'both'
    controls.append('both_predictors_and_unresolved_gene_preserved_without_copy_admission')
    derived = dict(kind='separate_projection', exact_derived_translation=True, shared_original_genome_dependency=True)
    empty = dict(product, original_cds_targets=[], missing_target_reason='no_publisher_cds')
    row = product_record(empty, [], derived, availability)
    assert row['joined_status'] == DERIVED and not row['original_target_evidence']
    assert row['separate_derived_evidence']['shared_original_genome_dependency']
    controls.append('derived_span_does_not_replace_missing_original_target')
    assert product_record(empty, [], None, availability)['joined_status'] == MISSING
    controls.append('missing_original_and_derived_sources_retained')
    alt = dict(empty, selected_representative=False)
    row = product_record(alt, [], None, None)
    assert row['availability_scope'] == 'alternative_product_not_assigned_by_this_baseline'
    assert row['representative_availability'] is None
    controls.append('alternative_not_assigned_is_not_a_claim_of_no_model')
    for name, changed in [('model_sequence', dict(availability, sequence_sha256='z' * 64)),
                          ('model_length', dict(availability, length=4)),
                          ('model_identity', dict(availability, protein_id='foreign'))]:
        try:
            product_record(product, [joined], None, changed)
        except AssertionError:
            controls.append('reject_' + name)
        else:
            raise AssertionError(name)
    duplicated = copy.deepcopy(product)
    second = dict(joined, ordinal=2, cds_id='dna2')
    duplicated['original_cds_targets'].append(dict(cds_id='dna2', ordinal=2, status=target['status']))
    assert product_record(duplicated, [joined, second], None, availability)['joined_status'] == REVIEW
    controls.append('multiple_targets_preserved_as_product_review')
    assert len(controls) == 17
    pins = {}
    for path in [Path(__file__), Path('scripts/coding_structure_source_coupling_v1.py')]:
        bind(pins, path)
    result = dict(status='passed_literal_coding_structure_source_coupling_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(), controls=controls, control_count=len(controls),
        source_hashes=pins, scientific_eligibility=False, gpu=False, new_predictions=0,
        scope='Literal source identity/header/hash/model binding, multiple/missing/derived evidence and no-admission controls. '
              'No real taxon subset, sequence translation, model prediction or biological inference.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
