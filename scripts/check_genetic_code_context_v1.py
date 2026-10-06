#!/usr/bin/env python3
"""Literal inheritance and code-context controls, without biological inference."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from genetic_code_context_v1 import ancestry, canonical, inherited_code, relationship, target_context, product_context
from reference_measurement_union_sources import bind


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--receipt', type=Path, required=True)
    args = p.parse_args()
    assert not args.receipt.exists()
    nodes = {1: (1, 1, 0, 0, 0), 2: (1, 12, 0, 4, 0), 3: (2, 12, 1, 4, 1),
             4: (3, 12, 1, 4, 1), 5: (2, 1, 1, 4, 1)}
    checks = []
    assert ancestry(4, nodes) == [4, 3, 2, 1]
    checks.append('complete_source_parent_path')
    assert canonical(20, {20: 21, 21: 4}) == 4
    checks.append('merged_taxid_chain')
    n, m = [inherited_code(4, nodes, c) for c in ['nuclear', 'mitochondrial']]
    assert n['code'] == 12 and m['code'] == 4 and n['source_taxid'] == m['source_taxid'] == 2
    assert n['recorded_values_agree'] and m['recorded_values_agree']
    checks.append('separate_compartment_inheritance')
    assert not inherited_code(5, nodes, 'nuclear')['recorded_values_agree']
    checks.append('inherited_record_value_disagreement_retained')
    context = dict(nuclear=n, mitochondrial=m)
    expected = {None: 'inherited_test_code_not_recorded', 0: 'inherited_test_code_unspecified',
        12: 'matches_snapshot_nuclear_code_compartment_unverified',
        4: 'matches_snapshot_mitochondrial_code_requires_compartment_verification',
        1: 'differs_from_snapshot_nuclear_and_mitochondrial_codes_requires_review'}
    for code, label in expected.items():
        assert relationship(code, context) == label
        checks.append(label)
    assert relationship(1, dict(nuclear=dict(code=1), mitochondrial=dict(code=1))) == 'matches_snapshot_nuclear_and_mitochondrial_codes_compartment_unverified'
    checks.append('same_code_does_not_resolve_compartment')
    context['taxon_id'] = 'literal'
    target = dict(ordinal=1, cds_id='cds', protein_id='protein', original_target_dna_sha256='source-dna',
                  joined_status='original_exact', translation_evidence=dict(inherited_status='exact_translation',
                  translation_table='1', code_provenance='table_1_assumption_no_explicit_code'))
    original = target_context(target, context)
    assert original['original_joined_status'] == 'original_exact'
    assert original['inherited_code_provenance'] == 'table_1_assumption_no_explicit_code'
    assert original['snapshot_code_relationship'] == expected[1] and not original['genetic_code_admission']
    checks.append('original_exact_and_code_assumption_retained_under_taxonomic_disagreement')
    product = dict(protein_id='protein', sequence_sha256='original-protein', selected_representative=True,
                   joined_status='original_exact', original_target_evidence=[target],
                   separate_derived_evidence=None, representative_availability=dict(availability='both'))
    joined = product_context(product, context)
    assert joined['availability'] == 'both' and joined['snapshot_code_classification'] == expected[1]
    assert not joined['genetic_code_admission'] and not joined['biological_codon_eligibility']
    checks.append('both_predictors_preserved_without_codon_admission')
    second = dict(target, ordinal=2, cds_id='other', translation_evidence=dict(target['translation_evidence'], translation_table='12'))
    multiple = product_context(dict(product, original_target_evidence=[target, second]), context)
    assert multiple['snapshot_code_classification'] == 'multiple_original_code_context_relationships_require_review'
    assert multiple['original_target_count'] == 2
    checks.append('multiple_original_target_code_relations_retained')
    supplement = dict(kind='separate_annotation_boundary_codon_projection', exact_derived_translation=True,
                      shared_original_genome_dependency=True)
    derived = product_context(dict(product, original_target_evidence=[], separate_derived_evidence=supplement), context)
    assert derived['original_target_count'] == 0
    assert derived['snapshot_code_classification'] == 'separate_derived_test:' + expected[1]
    assert derived['separate_derived_code_context']['shared_original_genome_dependency']
    checks.append('separate_derived_test_does_not_become_original_target')
    alternate = product_context(dict(product, selected_representative=False, representative_availability=None), context)
    assert alternate['availability'] == 'alternative_not_assigned' and not alternate['selected_representative']
    checks.append('alternative_baseline_scope_retained')
    missing = product_context(dict(product, original_target_evidence=[]), context)
    assert missing['snapshot_code_classification'] == 'no_original_or_separate_derived_test_code'
    checks.append('missing_code_source_retained')
    for name, action in [('merged_cycle', lambda: canonical(1, {1: 2, 2: 1})),
                         ('parent_cycle', lambda: ancestry(1, {1: (2, 1, 0, 0, 0), 2: (1, 1, 0, 0, 0)})),
                         ('inheritance_cycle', lambda: inherited_code(1, {1: (1, 1, 1, 0, 0)}, 'nuclear')),
                         ('missing_parent', lambda: ancestry(7, nodes))]:
        try:
            action()
        except AssertionError:
            checks.append('reject_' + name)
        else:
            raise AssertionError(name)
    pins = {}
    for path in [Path(__file__), Path('scripts/genetic_code_context_v1.py')]:
        bind(pins, path)
    result = dict(status='passed_literal_genetic_code_context_inheritance_controls',
        checked_utc=datetime.now(timezone.utc).isoformat(), controls=checks, control_count=len(checks),
        source_hashes=pins, scientific_eligibility=False, genetic_code_admission=False,
        gpu=False, new_predictions=0,
        scope='Literal taxonomy/code inheritance and distinct nuclear/mitochondrial source relationship tests; '
              'no biological species/code/compartment admission, translation or taxon pilot.')
    with args.receipt.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
