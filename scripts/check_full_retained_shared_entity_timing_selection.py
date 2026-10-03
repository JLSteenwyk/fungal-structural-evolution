#!/usr/bin/env python3
"""Exercise eligible full-grid timing selection separately from source qualification."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

from full_expanded_model_design_sources import OUTCOMES, digest
from full_retained_shared_entity_fit_sources import cases
from full_retained_shared_entity_timing import groups, READY
from run_ortholog_pair_guide_comparison import sha


def fixture():
    rng = np.random.default_rng(20261002)
    source = dict(fit_contract='synthetic-retained-selection-contract',retained_contract='synthetic-retained-contract')
    plan = dict(methods=['ml', 'reml'], loading_modes=['signed','unsigned'],trees=['tree-' + str(i) for i in range(5)])
    cohort = dict(cohort_id='synthetic-cohort', ordered_case_ids_sha256=digest(list(range(16))))
    entries = []
    for i in range(30):
        columns = 1 + i % 3
        matrix = np.column_stack([np.ones(16), rng.normal(size=(16, columns - 1))])
        design = dict(design_id='design-' + str(i), raw_design_sha256=digest(matrix.tolist()),
            active_column_indices=list(range(columns)), exactly_zero_column_indices=[],
            predictor_columns=['intercept', 'sequence', 'sequence_squared'][:columns])
        responses = {outcome: (dict(fit_input_id=digest([i, outcome]), outcome=outcome,
            records=16, response_sha256=digest([i, outcome, 'response']),
            disposition='ready_for_working_covariance_fit'), rng.normal(size=16))
            for outcome in OUTCOMES}
        # Synthetic declared qualifications isolate census/selection bookkeeping.
        # Numerical qualification is tested by the separate six-variance probes.
        audits = {(mode, tree): dict(audit_id=digest([i, mode, tree]), disposition=READY,
            numerical_audit=dict(normalized_design_condition_number=float(1 + i % 5)),
            retained_kernel_names=['residual','background_node','model_pair','family_intercept','species'])
            for mode in ['signed', 'unsigned'] for tree in plan['trees']}
        audits={key:(audit,dict(audit_id='original-'+audit['audit_id']),dict(explicit_synthetic_certificate=True)) for key,audit in audits.items()}
        entries.append((design, matrix, responses, audits))
    return source, plan, cohort, entries


def check(source, plan, cohort, entries):
    census, selected, counts = groups(source, plan, cohort, entries)
    expected = {}
    identities = {}
    for identity, matrix, response,_,_,_ in cases(source, plan, cohort, entries):
        identities[identity['candidate_id']] = identity
        key = tuple(identity[k] for k in ['loading_mode', 'tree', 'method', 'outcome'])
        if identity['source_combined_disposition'] != READY:
            continue
        entry = next(e for e in entries if e[0]['design_id'] == identity['design_id'])
        condition = entry[3][identity['loading_mode'], identity['tree']][0]['numerical_audit']['normalized_design_condition_number']
        expected.setdefault(key, []).append((matrix.shape[1], condition, identity['candidate_id']))
    assert len(census) == len(identities) == 1200
    assert len(selected) == len(expected)
    for group_id, representative in selected.items():
        values = expected[representative['key']]
        assert representative['rank'] == max(values)
        assert counts[group_id] == len(values)
    for row in census:
        original = identities[row['candidate_id']]
        assert row['identity_sha256'] == digest(original)
        assert row['source_disposition'] == original['source_combined_disposition']
        assert row['scientific_eligibility'] is False
        if original['source_combined_disposition'] == READY:
            assert row['benchmark_candidate_id'] == selected[row['group_id']]['identity']['candidate_id']
        else:
            assert row['group_id'] is None and row['benchmark_candidate_id'] is None
    return census, selected, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source, plan, cohort, entries = fixture()
    census, selected, counts = check(source, plan, cohort, entries)
    assert len(selected) == 40 and set(counts.values()) == {30}
    assert all(r['rank'][:2] == (3, 5.) for r in selected.values())
    reversed_result = check(source, plan, cohort, list(reversed(entries)))
    assert sorted(census, key=lambda r: r['candidate_id']) == sorted(reversed_result[0], key=lambda r: r['candidate_id'])
    assert {g: r['identity'] for g, r in selected.items()} == {g: r['identity'] for g, r in reversed_result[1].items()}
    mixed = copy.deepcopy(entries)
    for i, entry in enumerate(mixed):
        if i % 2 == 0:
            entry[2][OUTCOMES[0]][0]['disposition'] = 'constant_response_requires_review'
        if i % 3 == 0:
            entry[3]['signed', plan['trees'][0]][0]['disposition'] = 'synthetic_covariance_review'
    mixed_census, mixed_selected, mixed_counts = check(source, plan, cohort, mixed)
    eligible = sum(r['group_id'] is not None for r in mixed_census)
    assert 0 < eligible < 1200 and eligible == sum(mixed_counts.values())
    excluded = copy.deepcopy(entries)
    for entry in excluded:
        for audit in entry[3].values():
            audit[0]['disposition'] = 'synthetic_covariance_review'
    excluded_census, excluded_selected, excluded_counts = check(source, plan, cohort, excluded)
    assert not excluded_selected and not excluded_counts
    assert all(r['group_id'] is None for r in excluded_census)
    paths = [Path(__file__), Path('scripts/full_retained_shared_entity_timing.py'), Path('scripts/full_retained_shared_entity_fit_sources.py')]
    result = dict(status='passed_full_retained_shared_entity_timing_eligible_selection_contracts_v1',
        synthetic_candidate_census=1200, full_eligible_timing_groups=40,
        eligible_candidates_per_full_group=30, mixed_source_eligible_candidates=eligible,
        ranking_by_dimension_then_condition_then_candidate_checked=True,
        full_grid_order_invariance_checked=True, source_review_and_empty_groups_retained=True,
        source_hashes={str(p): sha(p) for p in paths}, scientific_eligibility=False,
        scope='Synthetic complete 30-design/two-response/two-method/two-loading/five-tree selection contracts. Declared source qualification is synthetic; this does not independently qualify covariance, measure production timing, accept fits, or change any pinned queued source.')
    with args.output.open('x') as f:
        f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
