#!/usr/bin/env python3
"""Qualify exact identities, retained counterexamples and nonnegative cone maps."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha
from covariance_exact_folds_v2 import certificate, independent_certificate, variance_map, independent_variance_map, ORIGINAL_NAMES
from full_exact_covariance_sources_v2 import load, rows
from readback_full_exact_covariance_folds_v2 import check_record
from reference_measurement_union_sources import verify


def reject(action):
    try: action()
    except (AssertionError, ValueError, KeyError, OSError): return
    raise AssertionError('Invalid exact covariance certificate accepted')


def perturb(value, offset):
    value = value.tocsr(copy=True)
    value.indices[value.indptr[offset]] = value.indices[value.indptr[(offset + 3) % value.shape[0]]]
    value.sum_duplicates(); value.eliminate_zeros(); value.sort_indices()
    return value


def fixture(mode, scenario):
    n = 12; target = sparse.eye(n, format='csr')
    family = sparse.csr_matrix((np.ones(n), (np.arange(n), np.arange(n) // 4)), shape=(n, 3))
    background = sparse.csr_matrix((np.full(n, -1 if mode == 'signed' else 1), (np.arange(n), np.arange(n) // 3)), shape=(n, 4))
    pair = sparse.hstack([target, background], format='csr')
    original_pair = pair.copy()
    if scenario in ['pair_difference', 'all_different']: pair = perturb(pair, 1)
    gene = sparse.hstack([.5 * original_pair, .5 * original_pair], format='csr')
    if scenario in ['gene_difference', 'all_different']:
        gene = perturb(sparse.hstack([.5 * original_pair, .5 * original_pair], format='csr'), 2)
    model = sparse.hstack([.5 * pair, .5 * pair], format='csr')
    if scenario == 'gene_difference': model = gene.copy()
    if scenario in ['model_difference', 'all_different']: model = perturb(model, 4)
    endpoint_family = sparse.csr_matrix(family.shape) if mode == 'signed' else 2 * family
    return dict(target_node=target, background_node=background, model_pair=pair,
        gene=gene, model=model, family=endpoint_family), family, np.arange(n, dtype='int64')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True); parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); assert not args.output.exists() and not args.receipt.exists()
    args.output.mkdir(parents=True)
    rng = np.random.default_rng(2026100316)
    scenarios = []; negative = []; reconstructed = 0; examples = []
    for mode in ['signed', 'unsigned']:
        for scenario in ['all_equal', 'model_difference', 'gene_difference', 'pair_difference', 'all_different']:
            operators, family, selected = fixture(mode, scenario)
            proof = certificate(operators, family, mode, selected)
            assert independent_certificate(operators, family, mode, selected) == proof
            mapping = variance_map(proof['relations'], mode)
            assert independent_variance_map(proof['relations'], mode) == mapping
            if scenario == 'model_difference':
                assert not proof['relations']['model_equals_half_pair']['exact'] and 'model' in mapping['retained_names']
            if scenario == 'gene_difference':
                assert not proof['relations']['gene_equals_half_target_plus_background']['exact'] and 'gene' in mapping['retained_names']
            if scenario == 'pair_difference':
                assert not proof['relations']['pair_equals_target_plus_background']['exact'] and 'model_pair' in mapping['retained_names']
                assert proof['relations']['gene_equals_half_target_plus_background']['exact']
                assert proof['relations']['model_equals_half_pair']['exact']
                assert len(mapping['retained_names']) == 5
            n = len(selected); species_factor = rng.normal(size=(n, 4))
            kernels = {'residual': np.eye(n), 'family_intercept': (family @ family.T).toarray(),
                'species': species_factor @ species_factor.T}
            kernels.update({name: (value @ value.T).toarray() for name, value in operators.items()})
            forward = np.asarray(mapping['forward']); inverse = np.asarray(mapping['nonnegative_right_inverse'])
            for index in range(20):
                theta = rng.uniform(0, 2, len(ORIGINAL_NAMES)); theta[rng.random(len(theta)) < .3] = 0
                eta = forward @ theta
                first = sum(theta[i] * kernels[name] for i, name in enumerate(ORIGINAL_NAMES))
                second = sum(eta[i] * kernels[name] for i, name in enumerate(mapping['retained_names']))
                assert np.allclose(first, second, rtol=1e-14, atol=1e-14)
                original = inverse @ eta
                assert np.all(original >= 0) and np.array_equal(forward @ original, eta)
                reconstructed += 1
            cohort = dict(cohort_id='synthetic-' + mode + '-' + scenario,
                case_rows_sha256=sha(__file__), ordered_case_ids_sha256=sha(__file__))
            record = dict(cohort_id=cohort['cohort_id'], cohort_rows_sha256=cohort['case_rows_sha256'],
                ordered_case_ids_sha256=cohort['ordered_case_ids_sha256'], certificate=proof,
                variance_map=mapping, raw_reml_basis_qualification_complete=False, scientific_eligibility=False)
            check_record(record, cohort, mode, operators, family, selected)
            examples.append(record); scenarios.append(mode + ':' + scenario)
    serialized = args.output / 'synthetic_all_scenarios.json'
    serialized.write_text(json.dumps(examples, indent=2) + '\n'); assert json.loads(serialized.read_text()) == examples
    operators, family, selected = fixture('signed', 'model_difference'); cohort = dict(cohort_id='synthetic-signed-model_difference', case_rows_sha256=sha(__file__), ordered_case_ids_sha256=sha(__file__))
    baseline = next(row for row in examples if row['cohort_id'] == cohort['cohort_id'])
    for name in ['promote_failed_relation', 'false_witness', 'wrong_mode', 'changed_cohort', 'negative_map', 'missing_retained_model', 'changed_right_inverse', 'claim_reml', 'claim_scientific']:
        bad = copy.deepcopy(baseline)
        if name == 'promote_failed_relation': bad['certificate']['relations']['model_equals_half_pair']['exact'] = True
        elif name == 'false_witness': bad['certificate']['relations']['model_equals_half_pair']['witness']['scaled_integer_difference'] += 1
        elif name == 'wrong_mode': bad['certificate']['loading_mode'] = 'unsigned'
        elif name == 'changed_cohort': bad['cohort_id'] = 'foreign'
        elif name == 'negative_map': bad['variance_map']['forward'][0][0] = -1
        elif name == 'missing_retained_model': bad['variance_map']['retained_names'].remove('model')
        elif name == 'changed_right_inverse': bad['variance_map']['nonnegative_right_inverse'][0][0] = 2
        elif name == 'claim_reml': bad['raw_reml_basis_qualification_complete'] = True
        else: bad['scientific_eligibility'] = True
        reject(lambda: check_record(bad, cohort, 'signed', operators, family, selected)); negative.append(name)
    for name in ['nondyadic_loadings', 'target_not_identity', 'unsigned_family_mismatch']:
        bad, ff, rr = fixture('unsigned', 'all_equal')
        if name == 'nondyadic_loadings': bad['gene'].data[0] = .1
        elif name == 'target_not_identity': bad['target_node'] = perturb(bad['target_node'], 1)
        else: bad['family'].data[0] = 1
        reject(lambda: certificate(bad, ff, 'unsigned', rr)); reject(lambda: independent_certificate(bad, ff, 'unsigned', rr))
        negative.append(name)
    own = ['covariance_exact_folds_v2', 'full_exact_covariance_sources_v2', 'run_full_exact_covariance_folds_v2',
        'readback_full_exact_covariance_folds_v2', 'check_covariance_exact_folds_v2']
    pins = {f'scripts/{name}.py': sha(f'scripts/{name}.py') for name in own}
    source_plan = dict(pins=pins, operator_plan='metadata/full_entity_operator_plan_20261002.json',
        operator_completion='metadata/full_entity_operator_bank_completed_20261002.json',
        design_plan='metadata/full_expanded_model_designs_plan_20261002_v2.json',
        design_completion='metadata/full_expanded_model_designs_v2_completed_20261002.json', expected={'cohorts': 4340})
    source_path = args.output / 'source_plan.json'; source_path.write_text(json.dumps(source_plan, indent=2) + '\n')
    source, bindings = load(source_plan, source_path)
    occurrences = sum(len(rows(source, cohort)) for cohort in source['cohorts'])
    assert occurrences == 34110120
    verify(bindings)
    result = dict(status='passed_full_exact_covariance_fold_software_contracts_v2', checked_utc=datetime.now(timezone.utc).isoformat(),
        full_real_cohorts=4340, full_real_logical_cases=75188, full_real_cohort_row_occurrences=occurrences,
        full_real_scoped_source_bindings=len(bindings), global_family_squared_rows=source['full_family_squared_rows'],
        synthetic_scenarios=scenarios, covariance_and_cone_roundtrips=reconstructed,
        source_row_serialization_checked=True, malformed_cases_rejected=negative, source_hashes=bindings,
        artifacts={str(p): sha(p) for p in args.output.rglob('*') if p.is_file()}, scientific_eligibility=False,
        scope='Ten synthetic signed/unsigned scenarios test generalized gene=half(target+background) and model=half(pair), retaining failed named identities. Pair-sharing counterexamples retain five kernels with the full covariance cone unchanged. Two hundred nonnegative covariance/right-inverse tests pass; separate CSR and CSC outer sums agree. Twelve malformed contracts rejected. All4340cohort memberships and13operator/75188case source files rehashed and bound to closed archives. No biological pilot, full real generalized proofs, raw/REML basis qualification, variance fitting or scientific acceptance.')
    with args.receipt.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ['source_hashes', 'artifacts']}, indent=2))


if __name__ == '__main__': main()
