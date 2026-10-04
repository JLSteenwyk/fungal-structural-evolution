"""Complete original designs/responses bound to closed four-control audits.

Every original cohort, design, response, policy, mode and tree remains present.
Parallel producer and reader checkpoints are additional mandatory artifacts.
Scoped source hashes and the full consumed numerical exports are verified;
unconsumed historical parent archives are not replayed as new science.
"""
import itertools
import json
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import closed_subset
from full_expanded_model_design_sources import AXES, DEGREES, OUTCOMES, array_digest, digest, fit_id
from full_expanded_model_input_sources import ORDERS
from full_weighted_covariance_qualification import PRODUCER, READER, SUMMARY, expected_record, source_inputs
from full_weighted_covariance_sources_v2 import MODES, POLICIES, basis, cohort_controls, design_matrix, jsonl, retained_operators
from reference_measurement_union_sources import bind, verify

SCHEMA = 'full-four-control-shared-entity-fit-v1'


def load(plan, path):
    assert plan['methods'] == ['ml', 'reml'] and plan['loading_modes'] == MODES and plan['policies'] == POLICIES
    qp = Path(plan['qualification_plan']); qplan = json.loads(qp.read_text())
    assert qplan['trees'] == plan['trees']
    source, bindings = source_inputs(qplan, qp); original = dict(bindings)
    root = Path(qplan['output']); rp = root / 'receipt.json'; rb = root / 'readback.json'
    producer = json.loads(rp.read_text()); reader = json.loads(rb.read_text())
    assert producer['status'] == PRODUCER and reader['status'] == READER
    assert producer['plan_sha256'] == reader['plan_sha256'] == sha(qp)
    assert producer['source_contract'] == reader['source_contract'] == source['contract']
    assert producer['source_hashes'] == original
    assert reader['producer_receipt_sha256'] == sha(rp)
    assert reader['artifacts'] == producer['artifacts']
    assert producer['working_model_fits_computed'] == reader['working_model_fits_computed'] == 0
    assert producer['scientific_eligibility'] is reader['scientific_eligibility'] is False
    manifest = json.loads((root / 'cohort_manifest.json').read_text())
    assert [m['cohort_id'] for m in manifest] == [c['cohort_id'] for c in source['cohorts']]
    names = {'stage_plan.json', 'cohort_manifest.json'} | {
        m[k] for m in manifest for k in ['audit_file', 'link_file']}
    from weighted_parallel_checkpoint_contracts_v1 import checkpoint_artifacts
    extra = checkpoint_artifacts(root, manifest, producer, reader)
    names.update(extra)
    assert set(producer['artifacts']) == names
    reader_names = set(reader['reader_checkpoint_artifacts'])
    completion = closed_subset(plan['qualification_completion'],
        'complete_verified_full_four_control_covariance_numerical_qualification_v1',
        [qp, rp, rb, *[root / n for n in sorted(names | reader_names)]], bindings)
    assert completion['producer_receipt'] == str(rp) and completion['producer_receipt_sha256'] == sha(rp)
    assert completion['independent_readback'] == str(rb) and completion['independent_readback_sha256'] == sha(rb)
    for k in SUMMARY: assert producer[k] == reader[k] == completion[k], k
    assert json.loads((root / 'stage_plan.json').read_text()) == dict(
        schema='full-four-control-covariance-numerical-qualification-v1', plan_sha256=sha(qp), source_contract=source['contract'])
    for m in manifest:
        for key, h in [('audit_file', 'audit_sha256'), ('link_file', 'link_sha256')]:
            assert bindings[str(root / m[key])] == m[h] == producer['artifacts'][m[key]]
    for p, h in plan['pins'].items(): bind(bindings, p, h)
    bind(bindings, path)
    assert plan['expected']['candidate_rows'] == completion['designs'] * 2 * 2 * 2 * 5 * 4
    assert plan['expected']['setting_fit_links'] == completion['settings'] * 2 * 2 * 5 * 4
    source.update(numerical_root=root, numerical_manifest=manifest, numerical_completion=completion)
    source['fit_contract'] = digest(dict(schema=SCHEMA,
        numerical_contract=source['contract'], numerical_completion_sha256=sha(plan['qualification_completion']),
        source_hashes=bindings, methods=plan['methods'], optimizer=plan['optimizer'],
        independent_audit=plan['independent_audit'], independent_backend=plan['independent_backend']))
    verify(bindings)
    return source, bindings


def cohorts(source, plan):
    designs = jsonl(source['root'] / 'unique_designs.jsonl')
    fits = jsonl(source['root'] / 'unique_fit_inputs.jsonl')
    for cohort, control, manifest in zip(source['cohorts'], source['controls'], source['numerical_manifest']):
        rows, diagonals, record = cohort_controls(source, cohort, control)
        entries = []; seen = set()
        for _ in range(30):
            d = next(designs); assert d['cohort_id'] == cohort['cohort_id']
            seen.add((d['order_contrast'], d['sequence_axis'], d['degree']))
            matrix = design_matrix(source, cohort, rows, d); responses = {}
            for outcome in OUTCOMES:
                f = next(fits); y = source['arrays'][cohort['mask'], d['order_contrast']][outcome][rows]
                assert f['cohort_id'] == cohort['cohort_id'] and f['design_id'] == d['design_id'] and f['outcome'] == outcome
                assert np.isfinite(y).all() and f['records'] == len(rows)
                assert f['trees'] == plan['trees']
                assert f['response_min'] == float(y.min()) and f['response_max'] == float(y.max())
                assert f['response_sha256'] == array_digest(y, '<f8')
                assert f['fit_input_id'] == fit_id(d['design_id'], outcome, f['response_sha256'])
                expected = (d['disposition'] if d['disposition'] != 'full_rank_design' else
                    'constant_response_requires_review' if len(set(y.tolist())) == 1 else 'ready_for_working_covariance_fit')
                assert f['disposition'] == expected
                responses[outcome] = f, y
            entries.append((d, matrix, responses, {}))
        assert seen == set(itertools.product(ORDERS, AXES, DEGREES))
        exports = jsonl(source['numerical_root'] / manifest['audit_file'])
        for mode, tree in itertools.product(MODES, plan['trees']):
            routes = [basis(source, cohort, mode, diagonal) for diagonal in diagonals]
            for d, matrix, responses, audits in entries:
                for j, policy in enumerate(POLICIES):
                    saved = next(exports)
                    expected = expected_record(source, cohort, d, routes[j], diagonals[j], control, mode, tree, policy)
                    assert all(saved[k] == v for k, v in expected.items())
                    audits[mode, tree, policy] = saved, routes[j], diagonals[j]
        assert next(exports, None) is None
        yield cohort, rows, entries
    assert next(designs, None) is None and next(fits, None) is None


def identity(source, cohort, design, fit, audit, mode, tree, policy, method):
    combined = fit['disposition'] if fit['disposition'] != 'ready_for_working_covariance_fit' else audit['disposition']
    record = dict(candidate_id=digest([SCHEMA, source['fit_contract'], fit['fit_input_id'], mode, tree, policy, method]),
        source_contract=source['fit_contract'], cohort_id=cohort['cohort_id'], design_id=design['design_id'],
        fit_input_id=fit['fit_input_id'], outcome=fit['outcome'], loading_mode=mode, tree=tree,
        control_policy=policy, method=method, covariance_audit_id=audit['audit_id'], covariance_audit_sha256=digest(audit),
        source_fit_disposition=fit['disposition'], source_covariance_disposition=audit['disposition'],
        source_combined_disposition=combined, records=fit['records'], response_sha256=fit['response_sha256'],
        raw_design_sha256=design['raw_design_sha256'], ordered_case_ids_sha256=cohort['ordered_case_ids_sha256'],
        active_column_indices=design['active_column_indices'], exactly_zero_column_indices=design['exactly_zero_column_indices'],
        coefficient_columns=[design['predictor_columns'][i] for i in design['active_column_indices']],
        scientific_eligibility=False, nonuniform_weighting_accepted=False, component_variance_attribution_accepted=False)
    record.update({k: audit[k] for k in ['diagonal_sha256', 'control_record_sha256',
        'exact_certificate_sha256', 'retained_kernel_names']})
    return record


def cases(source, plan, cohort, entries):
    for design, matrix, responses, audits in entries:
        for mode, tree, policy in itertools.product(MODES, plan['trees'], POLICIES):
            audit, route, diagonal = audits[mode, tree, policy]
            for outcome, method in itertools.product(OUTCOMES, plan['methods']):
                fit, y = responses[outcome]
                record = identity(source, cohort, design, fit, audit, mode, tree, policy, method)
                yield record, matrix, y, audit, route, diagonal


def operators_for(source, rows, audit):
    return retained_operators(source, rows, audit['loading_mode'], audit['retained_kernel_names'])
