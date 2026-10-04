#!/usr/bin/env python3
"""Exercise actual refits and serialized independent readback over all 64 axes."""
import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha
from check_retained_shared_entity_candidates import inputs
from check_shared_entity_likelihood import dense
from check_weighted_shared_entity_simulation import model_for, rejected
from full_expanded_model_design_sources import digest
from full_weighted_fit_admission import FULL
from full_weighted_fit_exports import atomic
from reference_measurement_union_sources import bind, verify
from refit_weighted_shared_entity_simulation import GaussianRefit, original_identity
from weighted_shared_entity_simulation import coverage_counts

PLAN = 'metadata/full_weighted_shared_entity_fit_draft_plan_20261004_v1.json'
PLAN_HASH = '475a58c6abd9ca7e85120980f6f34718d389de28eca341dce9c6c094cb1ef36f'


def context_arguments(case):
    source, rows, x, y, operators, *_ = inputs(case['pair_exception'], case['loading_mode'])
    if case['policy'] != 'uniform':operators = {'target_node': sparse.eye(len(rows), format='csr'), **operators}
    if case['outcome'] == 'native_tm_dissimilarity_delta':y = np.tanh(y)+.15*x[:, 1]
    identity = original_identity(case['candidate']);audit = deepcopy(case['audit'])
    route = dict(names=audit['retained_kernel_names'], certificate_sha256=audit['exact_certificate_sha256'],
        route=audit['basis_route'], exact_uniform_one=audit['residual_diagonal_is_exact_uniform_one'])
    return dict(plan_path=PLAN, expected_plan_sha256=PLAN_HASH, source=source, rows=rows,
        identity=identity, expected_identity_sha256=digest(identity), design=x, original_response=y,
        operators=operators, audit=audit, route=route, diagonal=case['diagonal'],
        expected_model_contract=model_for(case).input_contract)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True);p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args();assert not a.receipt.exists();a.output.mkdir(parents=True, exist_ok=False)
    bindings = {}
    for parent in ['metadata/weighted_shared_entity_simulation_software_transport_20261004_v1.json',
                   'metadata/weighted_shared_entity_simulation_original_terminal_payload_20261004_v1.json']:
        saved = json.loads(Path(parent).read_text());verify(saved['source_hashes'])
        for name, d in saved['source_hashes'].items():bind(bindings, name, d)
        bind(bindings, parent)
    case_path = Path('data/software_audits/weighted-shared-entity-candidates-20261004-v1/actual_weighted_candidate_cases.json')
    cases = json.loads(case_path.read_text());assert len(cases) == 64;bind(bindings, case_path)
    assert sha(PLAN) == PLAN_HASH
    full = json.loads(Path(PLAN).read_text());verify(full['pins']);assert full['expected'] == FULL;bind(bindings, PLAN)
    for name in ['refit_weighted_shared_entity_simulation', 'check_weighted_shared_entity_simulation_refits',
                 'full_weighted_fit_admission', 'check_shared_entity_likelihood']:
        bind(bindings, 'scripts/'+name+'.py')
    # Shape/hash/semantic alterations operate on private copies, never parents.
    baseline = context_arguments(cases[0]);source_rejections = []
    for name in ['plan_hash', 'identity_hash', 'response_hash', 'model_contract', 'rotated_design',
                 'row_order', 'repeated_rows', 'fractional_rows', 'outside_rows', 'coefficient_mapping',
                 'science_promotion', 'nonuniform_promotion', 'source_gram', 'route', 'diagonal']:
        args = deepcopy(baseline)
        if name == 'plan_hash':args['expected_plan_sha256'] = 'foreign'
        elif name == 'identity_hash':args['expected_identity_sha256'] = 'foreign'
        elif name == 'response_hash':args['original_response'][0] += .1
        elif name == 'model_contract':args['expected_model_contract'] = 'foreign'
        elif name == 'rotated_design':args['design'][:, [1, 2]] = args['design'][:, [2, 1]]
        elif name == 'row_order':args['rows'] = args['rows'][::-1]
        elif name == 'repeated_rows':args['rows'][0] = args['rows'][1]
        elif name == 'fractional_rows':args['rows'] = args['rows'].astype(float)
        elif name == 'outside_rows':args['rows'][0] = len(args['source']['labels'])
        elif name == 'coefficient_mapping':args['identity']['coefficient_columns'] = ['one']
        elif name == 'science_promotion':args['identity']['scientific_eligibility'] = True
        elif name == 'nonuniform_promotion':args['identity']['nonuniform_weighting_accepted'] = True
        elif name == 'source_gram':
            args['audit']['numerical_audit']['raw_gram'][0][0] += 1
            args['identity']['covariance_audit_sha256'] = digest(args['audit'])
        elif name == 'route':args['route']['names'] = args['route']['names'][::-1]
        elif name == 'diagonal':args['diagonal'][0] += .1
        if name in ['coefficient_mapping', 'science_promotion', 'nonuniform_promotion', 'source_gram']:
            args['expected_identity_sha256'] = digest(args['identity'])
        # Source Gram is checked by the frozen production backend on refit.
        if name == 'source_gram':
            def operation():
                changed = GaussianRefit(**args)
                y, g, i = changed.generate([.2, -.3, .4], 1., np.zeros(len(changed.model.parameter_names)), 20261004, 'bad', 0)
                return changed.produce(y, g, i)
            rejected(operation)
        else:rejected(lambda: GaussianRefit(**args))
        source_rejections.append(name)
    primary_statuses = Counter();readback_statuses = Counter();refit_statuses = Counter()
    basis_dimensions = set();maximum_dense_objective_error = maximum_dense_beta_error = 0.
    record_paths = [];contexts = [];generating_truth = np.array([.2, -.3, .4]);fits = 0
    for number, case in enumerate(cases):
        args = context_arguments(case);context = GaussianRefit(**args);contexts.append(context)
        basis_dimensions.add(len(context.model.parameter_names)+1)
        directory = a.output/('case-'+str(number));directory.mkdir()
        atomic(directory/'source.json', dict(identity=context.identity, audit=context.audit, route=context.route,
            ordered_original_rows_sha256=context.original_rows_sha256, upstream_model_contract=context.model.input_contract,
            software_fixture_only=True))
        np.savez_compressed(directory/'inputs.npz', design=context.model.design, diagonal=context.model.diagonal,
            factor=context.model.factor, original_response=args['original_response'])
        for scenario in ['zero_components', 'positive_components']:
            ratios = np.zeros(len(context.model.parameter_names)) if scenario == 'zero_components' else np.full(len(context.model.parameter_names), .7)
            y, generation, identity = context.generate(generating_truth, 1., ratios, 20261004, scenario, 0)
            prefix = directory/scenario;prefix.mkdir()
            np.save(prefix/'response.npy', y)
            atomic(prefix/'generation.json', generation);atomic(prefix/'simulated_identity.json', identity)
            primary = context.produce(y, generation, identity);fits += 1
            atomic(prefix/'primary.json', primary)
            record = context.read(y, generation, identity, primary)
            atomic(prefix/'refit.json', record);record_paths.append(prefix/'refit.json')
            # Fresh serialized readback uses the frozen independent math. The
            # producer is rerun only if that reader needs failure reproduction.
            serialized = json.loads((prefix/'refit.json').read_text())
            assert context.replay(serialized) == serialized
            summary = context.coverage([serialized], 1);atomic(prefix/'accounting.json', summary)
            assert summary['attempted'] == 1 and summary['unresolved']+summary['independently_checked'] == 1
            primary_statuses[primary['disposition']] += 1
            readback_statuses[record['independent']['disposition']] += 1;refit_statuses[record['status']] += 1
            fit = primary['fit']
            if fit is not None and 'variance_ratios' in fit:
                reference = dense(context.source['labels'], context.model.incidence, context.model.factor,
                    context.model.diagonal, context.model.design, y, np.asarray(fit['variance_ratios']), case['method'])
                for key in ['negative_profiled_likelihood', 'beta', 'conditional_beta_covariance', 'profiled_scale']:
                    np.testing.assert_allclose(fit[key], reference[key], rtol=3e-9, atol=3e-10)
                maximum_dense_objective_error = max(maximum_dense_objective_error, abs(fit['negative_profiled_likelihood']-reference['negative_profiled_likelihood']))
                maximum_dense_beta_error = max(maximum_dense_beta_error, float(np.max(abs(np.asarray(fit['beta'])-reference['beta']))))
            print('actual_gaussian_refits', fits, flush=True)
    assert fits == 128 and len(contexts) == 64 and basis_dimensions == {4, 5, 6}
    axes = {(c['pair_exception'], c['loading_mode'], c['policy'], c['outcome'], c['method']) for c in cases}
    assert len(axes) == 64
    record = json.loads(record_paths[0].read_text());context = contexts[0];replay_rejections = []
    for name in ['generation_seed', 'generation_truth', 'response_hash', 'calibration_contract', 'interval_reference',
                 'method', 'primary_objective', 'primary_start', 'primary_identity', 'primary_budget',
                 'independent_disposition', 'refit_status', 'lost_review', 'science_promotion', 'invented_failure']:
        changed = deepcopy(record)
        if name == 'generation_seed':changed['generation']['master_seed'] += 1
        elif name == 'generation_truth':changed['generation']['scenario']['beta'][0] += .1
        elif name == 'response_hash':changed['generation']['response_sha256'] = 'foreign'
        elif name == 'calibration_contract':changed['calibration_contract'] = 'foreign'
        elif name == 'interval_reference':changed['calibration_configuration']['interval_definition']['reference'] = 't'
        elif name == 'method':changed['calibration_configuration']['method'] = 'reml'
        elif name == 'primary_objective':changed['primary']['fit']['negative_profiled_likelihood'] += 1
        elif name == 'primary_start':changed['primary']['fit']['starts'][1]['initial_coordinates'][0] += .1
        elif name == 'primary_identity':changed['primary']['original_response_sha256'] = 'foreign'
        elif name == 'primary_budget':changed['primary']['fit']['starts'][0]['search_evaluations'] = 501
        elif name == 'independent_disposition':changed['independent']['disposition'] = 'invented'
        elif name == 'refit_status':changed['status'] = 'invented'
        elif name == 'lost_review':changed['independent'] = {}
        elif name == 'science_promotion':changed['scientific_eligibility'] = True
        else:changed['primary'].update(fit=None, disposition='shared_entity_fit_error_requires_review', error_type='ValueError', error_message='invented')
        rejected(lambda: context.replay(changed));replay_rejections.append(name)
    checked = next((json.loads(p.read_text()) for p in record_paths
        if json.loads(p.read_text())['status'] == 'refit_independently_checked_pending_calibration'), None)
    interval_rejections = []
    if checked is not None:
        actual_context = next(c for c in contexts if c.calibration_contract == checked['calibration_contract'])
        for name in ['coverage', 'standard_error', 'interval_lower', 'generating_beta']:
            changed = deepcopy(checked)
            if name == 'coverage':changed['nominal_interval_covers_generating_beta'][0] = not changed['nominal_interval_covers_generating_beta'][0]
            else:changed['nominal_interval_diagnostics'][{'standard_error':'standard_errors', 'interval_lower':'nominal_lower', 'generating_beta':'generating_beta'}[name]][0] += .1
            rejected(lambda: actual_context.replay(changed));interval_rejections.append(name)
    # Explicitly injected software-only exception, never counted as actual
    # numerical calibration. Assertions/invariant corruption remain fatal.
    y = np.load(record_paths[0].parent/'response.npy')
    with patch('refit_weighted_shared_entity_simulation.numeric', side_effect=ValueError('injected software-only reader exception')):
        error_record = context.read(y, record['generation'], record['simulated_identity'], record['primary'])
    assert error_record['status'] == 'refit_error_requires_review' and 'nominal_interval_covers_generating_beta' not in error_record
    error_summary = coverage_counts([error_record], 3, 1)
    assert error_summary['unresolved'] == error_summary['attempted'] == 1
    rejected(lambda: context.replay(error_record))
    atomic(a.output/'injected-error-accounting.json', dict(record=error_record, accounting=error_summary,
        injected_software_branch_only=True, actual_numerical_calibration=False))
    generating_rejections = []
    for name, ratios in [('outside_scaled_box', np.full(len(context.model.parameter_names), 1e20)),
                         ('negative_ratio', np.full(len(context.model.parameter_names), -1.))]:
        rejected(lambda: context.generate(generating_truth, 1., ratios, 20261004, name, 0));generating_rejections.append(name)
    accounting_rejections = []
    for name, records, replicates in [('missing_replicate', [record], 2), ('duplicate_replicate', [record, record], 2),
                                    ('foreign_context', [json.loads(record_paths[-1].read_text())], 1)]:
        rejected(lambda: context.coverage(records, replicates));accounting_rejections.append(name)
    for path in a.output.rglob('*'):
        if path.is_file():bind(bindings, path)
    verify(bindings)
    result = dict(status='passed_actual_weighted_gaussian_refits_and_frozen_readback_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(), original_software_model_cases=64,
        original_axis_combinations=64, basis_dimensions=sorted(basis_dimensions), generating_scenarios=2,
        actual_software_refits=128, serialized_independent_refit_replays=128,
        primary_status_counts=dict(primary_statuses), independent_status_counts=dict(readback_statuses),
        refit_status_counts=dict(refit_statuses), dense_selected_point_comparisons=128,
        maximum_absolute_dense_objective_error=maximum_dense_objective_error,
        maximum_absolute_dense_beta_error=maximum_dense_beta_error,
        source_rejections=source_rejections, replay_rejections=replay_rejections,
        interval_rejections=interval_rejections, generating_rejections=generating_rejections,
        accounting_rejections=accounting_rejections, injected_reader_error_branch_retained=True,
        interval_definition=context.interval, original_full_expected=FULL,
        full_production_source_admitted=False, full_data_simulation_launched=False,
        empirical_coverage_accepted=False, scientific_eligibility=False,
        production_biological_fits_computed=0, gpu=False, new_cost_usd=0, source_hashes=bindings,
        scope='128 actual unmocked frozen mathematical refits over 64 synthetic original model axes, '
            'with serialized independent checks, selected-point dense references and all reviews retained. '
            'Single replicate per model/scenario exercises plumbing, not coverage calibration. Exception '
            'injection is a separate explicit software fixture. No full source admission, biological '
            'calibration, global joint-null/multiple-testing analysis or aim completion claimed.')
    atomic(a.receipt, result)
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__':main()
