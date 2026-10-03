#!/usr/bin/env python3
"""Dense/score/optimizer contracts for both exact retained covariance models."""
import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha
from check_full_reduced_covariance_qualification import fixture, rejected
from check_full_shared_entity_fits import settings
from check_shared_entity_likelihood import dense
from full_expanded_model_design_sources import digest
from readback_retained_shared_entity_candidate import numeric
from reduced_covariance_basis import reduce_record
from retained_shared_entity_candidate import READY, candidate, backend_guard
from shared_entity_likelihood import SharedEntityLikelihood


def inputs(exception, mode):
    seed = 100 + int(exception)
    old, cert = fixture(exception, mode, seed)
    n = 40; labels = np.repeat(np.arange(8), 5)
    target = sparse.eye(n, format='csr')
    background = sparse.csr_matrix((np.ones(n), (np.arange(n), labels * 3 + (np.arange(n) % 5) // 2)),
                                   shape=(n, n))
    index = np.arange(n)
    if exception:
        index[1] = index[0]
    pair_target = sparse.csr_matrix((np.ones(n), (np.arange(n), index)), shape=(n, n))
    pair = sparse.hstack([pair_target, background], format='csr')
    intercept = sparse.csr_matrix((np.ones(n), (np.arange(n), labels)), shape=(n, 8))
    rng = np.random.default_rng(seed)
    factor = rng.normal(size=(n, 6)); x = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
    y = .6 * (background @ rng.normal(size=n)) + .3 * (intercept @ rng.normal(size=8)) + rng.normal(size=n)
    if exception:
        y += .4 * (pair @ rng.normal(size=pair.shape[1]))
    source = dict(labels=labels, factors={'mafft_guide': factor})
    saved = reduce_record(old, cert, 'explicit-synthetic-retained-fit-source')
    assert saved['disposition'] == READY
    bank = dict(background_node=background, model_pair=pair, family_intercept=intercept)
    operators = {k: bank[k] for k in saved['retained_kernel_names'][1:-1]}
    return source, np.arange(n), x, y, operators, saved, old, cert


def identity(audit, method, fit_disposition='ready_for_working_covariance_fit'):
    combined = audit['disposition'] if fit_disposition == 'ready_for_working_covariance_fit' else fit_disposition
    return dict(candidate_id=digest(['synthetic-retained-candidate', audit['audit_id'], method, fit_disposition]),
                covariance_audit_id=audit['audit_id'], covariance_audit_sha256=digest(audit),
                source_covariance_disposition=audit['disposition'], source_fit_disposition=fit_disposition,
                source_combined_disposition=combined, loading_mode=audit['loading_mode'], tree=audit['tree'],
                records=audit['records'], method=method, scientific_eligibility=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    args = p.parse_args(); args.output.mkdir(exist_ok=False)
    assert not args.receipt.exists()
    plan = settings(); plan['independent_backend'] = 'component_spectral_v1'
    fitted_counts = Counter(); read_counts = Counter(); dense_cases = derivative_cells = 0
    max_objective = max_gradient = max_derivative = 0.; rejects = []; fit_cases = []
    exclusions = 0; envelope_differences = 0
    for exception, mode in itertools.product([False, True], ['signed', 'unsigned']):
        source, rows, x, y, operators, audit, old, cert = inputs(exception, mode)
        snapshot = deepcopy(audit)
        guard = backend_guard(source, rows, x, operators, audit)
        assert audit == snapshot
        envelope_differences += int(guard['raw_roundoff_envelope'] != audit['numerical_audit']['raw_roundoff_envelope'])
        likelihood = SharedEntityLikelihood(source['labels'], operators, source['factors'][audit['tree']], np.ones(40), x, y)
        for method in ['ml', 'reml']:
            for scale in [0., .7, 100.]:
                ratios = np.full(len(operators) + 1, scale)
                value = likelihood.evaluate(ratios, method)
                reference = dense(source['labels'], operators, source['factors'][audit['tree']], np.ones(40), x, y, ratios, method)
                for key in reference:
                    np.testing.assert_allclose(value[key], reference[key], rtol=3e-9, atol=3e-10)
                max_objective = max(max_objective, abs(value['negative_profiled_likelihood'] - reference['negative_profiled_likelihood']))
                max_gradient = max(max_gradient, float(np.max(abs(value['gradient'] - reference['gradient']))))
                for j in range(len(ratios)):
                    h = 1e-4 * max(1., scale)
                    def evaluate(step):
                        changed = ratios.copy(); changed[j] += step * h
                        return likelihood.evaluate(changed, method, False)['negative_profiled_likelihood']
                    derivative = ((-25*evaluate(0)+48*evaluate(1)-36*evaluate(2)+16*evaluate(3)-3*evaluate(4))/(12*h)
                                  if scale == 0. else (evaluate(-2)-8*evaluate(-1)+8*evaluate(1)-evaluate(2))/(12*h))
                    np.testing.assert_allclose(value['gradient'][j], derivative, rtol=3e-5, atol=3e-7)
                    max_derivative = max(max_derivative, abs(value['gradient'][j] - derivative)); derivative_cells += 1
                dense_cases += 1
            record = identity(audit, method)
            fitted = candidate(source, plan, rows, record, x, y, operators, audit, old, cert)
            replay = numeric(source, plan, rows, record, x, y, fitted, operators, audit, old, cert)
            assert fitted['numerical_attempted'] is True
            assert fitted['fit'] is not None and 'variance_ratios' in fitted['fit'], fitted
            assert fitted['fit']['parameter_names'] == audit['retained_kernel_names'][1:]
            assert fitted['fit']['source_covariance_qualification'] == guard
            assert audit == snapshot
            fitted_counts[fitted['disposition']] += 1; read_counts[replay['disposition']] += 1
            fit_cases.append(dict(pair_exception=exception, loading_mode=mode, method=method,
                                  candidate=fitted, independent_audit=replay,
                                  inherited_audit=audit, original_audit=old, exact_certificate=cert))
            # Budgets remain enforced; failed searches are reproduced without promotion.
            exhausted = deepcopy(plan); exhausted['optimizer'].update(max_iterations=1, max_evaluations=1)
            failure = candidate(source, exhausted, rows, record, x, y, operators, audit, old, cert)
            replay_failure = numeric(source, exhausted, rows, record, x, y, failure, operators, audit, old, cert)
            assert replay_failure['disposition'] == 'original_numerical_failure_reproduced_requires_review'
            constant = identity(audit, method, 'constant_response_requires_review')
            excluded = candidate(source, plan, rows, constant, x, np.full(40, .125), operators, audit, old, cert)
            assert excluded['numerical_attempted'] is False
            assert numeric(source, plan, rows, constant, x, np.full(40, .125), excluded,
                           operators, audit, old, cert)['disposition'] == 'source_review_or_exclusion_retained'
            exclusions += 1
            for status in ['empty_design', 'rank_deficient_design', 'numerical_covariance_qualification_requires_review']:
                original = deepcopy(old); original.update(disposition=status, numerical_audit=None)
                nonready = reduce_record(original, cert, audit['retained_source_contract'])
                nonready_record = identity(nonready, method)
                exported = candidate(source, plan, rows, nonready_record, x, y, operators, nonready, original, cert)
                assert numeric(source, plan, rows, nonready_record, x, y, exported,
                               operators, nonready, original, cert)['disposition'] == 'source_review_or_exclusion_retained'
                assert exported['numerical_attempted'] is False; exclusions += 1
            for name in ['objective', 'variance', 'seed', 'backend_gram', 'invented_failure']:
                changed = deepcopy(fitted)
                if name == 'objective': changed['fit']['negative_profiled_likelihood'] += 1.
                elif name == 'variance': changed['fit']['variance_components'][0] += 1.
                elif name == 'seed': changed['fit']['starts'][1]['initial_coordinates'][0] += .1
                elif name == 'backend_gram': changed['fit']['source_covariance_qualification']['raw_gram'][0][0] += 1.
                else: changed.update(fit=None, disposition='shared_entity_fit_error_requires_review', error_type='ArithmeticError', error_message='invented')
                rejected(lambda: numeric(source, plan, rows, record, x, y, changed, operators, audit, old, cert))
                rejects.append([exception, mode, method, name])
        # Rehash altered source rows so semantic checks, not only hashes, reject them.
        for name in ['narrowed_inherited_error', 'changed_source_gram', 'lost_pair_or_background', 'promoted_attribution']:
            changed = deepcopy(audit)
            if name == 'narrowed_inherited_error': changed['numerical_audit']['raw_roundoff_envelope'][0][0] = 0.
            elif name == 'changed_source_gram': changed['numerical_audit']['raw_gram'][0][0] += 1.
            elif name == 'lost_pair_or_background': changed['retained_kernel_names'].pop(1)
            else: changed['component_variance_attribution_accepted'] = True
            rejected(lambda: candidate(source, plan, rows, identity(changed, 'reml'), x, y, operators, changed, old, cert))
            rejects.append([exception, mode, name])
    assert len(fit_cases) == 8 and dense_cases == 24 and derivative_cells == 84 and exclusions == 32
    assert len(rejects) == 56 and envelope_differences == 4
    # This gate validates the numerical bridge, not a full fitting grid or runtime estimate.
    artifact = args.output / 'retained_candidate_cases.json'
    with artifact.open('x') as f: json.dump(fit_cases, f, indent=2, allow_nan=False); f.write('\n')
    paths = [Path(__file__), Path('scripts/retained_shared_entity_candidate.py'),
             Path('scripts/readback_retained_shared_entity_candidate.py'), Path('scripts/fit_shared_entity_likelihood.py'),
             Path('scripts/readback_full_shared_entity_fits_v2.py'), Path('scripts/reduced_covariance_basis.py'),
             Path('scripts/check_full_reduced_covariance_qualification.py'), Path('scripts/check_shared_entity_likelihood.py'),
             Path('scripts/shared_entity_likelihood.py'), Path('scripts/independent_shared_entity_likelihood_fast.py'),
             Path('scripts/independent_shared_entity_likelihood.py'), Path('scripts/independent_shared_entity_optimizer.py'),
             Path('scripts/covariance_basis_independent.py'), Path('scripts/readback_full_covariance_qualification.py')]
    result = dict(status='passed_retained_covariance_candidate_dense_score_and_independent_fit_contracts',
                  checked_utc=datetime.now(timezone.utc).isoformat(), dense_ml_reml_cases=dense_cases,
                  finite_difference_score_cells=derivative_cells, numerical_candidate_cases=len(fit_cases),
                  candidate_status_counts=dict(fitted_counts), independent_status_counts=dict(read_counts),
                  source_reviews_or_constants_retained=exclusions, altered_cases_rejected=rejects,
                  inherited_and_fresh_envelopes_differing_scenarios=envelope_differences,
                  inherited_envelopes_never_replaced=True, strict_source_and_additional_backend_guards_required=True,
                  maximum_absolute_dense_objective_error=max_objective, maximum_absolute_dense_gradient_error=max_gradient,
                  maximum_absolute_finite_difference_error=max_derivative, source_hashes={str(p): sha(p) for p in paths},
                  artifacts={str(artifact): sha(artifact)}, scientific_eligibility=False,
                  full_grid_fitting_implemented=False, production_fitting_launched=False,
                  scope='Four/five exact retained kernels, signed/unsigned and ML/REML; eight actual synthetic '
                  'multi-start candidates with independent latent/spectral/start/curvature/search audits as applicable. '
                  'Twenty-four dense covariance/score cases including zero boundaries,84 derivative cells, eight '
                  'exhausted-budget replays,32 source/constant exclusions and56 altered cases. Original principal '
                  'error envelopes and guards retained separately from fresh backend guards. No fungal pilot, '
                  'closed production source proof, full-grid fitting workflow, runtime estimate or biological acceptance.')
    with args.receipt.open('x') as f: json.dump(result, f, indent=2, allow_nan=False); f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts', 'altered_cases_rejected']}, indent=2))


if __name__ == '__main__':
    main()
