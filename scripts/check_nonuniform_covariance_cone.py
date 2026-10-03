#!/usr/bin/env python3
"""Dense cone roundtrips and nonuniform ML/REML/score software qualification."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha
from check_retained_shared_entity_candidates import inputs
from check_shared_entity_likelihood import dense
from covariance_exact_folds_v2 import ORIGINAL_NAMES
from independent_nonuniform_covariance_cone import mapping as independent_mapping, readback
from independent_shared_entity_likelihood_fast import ComponentSpectralLikelihood
from nonuniform_covariance_cone import mapping, record, validate_diagonal
from shared_entity_likelihood import SharedEntityLikelihood


def rejected(action):
    try:
        action()
    except (AssertionError, ValueError, KeyError):
        return
    raise AssertionError('Malformed nonuniform cone certificate accepted')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); args.output.mkdir(exist_ok=False)
    assert not args.receipt.exists()
    rng = np.random.default_rng(202610034)
    roundtrips = scores = derivatives = identity_rejections = 0
    max_covariance = max_objective = max_gradient = max_derivative = 0.
    exports = []
    alterations = ['fold_target_into_residual', 'lost_exception_pair', 'wrong_gene_coefficient',
        'wrong_family_fold', 'negative_map', 'wrong_right_inverse', 'promote_weighting',
        'claim_numerical_qualification', 'invent_diagonal', 'change_parent_sha',
        'change_rows_sha', 'change_mode', 'drop_target', 'promote_science']
    for exception, mode in itertools.product([False, True], ['signed', 'unsigned']):
        source, rows, x, y, ops, audit, old, parent = inputs(exception, mode)
        saved = record(parent, 'explicit-synthetic-nonuniform-cone')
        assert mapping(parent) == independent_mapping(parent); readback(saved, parent, saved['source_contract'])
        names = saved['variance_map']['retained_names']; n = len(rows)
        target = sparse.eye(n, format='csr')
        if exception:
            pair = ops['model_pair']
        else:
            pair = sparse.hstack([target, ops['background_node']], format='csr')
        bank = dict(target_node=target, background_node=ops['background_node'],
                    model_pair=pair, family_intercept=ops['family_intercept'])
        selected = {name: bank[name] for name in names[1:-1]}
        factor = source['factors'][audit['tree']]
        identity = np.eye(n); b = (bank['background_node'] @ bank['background_node'].T).toarray()
        p = (pair @ pair.T).toarray(); family = (bank['family_intercept'] @ bank['family_intercept'].T).toarray()
        patterns = [np.ones(n), np.full(n, 2.), np.where(np.arange(n) % 2, .5, 2.),
                    np.linspace(.125, 3., n), 1. + np.arange(n) * 2.**-40]
        forward = np.asarray(saved['variance_map']['forward'])
        inverse = np.asarray(saved['variance_map']['nonnegative_right_inverse'])
        for d in patterns:
            validate_diagonal(d, n)
            original_kernels = [np.diag(d), identity, b, p, .5 * (identity + b), .5 * p,
                                np.zeros((n, n)) if mode == 'signed' else 4 * family,
                                family, factor @ factor.T]
            retained_kernels = [original_kernels[ORIGINAL_NAMES.index(name)] for name in names]
            for v in [*np.eye(9), *rng.uniform(0., 5., (8, 9))]:
                expected = sum(amount * k for amount, k in zip(v, original_kernels))
                actual = sum(amount * k for amount, k in zip(forward @ v, retained_kernels))
                np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)
                max_covariance = max(max_covariance, float(np.max(abs(actual - expected))))
                roundtrips += 1
            for v in [*np.eye(len(names)), *rng.uniform(0., 5., (8, len(names)))]:
                np.testing.assert_array_equal(forward @ (inverse @ v), v)
                expected = sum(amount * k for amount, k in zip(v, retained_kernels))
                actual = sum(amount * k for amount, k in zip(inverse @ v, original_kernels))
                np.testing.assert_array_equal(actual, expected); roundtrips += 1
            primary = SharedEntityLikelihood(source['labels'], selected, factor, d, x, y)
            spectral = ComponentSpectralLikelihood(source['labels'], selected, factor, d, x, y, column_batch=3)
            for method, scale in itertools.product(['ml', 'reml'], [0., .7, 100.]):
                ratios = np.full(len(names) - 1, scale)
                actual = primary.evaluate(ratios, method)
                independent = spectral.evaluate(ratios, method)
                reference = dense(source['labels'], selected, factor, d, x, y, ratios, method)
                for key in reference:
                    np.testing.assert_allclose(actual[key], reference[key], rtol=2e-8, atol=2e-9)
                    np.testing.assert_allclose(independent[key], reference[key], rtol=2e-8, atol=2e-9)
                max_objective = max(max_objective, abs(actual['negative_profiled_likelihood'] - reference['negative_profiled_likelihood']))
                max_gradient = max(max_gradient, float(np.max(abs(actual['gradient'] - reference['gradient']))))
                permutation = rng.permutation(n)
                permuted = SharedEntityLikelihood(source['labels'][permutation],
                    {name: value[permutation] for name, value in selected.items()},
                    factor[permutation], d[permutation], x[permutation], y[permutation]).evaluate(ratios, method)
                for key in reference:
                    np.testing.assert_allclose(permuted[key], reference[key], rtol=2e-8, atol=2e-9)
                if scale != 100.:
                    for j in range(len(ratios)):
                        h = 1e-4
                        def point(distance):
                            varied = ratios.copy(); varied[j] += distance * h
                            return primary.evaluate(varied, method, False)['negative_profiled_likelihood']
                        numerical = ((-25 * point(0) + 48 * point(1) - 36 * point(2) + 16 * point(3) - 3 * point(4)) / (12 * h)
                                     if scale == 0. else (point(-2) - 8 * point(-1) + 8 * point(1) - point(2)) / (12 * h))
                        np.testing.assert_allclose(numerical, actual['gradient'][j], rtol=3e-5, atol=3e-7)
                        max_derivative = max(max_derivative, abs(numerical - actual['gradient'][j])); derivatives += 1
                scores += 1
            exports.append(dict(certificate=saved, diagonal=d.tolist(), pair_exception=exception, loading_mode=mode))
        for name in alterations:
            changed = deepcopy(saved); m = changed['variance_map']
            if name == 'fold_target_into_residual':
                m['forward'][0][1] = 1.; m['forward'][1][1] = 0.
            elif name == 'lost_exception_pair':
                m['retained_names'] = [k for k in names if k != ('model_pair' if exception else 'background_node')]
            elif name == 'wrong_gene_coefficient': m['forward'][1][4] = .25
            elif name == 'wrong_family_fold': m['forward'][-2][6] = 1.
            elif name == 'negative_map': m['forward'][1][4] = -.5
            elif name == 'wrong_right_inverse': m['nonnegative_right_inverse'][1][1] = 0.
            elif name == 'promote_weighting': m['nonuniform_weighting_accepted'] = True
            elif name == 'claim_numerical_qualification': m['raw_reml_basis_qualification_complete'] = True
            elif name == 'invent_diagonal': m['residual_diagonal_prepared'] = True
            elif name == 'change_parent_sha': changed['parent_exact_certificate_sha256'] = 'foreign'
            elif name == 'change_rows_sha': changed['cohort_rows_sha256'] = 'foreign'
            elif name == 'change_mode': changed['loading_mode'] = 'foreign'
            elif name == 'drop_target': m['retained_names'].remove('target_node')
            else: changed['scientific_eligibility'] = True
            rejected(lambda: readback(changed, parent, saved['source_contract'])); identity_rejections += 1
        for diagonal in [np.zeros(n), np.full(n, -1.), np.full(n, np.nan), np.full(n, np.inf), np.ones(n + 1)]:
            rejected(lambda: validate_diagonal(diagonal, n))
    assert (len(exports), roundtrips, scores, derivatives, identity_rejections) == (20, 610, 120, 360, 56)
    artifact = args.output / 'nonuniform_cone_cases.json'; artifact.write_text(json.dumps(exports, indent=2) + '\n')
    modules = ['check_nonuniform_covariance_cone', 'nonuniform_covariance_cone',
        'independent_nonuniform_covariance_cone', 'reduced_covariance_basis', 'covariance_exact_folds_v2',
        'check_retained_shared_entity_candidates', 'check_shared_entity_likelihood',
        'shared_entity_likelihood', 'shared_entity_covariance', 'independent_shared_entity_likelihood_fast',
        'independent_shared_entity_likelihood']
    result = dict(status='passed_positive_diagonal_covariance_cone_dense_likelihood_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(), synthetic_diagonal_cases=20,
        covariance_cone_roundtrips=610, ml_reml_dense_spectral_permutation_cases=120,
        finite_difference_score_cells=360, malformed_certificate_cases_rejected=56,
        invalid_diagonal_cases_rejected=20, maximum_covariance_error=max_covariance,
        maximum_dense_objective_error=max_objective, maximum_dense_gradient_error=max_gradient,
        maximum_finite_difference_error=max_derivative,
        source_hashes={'scripts/' + name + '.py': sha('scripts/' + name + '.py') for name in modules},
        artifacts={str(artifact): sha(artifact)}, scientific_eligibility=False,
        scope='Synthetic exact cone and likelihood/score checks for both modes and both pair identities, '
              'five diagonal patterns including constant, heterogeneous and near-uniform values. '
              'Target I stays distinct from residual D, including intentionally dependent constant cases. '
              'No weighted source preparation, raw/REML identifiability qualification, optimizer fit, '
              'confidence-to-variance calibration or biological acceptance. Full real certificate proof '
              'and actual weighted design/provenance/qualification/timing remain separate stages.')
    with args.receipt.open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}))


if __name__ == '__main__':
    main()
