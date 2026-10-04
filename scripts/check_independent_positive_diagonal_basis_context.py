#!/usr/bin/env python3
"""Qualify all cached latent D calculations against separate and dense results."""
import argparse
from datetime import datetime, timezone
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha
from check_positive_diagonal_kernel_products import dense, rejected
from check_retained_shared_entity_candidates import inputs
from independent_positive_diagonal_kernel_products import PositiveDiagonalKernelProducts
from independent_positive_diagonal_basis_context import IndependentPositiveDiagonalBasisContext
from positive_diagonal_basis_context import PositiveDiagonalBasisContext
from reduced_covariance_basis import independent_diagnostics


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True); a = p.parse_args()
    a.output.mkdir(exist_ok=False); assert not a.receipt.exists()
    rng = np.random.default_rng(202610041); exports = []; reordered = permutations = transforms = isolation = invalid = 0
    maximum_raw = maximum_projected = 0.; reviews = 0
    for exception, mode in itertools.product([False, True], ['signed', 'unsigned']):
        src, rows, x, y, ops, audit, old, parent = inputs(exception, mode)
        n = len(rows); labels = src['labels']; factor = src['factors'][audit['tree']].copy()
        operators = {'target_node': sparse.eye(n, format='csr'), 'background_node': ops['background_node']}
        if exception: operators['model_pair'] = ops['model_pair']
        operators['family_intercept'] = ops['family_intercept']
        cached = IndependentPositiveDiagonalBasisContext(labels, operators, factor).design(x)
        primary = PositiveDiagonalBasisContext(labels, operators, factor).design(x)
        patterns = [np.ones(n), np.full(n, 2.), np.where(np.arange(n) % 2, .5, 2.),
            np.linspace(.125, 3., n), 1. + np.arange(n) * 2.**-40, np.geomspace(1e-4, 1e4, n)]
        first = []
        for j, d in enumerate(patterns):
            actual = cached.audit(d)
            separate = PositiveDiagonalKernelProducts(labels, operators, d).tree(factor).project(x)
            component = primary.audit(d); reference = dense(operators, factor, d, x)
            for k in ['raw', 'projected', 'raw_error', 'projected_error', 'diagonal_core', 'diagonal_image_energy']:
                np.testing.assert_allclose(actual[k], separate[k], rtol=3e-9, atol=2e-8)
            for kind, field, expected in [('raw', 'raw_gram', reference[0]), ('projected', 'projected_gram', reference[1])]:
                np.testing.assert_allclose(actual[kind], expected, rtol=3e-9, atol=2e-8)
                np.testing.assert_allclose(actual[kind], component[field], rtol=3e-9, atol=2e-8)
                bound = actual[kind + '_error'] + np.asarray(component['raw_roundoff_envelope' if kind == 'raw' else 'projected_roundoff_envelope'])
                assert np.all(abs(actual[kind] - expected) <= bound + 64 * np.finfo(float).tiny)
                error = float(np.max(abs(actual[kind] - expected)))
                if kind == 'raw': maximum_raw = max(maximum_raw, error)
                else: maximum_projected = max(maximum_projected, error)
                status = independent_diagnostics(actual[kind], actual[kind + '_error'], cached.parent.bank.names)
                assert status['disposition'] == component['raw_diagnostics' if kind == 'raw' else 'reml_diagnostics']['disposition']
                if j in [0, 1]:
                    assert status['disposition'] == 'dependent_covariance_bases_require_review'; reviews += 1
            permutation = rng.permutation(n)
            changed = IndependentPositiveDiagonalBasisContext(labels[permutation], {k:z[permutation] for k,z in operators.items()}, factor[permutation]).design(x[permutation]).audit(d[permutation])
            for k in ['raw', 'projected']: np.testing.assert_allclose(changed[k], actual[k], rtol=3e-9, atol=2e-8)
            permutations += 1
            t = rng.normal(size=(x.shape[1], x.shape[1])) + 4 * np.eye(x.shape[1])
            changed = IndependentPositiveDiagonalBasisContext(labels, operators, factor).design(x @ t).audit(d)
            np.testing.assert_allclose(changed['projected'], actual['projected'], rtol=3e-9, atol=2e-8); transforms += 1
            first.append(actual)
            exports.append(dict(pair_exception=exception, loading_mode=mode, diagonal_case=j,
                actual={k:v.tolist() if isinstance(v, np.ndarray) else v for k,v in actual.items()}))
        for j in rng.permutation(6):
            again = cached.audit(patterns[j])
            for k in first[j]:
                if isinstance(first[j][k], np.ndarray): assert np.array_equal(again[k], first[j][k])
                else: assert again[k] == first[j][k]
            reordered += 1
        # Caller-owned factor, X, incidence, returned results and D cannot
        # silently change any already constructed nonresidual cache.
        factor[:] = np.nan; x[:] = np.nan
        for z in operators.values(): z.data[:] = 0
        returned = cached.audit(patterns[2]); returned['raw'][:] = np.nan
        d = patterns[2].copy(); result = cached.audit(d); d[:] = np.nan
        for k in ['raw', 'projected']: assert np.array_equal(cached.audit(patterns[2])[k], first[2][k])
        isolation += 5
        for d in [np.zeros(n), -np.ones(n), np.full(n, np.nan), np.full(n, np.inf), np.ones(n + 1)]:
            rejected(lambda: cached.audit(d)); invalid += 1
        # Fresh valid inputs for exact named uniform routes, including q4/q5.
        src, rows, x, y, ops, audit, old, parent = inputs(exception, mode)
        uniform_ops = {'background_node': ops['background_node']}
        if exception: uniform_ops['model_pair'] = ops['model_pair']
        uniform_ops['family_intercept'] = ops['family_intercept']
        f = src['factors'][audit['tree']]
        actual = IndependentPositiveDiagonalBasisContext(src['labels'], uniform_ops, f).design(x).audit(np.ones(n))
        reference = dense(uniform_ops, f, np.ones(n), x)
        np.testing.assert_allclose(actual['raw'], reference[0], rtol=3e-9, atol=2e-8)
        np.testing.assert_allclose(actual['projected'], reference[1], rtol=3e-9, atol=2e-8)
        exports.append(dict(pair_exception=exception, loading_mode=mode, diagonal_case='named_uniform_basis',
            actual={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in actual.items()}))
        zero = {**uniform_ops, 'explicit_zero': sparse.csr_matrix((n, 0))}
        result = IndependentPositiveDiagonalBasisContext(src['labels'], zero, f).design(x).audit(np.ones(n))
        assert independent_diagnostics(result['raw'], result['raw_error'], ['residual', *zero, 'species'])['disposition'] == 'unresolved_kernel_norm_requires_review'
        for bad in [x[:, :1] * 0, np.column_stack([x, x[:, 0]]), np.ones((n, n)), x[:-1], np.full_like(x, np.nan)]:
            rejected(lambda: IndependentPositiveDiagonalBasisContext(src['labels'], uniform_ops, f).design(bad)); invalid += 1
    artifact = a.output / 'independent_cached_diagonal_cases.json'
    with artifact.open('x') as f: json.dump(exports, f, sort_keys=True, allow_nan=False)
    modules = ['independent_positive_diagonal_basis_context', 'check_independent_positive_diagonal_basis_context',
        'independent_positive_diagonal_kernel_products', 'positive_diagonal_basis_context', 'covariance_basis_context',
        'covariance_basis_audit', 'reduced_covariance_basis', 'covariance_exact_folds_v2',
        'check_positive_diagonal_kernel_products', 'check_retained_shared_entity_candidates', 'run_positive_diagonal_kernel_software_stage']
    result = dict(status='passed_independent_cached_positive_diagonal_basis_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(), general_diagonal_cases=24, named_uniform_basis_cases=4,
        row_permutations=permutations, full_rank_design_transforms=transforms, reordered_calls=reordered,
        caller_input_and_result_isolation_checks=isolation, constant_diagonal_diagnostic_reviews=reviews,
        explicit_zero_norm_cases=4, invalid_diagonal_design_cases_rejected=invalid,
        maximum_dense_raw_error=maximum_raw, maximum_dense_projected_error=maximum_projected,
        unchanged_comparison_rtol=3e-9, unchanged_comparison_atol=2e-8,
        source_hashes={'scripts/' + m + '.py': sha('scripts/' + m + '.py') for m in modules},
        artifacts={str(artifact):sha(artifact)}, scientific_eligibility=False,
        full_real_data_weighted_qualification_complete=False,
        scope='Complete declared q4/q5/q6 mode/diagonal comparisons to separate latent calculations, component cache and explicit dense error contrasts. Fresh bounds, same strict comparisons, constant/near-uniform dependence and zero norm reviews retained. No saved uniform audit or measured real-data speedup, fit or biological acceptance.')
    with a.receipt.open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes', 'artifacts']}), flush=True)


if __name__ == '__main__': main()
