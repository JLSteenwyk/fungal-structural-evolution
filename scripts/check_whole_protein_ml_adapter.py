#!/usr/bin/env python3
"""Dense Gaussian validation for five whole-protein designs and two responses."""
import json
from pathlib import Path
import numpy as np
from fit_whole_protein_ml import fit_input
from check_matched_ml import dense
from screen_duplication_alignment_reuse import sha


def main():
    rng = np.random.default_rng(20260929)
    n = 48
    bg = np.repeat(np.arange(12), 4)
    fam = bg // 3
    factor = rng.normal(size=(n, 3)) * .2
    kernels = [(bg[:, None] == bg).astype(float), (fam[:, None] == fam).astype(float), factor @ factor.T]
    t, b = rng.uniform(.1, 2., size=(2, n))
    nuisance = rng.normal(size=(n, 3)) * [1., 10., .1]
    variants = dict(linear=(t-b)[:, None], quadratic=np.column_stack([t-b, t*t-b*b]),
                    cubic=np.column_stack([t-b, t*t-b*b, t**3-b**3]),
                    positive_log=(np.log(t)-np.log(b))[:, None], identity=rng.normal(size=(n, 1)))
    responses = dict(rmsd=rng.normal(size=n), tm=rng.normal(size=n)*.1)
    cases = []
    max_error = 0.
    for name, predictor in variants.items():
        design = np.column_stack([np.ones(n), predictor, nuisance])
        columns = ['response', 'intercept'] + ['x'+str(i) for i in range(design.shape[1]-1)]
        for outcome, y in responses.items():
            fitted = fit_input(bg, fam, factor, np.column_stack([y, design]), columns)
            for c in fitted['candidates']:
                reference = dense(design, y, kernels, np.expm1(c['theta']))
                error = abs(reference['negative_profiled_ml'] - c['objective'])
                max_error = max(max_error, error)
                np.testing.assert_allclose(reference['negative_profiled_ml'], c['objective'], rtol=1e-8, atol=1e-7)
            reference = dense(design, y, kernels, fitted['ratios'])
            np.testing.assert_allclose(reference['beta'], fitted['raw_unit_beta'], rtol=1e-7, atol=1e-8)
            np.testing.assert_allclose(reference['conditional_beta_covariance'], fitted['raw_unit_conditional_beta_covariance'], rtol=1e-7, atol=1e-8)
            assert fitted['variance_profile_denominator'] == n
            assert fitted['residual_degrees_of_freedom'] == n-design.shape[1]
            assert len(fitted['candidates']) == 22
            cases.append(dict(variant=name, outcome=outcome, status=fitted['status'], candidates=22))
            print('Checked whole-protein adapter', name, outcome, flush=True)
    invalid = np.column_stack([responses['rmsd'], np.ones(n), np.zeros(n)])
    try:
        fit_input(bg, fam, factor, invalid, ['response', 'intercept', 'constant'])
    except ValueError:
        pass
    else:
        raise AssertionError('Inactive covariate accepted')
    result = dict(status='passed_whole_protein_ml_adapter_dense_fixtures', cases=cases,
                  candidate_likelihoods_checked=220, maximum_objective_error=max_error,
                  source_hashes={p:sha(p) for p in ['scripts/fit_whole_protein_ml.py', 'scripts/fit_matched_ml.py', 'scripts/check_matched_ml.py', 'scripts/check_whole_protein_ml_adapter.py']},
                  scope='Synthetic implementation validation for all five design forms and both response scales. Dense references check raw coefficients, conditional covariance and every candidate likelihood. Not empirical calibration, full-data fitting or global-optimum proof.')
    with Path('metadata/whole_protein_ml_adapter_fixtures_20260929.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
