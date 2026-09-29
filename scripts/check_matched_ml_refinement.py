#!/usr/bin/env python3
"""Validate ML refinement candidates against a dense covariance calculation."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from threadpoolctl import threadpool_limits
from refine_matched_ml_analytic import refine
from screen_duplication_alignment_reuse import sha


def main():
    threadpool_limits(1)
    rng = np.random.default_rng(20260928)
    n = 48
    background = np.repeat(np.arange(12), 4)
    family = background // 3
    factor = rng.normal(size=(n, 4)) * .2
    u, v = rng.uniform(.05, .95, (2, n))
    y = .3 + .5 * (u - v) + rng.normal(size=n) * .6
    checked = 0
    maximum_error = 0.
    statuses = {}
    for degree in [1, 2, 3]:
        design = np.column_stack([np.ones(n)] + [u ** k - v ** k for k in range(1, degree + 1)])
        for zero_species in [False, True]:
            f = np.zeros_like(factor) if zero_species else factor
            result = refine(background, family, f, design, y, np.log1p([.2, .5, .1]))
            assert result['likelihood'] == 'ordinary_gaussian_ml'
            assert result['variance_profile_denominator'] == n
            assert len(result['candidates']) == 24
            assert result['checks']['original_objective_not_worsened']
            assert result['candidates'][0]['start'] == 'original_unoptimized_reference'
            assert not result['candidates'][0]['success']
            assert sum(c['start'] == 'original_parameters' for c in result['candidates']) == 1
            for candidate in result['candidates']:
                b, g, s = np.expm1(candidate['theta'])
                covariance = np.eye(n) + b * (background[:, None] == background) + g * (family[:, None] == family) + s * (f @ f.T)
                chol = cho_factor(covariance)
                rx = cho_solve(chol, design)
                beta = np.linalg.solve(design.T @ rx, design.T @ cho_solve(chol, y))
                residual = y - design @ beta
                q = residual @ cho_solve(chol, residual)
                objective = .5 * (2 * np.log(np.diag(chol[0])).sum() + n * (1 + np.log(2 * np.pi * q / n)))
                np.testing.assert_allclose(objective, candidate['objective'], rtol=1e-9, atol=1e-7)
                maximum_error = max(maximum_error, abs(objective - candidate['objective']))
                checked += 1
            best = min(result['candidates'], key=lambda c: (c['objective'], not c['success']))
            assert result['checks']['best_optimizer_success'] == best['success']
            assert result['log1p_ratios'] == best['theta']
            statuses[result['status']] = statuses.get(result['status'], 0) + 1
    for initial in [[-1, 0, 0], [0, float('nan'), 0], [0, 0], [0, 0, 100]]:
        try:
            refine(background, family, factor, design, y, initial)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid original parameter vector accepted')
    result = dict(status='passed_dense_ml_refinement_fixtures', fits=6, dense_candidates_checked=checked,
                  invalid_initial_vectors_rejected=4, maximum_dense_objective_error=maximum_error,
                  status_counts=statuses, source_hashes={p: sha(p) for p in [__file__, 'scripts/refine_matched_ml_analytic.py', 'scripts/cached_matched_ml.py', 'scripts/matched_ml_gradient.py']},
                  scope='Synthetic numerical validation across linear/quadratic/cubic designs and zero/nonzero species factors; no production flags resolved and no global-optimum or calibration claim.')
    with Path('metadata/matched_ml_refinement_checks_20260928.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
