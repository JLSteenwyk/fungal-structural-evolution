#!/usr/bin/env python3
"""Evaluate bounded Newton proposals for one unresolved interior-gradient case."""
import json
from pathlib import Path
import numpy as np
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from screen_duplication_alignment_reuse import sha


def main():
    pp = Path('metadata/refinement_curvature_diagnostic_plan_20260928.json')
    plan = json.loads(pp.read_text())
    source = json.loads(Path(plan['source_plan']).read_text())
    bindings = {str(pp): sha(pp), **plan['pins'], **source['pins']}
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    fit = json.loads(Path(source['fit']).read_text())['payload']
    with np.load(source['inputs']) as a:
        matrix, bg, family, indices = a['matrix'], a['background'], a['family'], a['pattern_rows']
    with np.load(source['factor']) as a:
        factor = a['factor'][indices]
    x = np.column_stack([np.ones(len(matrix)), matrix[:, 1:][:, fit['active_covariates']] / fit['covariate_scales']])
    y = matrix[:, 0]
    theta = np.asarray(fit['log1p_ratios'])
    upper = np.log1p(fit['maximum_ratio'])
    cache = CachedMatchedML(bg, family, factor, x, y)
    def gradient(t):
        return evaluate_gradient(cache, np.expm1(t))['log1p_ratio_gradient']
    def direct(t):
        assert (t >= 0).all() and (t <= upper).all()
        return profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(t)), x, y)['negative_profiled_ml']
    base, g = direct(theta), gradient(theta)
    trials = []
    for step in [1e-5, 1e-6]:
        columns = []
        for axis in range(3):
            offset = np.zeros(3)
            offset[axis] = step
            assert (theta - offset >= 0).all() and (theta + offset <= upper).all()
            columns.append((gradient(theta + offset) - gradient(theta - offset)) / (2 * step))
        raw = np.column_stack(columns)
        hessian = (raw + raw.T) / 2
        eigenvalues = np.linalg.eigvalsh(hessian)
        assert (eigenvalues > 0).all()
        delta = np.linalg.solve(hessian, -g)
        for fraction in [1., .5, .25, .125]:
            candidate = theta + fraction * delta
            objective, derivatives = direct(candidate), gradient(candidate)
            trials.append(dict(hessian_step=step, step_fraction=fraction,
                               hessian_eigenvalues=eigenvalues.tolist(), hessian_asymmetry=float(abs(raw - raw.T).max()),
                               delta=delta.tolist(), theta=candidate.tolist(), objective=objective,
                               objective_change=objective - base, gradient=derivatives.tolist(),
                               maximum_absolute_gradient=float(abs(derivatives).max()),
                               direct_objective_not_worsened=bool(objective <= base),
                               interior_gradient_below_existing_threshold=bool(abs(derivatives).max() <= 1e-3)))
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    result = dict(status='complete_single_case_curvature_diagnostic', source_hashes=bindings,
                  script_sha256=sha(__file__), base_theta=theta.tolist(), base_objective=base,
                  base_gradient=g.tolist(), trials=trials,
                  scope='One interior case; two finite-difference Hessians and eight bounded Newton proposals. Every direct objective and analytic gradient retained. No production update, global-optimum claim or replacement of all-face/start checks; proposed points need separate finite-difference confirmation.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
