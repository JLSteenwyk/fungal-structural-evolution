#!/usr/bin/env python3
"""Confirm diagnostic Newton proposals with direct finite differences."""
import json
from pathlib import Path
import numpy as np
from cached_matched_ml import profiled_ml
from matched_mixed_covariance import MatchedCovariance
from screen_duplication_alignment_reuse import sha


def main():
    rp = Path('results/model_validation/refinement-curvature-diagnostic-20260928-v1/receipt.json')
    r = json.loads(rp.read_text())
    bindings = {str(rp): sha(rp), **r['source_hashes']}
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    plan = json.loads(Path('metadata/refinement_gradient_diagnostic_plan_20260928.json').read_text())
    with np.load(plan['inputs']) as a:
        values, bg, family, indices = a['matrix'], a['background'], a['family'], a['pattern_rows']
    with np.load(plan['factor']) as a:
        factor = a['factor'][indices]
    active = values[:, 1:].max(axis=0) - values[:, 1:].min(axis=0) > 1e-12
    covariates = values[:, 1:][:, active]
    x = np.column_stack([np.ones(len(values)), covariates / covariates.std(axis=0)])
    y = values[:, 0]
    def objective(theta):
        assert np.isfinite(theta).all() and (theta >= 0).all() and (theta <= np.log1p(10000.)).all()
        return profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(theta)), x, y)['negative_profiled_ml']
    baseline = objective(np.asarray(r['base_theta']))
    np.testing.assert_allclose(baseline, r['base_objective'], rtol=1e-12, atol=1e-10)
    proposals = [t for t in r['trials'] if t['step_fraction'] == 1.]
    assert len(proposals) == 2
    checked = []
    for trial in proposals:
        theta = np.asarray(trial['theta'])
        value = objective(theta)
        assert value <= baseline
        np.testing.assert_allclose(value, trial['objective'], rtol=1e-12, atol=1e-10)
        for step in [1e-6, 1e-7]:
            gradient = []
            for axis in range(3):
                offset = np.zeros(3)
                offset[axis] = step
                gradient.append((objective(theta + offset) - objective(theta - offset)) / (2 * step))
            # The original acceptance criterion is unchanged; record each
            # step estimate instead of substituting a favorable step size.
            assert max(map(abs, gradient)) <= 1e-3
            checked.append(dict(hessian_step=trial['hessian_step'], finite_difference_step=step,
                                candidate_theta=theta.tolist(), objective=value,
                                finite_difference_gradient=gradient,
                                maximum_absolute_gradient=max(map(abs, gradient))))
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    result = dict(status='passed_direct_finite_difference_check_of_two_curvature_proposals',
                  source_hashes=bindings, checker_sha256=sha(__file__), proposals=2,
                  direct_objective_evaluations=27, checks=checked,
                  scope='Both fixed proposals checked at two finite-difference steps in all three interior directions. Existing1e-3 criterion retained. No production replacement, all-start acceptance, global-optimum or scientific inference claim.')
    with Path('metadata/refinement_curvature_proposals_checked_20260928.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
