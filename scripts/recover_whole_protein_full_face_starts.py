#!/usr/bin/env python3
"""Retain audited starts and recover their endpoints with bounded SLSQP."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import minimize
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_ml_gradient import evaluate_gradient
from matched_mixed_covariance import MatchedCovariance
from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    proof = json.loads(Path(plan['readback']).read_text())
    assert proof['status'] == 'passed_serialized_frozen_whole_protein_refinement_readback'
    bindings.update(proof['source_hashes'])
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    source = json.loads(Path(plan['source']).read_text())
    old = source['payload']
    assert old['status'] == 'ml_refinement_requires_review'
    assert old['checks']['all_full_face_starts_agree'] is False
    with np.load(plan['input'], allow_pickle=False) as a:
        matrix, bg, family, rows = a['matrix'], a['background'], a['family'], a['pattern_rows']
    with np.load(plan['factor'], allow_pickle=False) as a:
        factor = a['factor'][rows]
    scales = np.std(matrix[:, 2:], axis=0)
    np.testing.assert_array_equal(scales, source['covariate_scales'])
    assert np.all(matrix[:, 1] == 1)
    x = np.column_stack([matrix[:, 1], matrix[:, 2:] / scales])
    y = matrix[:, 0]
    cache = CachedMatchedML(bg, family, factor, x, y)
    upper = float(np.log1p(old['maximum_ratio']))
    def objective(theta):
        ratios = np.expm1(theta)
        value = profiled_ml(MatchedCovariance(bg, family, factor, 1., *ratios), x, y)['negative_profiled_ml']
        gradient = evaluate_gradient(cache, ratios)
        np.testing.assert_allclose(value, gradient['negative_profiled_ml'], rtol=1e-9, atol=1e-7)
        return value, gradient['log1p_ratio_gradient']
    def describe(theta):
        value, gradient = objective(theta)
        projected = np.where(theta <= 1e-7, np.minimum(gradient, 0),
                             np.where(theta >= upper - 1e-7, np.maximum(gradient, 0), gradient))
        return dict(theta=theta.tolist(), objective=float(value), gradient=gradient.tolist(),
                    projected_gradient=projected.tolist(), maximum_absolute_projected_gradient=float(max(abs(projected))))
    starts = [(i, c) for i, c in enumerate(old['candidates'])
              if all(c['active']) and c['start'] != 'original_unoptimized_reference']
    assert len(starts) == 4 and {c['start'] for _, c in starts} == {'0.05', '1.0', '20.0', 'original_parameters'}
    output = Path(plan['output'])
    output.mkdir(parents=True, exist_ok=False)
    recoveries = []
    began = time.perf_counter()
    for i, candidate in starts + [(None, dict(start='1.0_original_initialization', theta=[float(np.log(2))]*3))]:
        initial = np.asarray(candidate['theta'])
        before = describe(initial)
        fit = minimize(objective, initial, method='SLSQP', jac=True,
                       bounds=[(0., upper)]*3, options=dict(ftol=1e-12, maxiter=1000))
        after = describe(fit.x)
        row = dict(source_candidate_index=i, start=candidate['start'], initial=before, endpoint=after,
                   success=bool(fit.success), message=str(fit.message), iterations=int(fit.nit), evaluations=int(fit.nfev))
        recoveries.append(row)
        print(json.dumps(row), flush=True)
    # All original candidates remain eligible; a lower failed candidate cannot be hidden.
    choices = [(c['objective'], c['success'], 'original_candidate', i, np.array(c['theta']))
               for i, c in enumerate(old['candidates'])]
    choices += [(r['endpoint']['objective'], r['success'], 'recovery_endpoint', i, np.array(r['endpoint']['theta']))
                for i, r in enumerate(recoveries)]
    selected = min(choices, key=lambda v: (v[0], not v[1]))
    selected_summary = describe(selected[4])
    fitted = profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(selected[4])), x, y)
    conversion = np.r_[1., 1/scales]
    values = [r['endpoint']['objective'] for r in recoveries]
    checks = dict(all_recovery_optimizers_successful=all(r['success'] for r in recoveries),
                  all_recovery_projected_gradients_pass=all(r['endpoint']['maximum_absolute_projected_gradient'] <= 1e-3 for r in recoveries),
                  all_full_face_recovered_starts_agree=bool(np.ptp(values) <= 1e-5),
                  selected_optimizer_success=bool(selected[1]),
                  selected_projected_gradient_pass=selected_summary['maximum_absolute_projected_gradient'] <= 1e-3,
                  no_recovery_upper_bound_contact=all(max(r['endpoint']['theta']) < upper-1e-6 for r in recoveries),
                  selected_no_upper_bound_contact=bool(max(selected[4]) < upper-1e-6),
                  original_objective_not_worsened=bool(selected[0] <= old['negative_profiled_ml']),
                  recovered_starts_agree_with_selected=bool(max(abs(v-selected_summary['objective']) for v in values) <= 1e-5))
    verify()
    result = dict(status='full_face_start_recovery_passed_pending_readback' if all(checks.values()) else 'full_face_start_recovery_requires_review',
                  fit_input_id=source['fit_input_id'], tree=source['tree'], source=plan['source'],
                  original_candidates=old['candidates'], recoveries=recoveries,
                  selection=dict(kind=selected[2], index=selected[3], **selected_summary), checks=checks,
                  fitted={k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in fitted.items()},
                  raw_unit_beta=(fitted['beta']*conversion).tolist(),
                  raw_unit_conditional_beta_covariance=(fitted['conditional_beta_covariance']*np.outer(conversion, conversion)).tolist(),
                  source_hashes=bindings, optimizer=dict(method='SLSQP', ftol=1e-12, maxiter=1000, bounds=[0., upper]),
                  versions=dict(numpy=np.__version__, scipy=scipy.__version__), wall_seconds=time.perf_counter()-began,
                  scope='All four original full-face endpoints recovered separately, plus the original ratio-one initialization. Original 24 candidates retained and considered for selection. Numerical recovery only; source flags and production fits preserved. No global-optimum proof, production overlay or calibrated inference.')
    (output/'receipt.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'], checks=checks)), flush=True)


if __name__ == '__main__':
    main()
