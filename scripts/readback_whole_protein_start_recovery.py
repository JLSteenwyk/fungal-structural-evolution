#!/usr/bin/env python3
"""Replay recovery lineage, all likelihoods, GLS summaries and finite differences."""
import argparse
import json
import subprocess
import time
from pathlib import Path
import numpy as np
import psutil
from cached_matched_ml import CachedMatchedML
from matched_ml_gradient import evaluate_gradient
from matched_mixed_covariance import MatchedCovariance
from readback_selected_refinement_candidates import independent_fit
from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **config['pins']}
    launch = json.loads(Path(config['launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            process = psutil.Process(launch['pid'])
            if process.create_time() != launch['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            assert process.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(line.split('=', 1) for line in subprocess.check_output(
        ['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    assert sha(launch['plan']) == launch['plan_sha256']
    plan = json.loads(Path(launch['plan']).read_text())
    receipt_path = Path(plan['output'])/'receipt.json'
    saved = json.loads(receipt_path.read_text())
    bindings.update(saved['source_hashes'])
    bindings[str(receipt_path)] = sha(receipt_path)
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    source = json.loads(Path(plan['source']).read_text())
    old = source['payload']
    assert saved['source'] == plan['source'] and saved['fit_input_id'] == source['fit_input_id']
    assert saved['tree'] == source['tree'] == Path(plan['factor']).stem
    assert sha(plan['input']) == source['input_sha256']
    assert Path(plan['input']).stem == source['fit_input_id']
    assert saved['original_candidates'] == old['candidates'] and len(old['candidates']) == 24
    proof = json.loads(Path(plan['readback']).read_text())
    assert proof['status'] == 'passed_serialized_frozen_whole_protein_refinement_readback'
    assert proof['source_hashes'][plan['source']] == sha(plan['source'])
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
    assert saved['optimizer'] == dict(method='SLSQP', ftol=1e-12, maxiter=1000, bounds=[0., upper])
    def independent(theta):
        theta = np.asarray(theta)
        assert theta.shape == (3,) and np.isfinite(theta).all() and (theta >= 0).all() and (theta <= upper).all()
        return independent_fit(MatchedCovariance(bg, family, factor, 1., *np.expm1(theta)), x, y)
    maximum_error = 0.
    def replay(row):
        nonlocal maximum_error
        theta = np.asarray(row['theta'])
        value = independent(theta)['negative_profiled_ml']
        np.testing.assert_allclose(value, row['objective'], rtol=1e-9, atol=1e-7)
        maximum_error = max(maximum_error, abs(value-row['objective']))
        gradient = evaluate_gradient(cache, np.expm1(theta))['log1p_ratio_gradient']
        projected = np.where(theta <= 1e-7, np.minimum(gradient, 0),
                             np.where(theta >= upper-1e-7, np.maximum(gradient, 0), gradient))
        np.testing.assert_allclose(row['gradient'], gradient, rtol=1e-9, atol=1e-9)
        np.testing.assert_allclose(row['projected_gradient'], projected, rtol=1e-9, atol=1e-9)
        assert row['maximum_absolute_projected_gradient'] == float(max(abs(projected)))
    for c in old['candidates']:
        np.testing.assert_allclose(independent(c['theta'])['negative_profiled_ml'], c['objective'], rtol=1e-9, atol=1e-7)
    full = [(i, c) for i, c in enumerate(old['candidates']) if all(c['active']) and c['start'] != 'original_unoptimized_reference']
    expected = [(i, c['start'], c['theta']) for i, c in full] + [(None, '1.0_original_initialization', [float(np.log(2))]*3)]
    runs = saved['recoveries']
    assert len(expected) == len(runs) == 5
    for run, (index, start, theta) in zip(runs, expected):
        assert run['source_candidate_index'] == index and run['start'] == start
        np.testing.assert_array_equal(run['initial']['theta'], theta)
        assert isinstance(run['success'], bool) and run['message'] and run['evaluations'] > 0 and run['iterations'] > 0
        replay(run['initial'])
        replay(run['endpoint'])
    choices = [(c['objective'], c['success'], 'original_candidate', i, c['theta']) for i, c in enumerate(old['candidates'])]
    choices += [(r['endpoint']['objective'], r['success'], 'recovery_endpoint', i, r['endpoint']['theta']) for i, r in enumerate(runs)]
    best = min(choices, key=lambda v: (v[0], not v[1]))
    selection = saved['selection']
    assert selection['kind'] == best[2] and selection['index'] == best[3]
    np.testing.assert_array_equal(selection['theta'], best[4])
    replay(selection)
    fitted = independent(selection['theta'])
    for field, value in fitted.items():
        np.testing.assert_allclose(saved['fitted'][field], value, rtol=1e-8, atol=1e-8)
    assert saved['fitted']['variance_profile_denominator'] == len(y)
    assert saved['fitted']['residual_degrees_of_freedom'] == len(y)-x.shape[1]
    conversion = np.r_[1., 1/scales]
    np.testing.assert_allclose(saved['raw_unit_beta'], fitted['beta']*conversion, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(saved['raw_unit_conditional_beta_covariance'], fitted['conditional_beta_covariance']*np.outer(conversion, conversion), rtol=1e-9, atol=1e-9)
    values = [r['endpoint']['objective'] for r in runs]
    checks = dict(all_recovery_optimizers_successful=all(r['success'] for r in runs),
                  all_recovery_projected_gradients_pass=all(r['endpoint']['maximum_absolute_projected_gradient'] <= 1e-3 for r in runs),
                  all_full_face_recovered_starts_agree=max(values)-min(values) <= 1e-5,
                  selected_optimizer_success=best[1], selected_projected_gradient_pass=selection['maximum_absolute_projected_gradient'] <= 1e-3,
                  no_recovery_upper_bound_contact=all(max(r['endpoint']['theta']) < upper-1e-6 for r in runs),
                  selected_no_upper_bound_contact=max(selection['theta']) < upper-1e-6,
                  original_objective_not_worsened=best[0] <= old['negative_profiled_ml'],
                  recovered_starts_agree_with_selected=max(abs(v-selection['objective']) for v in values) <= 1e-5)
    assert checks == saved['checks']
    assert saved['status'] == ('full_face_start_recovery_passed_pending_readback' if all(checks.values()) else 'full_face_start_recovery_requires_review')
    # Differentiate the separate normal-equation likelihood, including one-sided
    # derivatives at a boundary. Compare two step sizes at every bad initial
    # endpoint and at the selected solution; do not differentiate the analytic score.
    points = [('selected', selection)] + [(r['start'], r['initial']) for r in runs if r['initial']['maximum_absolute_projected_gradient'] > 1e-3]
    finite = []
    for label, row in points:
        theta = np.array(row['theta'])
        for step in [1e-5, 1e-6]:
            values_fd = []
            for axis in range(3):
                direction = np.eye(3)[axis]*step
                if theta[axis] < step:
                    value = (-3*independent(theta)['negative_profiled_ml'] + 4*independent(theta+direction)['negative_profiled_ml'] - independent(theta+2*direction)['negative_profiled_ml'])/(2*step)
                elif theta[axis] > upper-step:
                    value = (3*independent(theta)['negative_profiled_ml'] - 4*independent(theta-direction)['negative_profiled_ml'] + independent(theta-2*direction)['negative_profiled_ml'])/(2*step)
                else:
                    value = (independent(theta+direction)['negative_profiled_ml']-independent(theta-direction)['negative_profiled_ml'])/(2*step)
                values_fd.append(value)
            np.testing.assert_allclose(values_fd, row['gradient'], rtol=1e-4, atol=1e-3)
            finite.append(dict(point=label, step=step, gradient=values_fd, maximum_absolute_error=float(max(abs(np.array(values_fd)-row['gradient'])))))
    verify()
    result = dict(status='passed_serialized_full_face_start_recovery_readback', numerical_recovery_passed=all(checks.values()),
                  fit_input_id=saved['fit_input_id'], tree=saved['tree'], checks=checks, original_candidates_replayed=24,
                  recovered_starts=5, maximum_objective_error=maximum_error, finite_difference_checks=finite,
                  source_hashes=bindings, producer_terminal_state=state,
                  scope='Every retained candidate, recovered-start lineage, endpoint, selection and GLS transform replayed. Analytic gradients independently compared with direct-likelihood finite differences. Covariance operator shared. Numerical recovery is not global optimality, full-grid acceptance, production integration or calibrated scientific inference.')
    with Path(config['output']).open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}), flush=True)


if __name__ == '__main__':
    main()
