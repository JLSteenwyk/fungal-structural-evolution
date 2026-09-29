#!/usr/bin/env python3
"""Replay all frozen-snapshot refinement candidates and output classifications."""
import csv
import itertools
import json
import subprocess
import time
from collections import Counter
from pathlib import Path
import numpy as np
import psutil
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_ml_gradient import evaluate_gradient
from matched_mixed_covariance import MatchedCovariance
from screen_duplication_alignment_reuse import sha


def main():
    lp = Path('metadata/flagged_polynomial_snapshot_refinement_launch_20260928.json')
    launch = json.loads(lp.read_text())
    bindings = {str(lp): sha(lp), launch['plan']: launch['plan_sha256']}
    while psutil.pid_exists(launch['pid']):
        try:
            p = psutil.Process(launch['pid'])
            if abs(p.create_time() - launch['created']) > .01 or p.status() == psutil.STATUS_ZOMBIE:
                break
            assert p.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(x.split('=', 1) for x in subprocess.check_output(['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    plan = json.loads(Path(launch['plan']).read_text())
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    rh = sha(root / 'receipt.json')
    assert receipt['status'] == 'complete_flagged_polynomial_snapshot_refinements_pending_independent_readback'
    bindings.update(receipt['source_hashes'])
    bindings.update({str(root / n): h for n, h in receipt['artifacts'].items()})
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    inputs = Path(plan['inputs'])
    recipes = {r['recipe']['fit_input_id']: r for r in json.loads((inputs / 'input_manifest.json').read_text())}
    expected = {(r['fit_input_id'], r['tree']): r for r in csv.DictReader((inputs / 'flagged_tree_fits.tsv').open(), delimiter='\t')}
    factors = {}
    for path in Path(plan['factors']).glob('*.npz'):
        with np.load(path) as arrays:
            factors[path.stem] = arrays['factor']
    manifest = json.loads((root / 'fit_manifest.json').read_text())
    seen, counts = set(), Counter()
    candidate_count, maximum_error = 0, 0.
    for item in manifest:
        key = item['fit_input_id'], item['tree']
        assert key in expected and key not in seen
        seen.add(key)
        path = root / item['path']
        assert sha(path) == item['sha256']
        saved = json.loads(path.read_text())
        original = expected[key]
        assert saved['source_fit_path'] == original['source_path'] and saved['source_fit_sha256'] == original['source_sha256']
        assert (saved['fit_input_id'], saved['tree']) == key and saved['polynomial_degree'] == int(original['polynomial_degree'])
        old = json.loads(Path(original['source_path']).read_text())['payload']
        recipe = recipes[key[0]]
        assert saved['input_sha256'] == recipe['sha256']
        with np.load(inputs / recipe['path']) as arrays:
            matrix = arrays['matrix']
            bg, family, factor = arrays['background'], arrays['family'], factors[key[1]][arrays['pattern_rows']]
        active = matrix[:, 1:].max(axis=0) - matrix[:, 1:].min(axis=0) > 1e-12
        scales = matrix[:, 1:][:, active].std(axis=0)
        x = np.column_stack([np.ones(len(matrix)), matrix[:, 1:][:, active] / scales])
        y = matrix[:, 0]
        p = saved['payload']
        np.testing.assert_array_equal(p['active_covariates'], active)
        np.testing.assert_array_equal(p['covariate_scales'], scales)
        candidates = p['candidates']
        grid = [([True] * 3, 'original_unoptimized_reference')]
        for face in itertools.product([False, True], repeat=3):
            for start in (['0.05', '1.0', '20.0'] if any(face) else ['exact_zero']):
                grid.append((list(face), start))
            if all(face):
                grid.append((list(face), 'original_parameters'))
        assert [(c['active'], c['start']) for c in candidates] == grid and len(candidates) == 24
        np.testing.assert_array_equal(candidates[0]['theta'], old['log1p_ratios'])
        assert not candidates[0]['success'] and p['maximum_ratio'] == old['maximum_ratio']
        upper = np.log1p(p['maximum_ratio'])
        for c in candidates:
            theta = np.asarray(c['theta'])
            assert theta.shape == (3,) and np.isfinite(theta).all() and (theta >= 0).all() and (theta <= upper).all()
            assert (theta[~np.asarray(c['active'])] == 0).all()
            direct = profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(theta)), x, y)
            error = abs(direct['negative_profiled_ml'] - c['objective'])
            np.testing.assert_allclose(direct['negative_profiled_ml'], c['objective'], rtol=1e-9, atol=1e-7)
            maximum_error = max(maximum_error, error)
            candidate_count += 1
        best = min(candidates, key=lambda c: (c['objective'], not c['success']))
        theta = np.asarray(best['theta'])
        np.testing.assert_array_equal(p['log1p_ratios'], theta)
        np.testing.assert_allclose(p['ratios'], np.expm1(theta), rtol=1e-12, atol=1e-12)
        grad = evaluate_gradient(CachedMatchedML(bg, family, factor, x, y), np.expm1(theta))['log1p_ratio_gradient']
        np.testing.assert_allclose(p['analytic_gradient'], grad, rtol=1e-9, atol=1e-9)
        projected = np.where(theta <= 1e-7, np.minimum(grad, 0), np.where(theta >= upper - 1e-7, np.maximum(grad, 0), grad))
        np.testing.assert_allclose(p['projected_gradient'], projected, rtol=1e-9, atol=1e-9)
        full = [c['objective'] for c in candidates if all(c['active']) and c['start'] != 'original_unoptimized_reference']
        checks = dict(best_optimizer_success=best['success'], projected_gradient_pass=bool(max(abs(projected)) <= 1e-3), upper_bound_contact=bool(any(theta >= upper - 1e-6)), all_full_face_starts_agree=bool(max(full) - min(full) <= 1e-5), original_objective_not_worsened=bool(best['objective'] <= candidates[0]['objective']))
        assert p['checks'] == checks
        passed = all(v for k, v in checks.items() if k != 'upper_bound_contact') and not checks['upper_bound_contact']
        assert p['status'] == item['status'] == ('ml_refinement_passed_numerical_checks' if passed else 'ml_refinement_requires_review')
        fitted = profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(theta)), x, y)
        for name in ['negative_profiled_ml', 'beta', 'profiled_scale', 'residual_quadratic', 'conditional_beta_covariance']:
            np.testing.assert_allclose(p[name], fitted[name], rtol=1e-8, atol=1e-8)
        conversion = np.r_[1., 1. / scales]
        np.testing.assert_allclose(p['raw_unit_beta'], fitted['beta'] * conversion, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(p['raw_unit_conditional_beta_covariance'], fitted['conditional_beta_covariance'] * conversion[:, None] * conversion[None, :], rtol=1e-10, atol=1e-10)
        assert p['variance_profile_denominator'] == len(y) and p['likelihood'] == 'ordinary_gaussian_ml'
        assert p['original_objective'] == candidates[0]['objective'] and p['objective_improvement'] == candidates[0]['objective'] - best['objective']
        counts[p['status']] += 1
        if len(seen) % 20 == 0:
            print('Replayed refinement outputs', len(seen), '/601', flush=True)
    assert seen == set(expected) and len(seen) == receipt['fits'] == 601
    assert dict(counts) == receipt['counts'] and candidate_count == 14424
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    assert sha(root / 'receipt.json') == rh
    result = dict(status='passed_full_flagged_polynomial_refinement_replay', source_receipt_sha256=rh, checker_sha256=sha(__file__), producer_terminal_state=state, fits=len(seen), candidate_likelihoods_replayed=candidate_count, maximum_objective_error=maximum_error, counts=dict(counts), scope='Separate-process replay of every candidate and output classification using the validated direct likelihood and analytic-gradient libraries. Shared mathematical libraries remain; this is not a second mathematical implementation, global-optimum proof, full-grid acceptance or inferential calibration.')
    with Path('metadata/flagged_polynomial_refinement_completed_readback_20260928.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
