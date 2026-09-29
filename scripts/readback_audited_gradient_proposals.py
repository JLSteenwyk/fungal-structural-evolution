#!/usr/bin/env python3
"""Replay every disposition and local proposal without the proposal helper."""
import json
import subprocess
import time
from collections import Counter
from pathlib import Path
import numpy as np
import psutil
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from screen_duplication_alignment_reuse import sha


def check_proposal(result, objective, gradient, theta, upper):
    theta = np.asarray(theta)
    if np.any(theta <= 1e-5) or np.any(theta >= upper - 1e-5):
        assert result == dict(status='not_applicable_near_boundary', proposals=[])
        return
    base, g = objective(theta), gradient(theta)
    np.testing.assert_allclose(result['original_objective'], base, rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(result['original_gradient'], g, rtol=1e-9, atol=1e-9)
    if np.linalg.norm(g, ord=np.inf) <= .001:
        assert result['status'] == 'no_gradient_correction_needed' and result['proposals'] == []
        return
    np.testing.assert_array_equal(result['original_theta'], theta)
    assert result['upper'] == upper and len(result['proposals']) == 2
    passed = []
    for step, trial in zip((1e-5, 1e-6), result['proposals']):
        assert trial['hessian_step'] == step
        offsets = np.eye(3) * step
        raw = np.array([(gradient(theta + d) - gradient(theta - d)) / (2 * step) for d in offsets]).T
        hessian = .5 * (raw + raw.T)
        eigenvalues = np.linalg.eigvalsh(hessian)
        np.testing.assert_allclose(trial['hessian'], hessian, rtol=1e-9, atol=1e-8)
        np.testing.assert_allclose(trial['eigenvalues'], eigenvalues, rtol=1e-9, atol=1e-8)
        if not np.isfinite(eigenvalues).all() or np.any(eigenvalues <= 0):
            assert not trial['passed'] and trial['reason'] == 'nonpositive_or_nonfinite_curvature'
            passed.append(False)
            continue
        point = theta - np.linalg.solve(hessian, g)
        np.testing.assert_allclose(trial['theta'], point, rtol=1e-12, atol=1e-12)
        if not np.isfinite(point).all() or np.any(point <= 1e-6) or np.any(point >= upper - 1e-6):
            assert not trial['passed'] and trial['reason'] == 'proposal_not_interior_for_direct_difference_checks'
            passed.append(False)
            continue
        value, derivative = objective(point), gradient(point)
        np.testing.assert_allclose(trial['objective'], value, rtol=1e-12, atol=1e-9)
        np.testing.assert_allclose(trial['gradient'], derivative, rtol=1e-9, atol=1e-9)
        assert len(trial['finite_differences']) == 2
        norms = []
        for h, saved in zip((1e-6, 1e-7), trial['finite_differences']):
            assert saved['step'] == h
            fd = np.array([(objective(point + d) - objective(point - d)) / (2 * h) for d in np.eye(3) * h])
            np.testing.assert_allclose(saved['gradient'], fd, rtol=1e-9, atol=1e-9)
            norms.append(np.linalg.norm(fd, ord=np.inf))
        checks = dict(objective_not_worsened=bool(value <= base), analytic_gradient_pass=bool(np.linalg.norm(derivative, ord=np.inf) <= .001), both_direct_gradient_checks_pass=bool(max(norms) <= .001))
        assert trial['checks'] == checks and trial['passed'] == all(checks.values())
        passed.append(all(checks.values()))
    assert result['status'] == ('two_local_proposals_passed' if all(passed) else 'local_polish_requires_review')


def main():
    pp = Path('metadata/audited_gradient_proposals_readback_plan_20260929.json')
    plan = json.loads(pp.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    launch = json.loads(Path(plan['launch']).read_text())
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
    producer = json.loads(Path(launch['plan']).read_text())
    root = Path(producer['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_audited_gradient_proposals_pending_independent_readback'
    bindings = {str(pp): sha(pp), **plan['pins'], **receipt['source_hashes'], str(root / 'receipt.json'): sha(root / 'receipt.json')}
    bindings.update({str(root / name): digest for name, digest in receipt['artifacts'].items()})
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    source = Path(producer['refinements'])
    manifest = json.loads((source / 'fit_manifest.json').read_text())
    expected = {(r['fit_input_id'], r['tree']): r for r in manifest}
    assert len(expected) == len(manifest) == 601
    eligible = dict(best_optimizer_success=True, projected_gradient_pass=False, upper_bound_contact=False, all_full_face_starts_agree=True, original_objective_not_worsened=True)
    seen, counts, proposals = set(), Counter(), 0
    for line in (root / 'dispositions.jsonl').read_text().splitlines():
        row = json.loads(line)
        key = row['fit_input_id'], row['tree']
        assert key in expected and key not in seen
        seen.add(key)
        item = expected[key]
        path = source / item['path']
        assert row['source_fit'] == str(path) and row['source_sha256'] == item['sha256'] == sha(path)
        saved = json.loads(path.read_text())
        p = saved['payload']
        assert row['original_status'] == p['status']
        if p['status'] == 'ml_refinement_passed_numerical_checks':
            assert row['disposition'] == 'already_passed_no_proposal' and 'proposal' not in row
        elif p['checks'] != eligible:
            assert row['disposition'] == 'other_review_flags_preserved' and row['original_checks'] == p['checks'] and 'proposal' not in row
        else:
            input_path = Path(producer['inputs']) / (key[0] + '.npz')
            assert sha(input_path) == saved['input_sha256']
            with np.load(input_path) as a:
                values, bg, family, indices = a['matrix'], a['background'], a['family'], a['pattern_rows']
            with np.load(Path(producer['factors']) / (key[1] + '.npz')) as a:
                factor = a['factor'][indices]
            active = np.ptp(values[:, 1:], axis=0) > 1e-12
            scales = np.std(values[:, 1:][:, active], axis=0)
            np.testing.assert_array_equal(active, p['active_covariates'])
            np.testing.assert_array_equal(scales, p['covariate_scales'])
            x = np.column_stack((np.ones(len(values)), values[:, 1:][:, active] / scales))
            y = values[:, 0]
            cache = CachedMatchedML(bg, family, factor, x, y)
            def objective(theta):
                return profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(theta)), x, y)['negative_profiled_ml']
            def gradient(theta):
                return evaluate_gradient(cache, np.expm1(theta))['log1p_ratio_gradient']
            check_proposal(row['proposal'], objective, gradient, p['log1p_ratios'], np.log1p(p['maximum_ratio']))
            assert row['disposition'] == row['proposal']['status']
            proposals += 1
        counts[row['disposition']] += 1
    assert seen == set(expected) and len(seen) == receipt['fits']
    assert dict(counts) == receipt['counts']
    verify()
    result = dict(status='passed_full_audited_gradient_proposal_readback', fits=len(seen), proposal_cases=proposals, counts=dict(counts), source_receipt_sha256=sha(root / 'receipt.json'), checker_sha256=sha(__file__), source_hashes=bindings, producer_terminal_state=state, scope='All dispositions and proposals replayed without the proposal helper. Direct likelihood and analytic-gradient libraries are shared. No fit replacement, global optimum or calibrated inference claim.')
    with Path(plan['proof']).open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'source_hashes'}), flush=True)


if __name__ == '__main__':
    main()
