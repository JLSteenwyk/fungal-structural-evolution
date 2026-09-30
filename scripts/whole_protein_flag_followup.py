"""Complete-scope source checks and separate numerical recovery helpers."""
import json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_ml_gradient import evaluate_gradient
from matched_mixed_covariance import MatchedCovariance
from refine_matched_ml_analytic import refine
from run_after_verified_dependencies import terminal_state
from screen_duplication_alignment_reuse import sha
from run_whole_protein_ml import digest

PASS = 'ml_candidate_passed_numerical_optimization_checks'
FLAG = 'ml_candidate_requires_optimization_review'
ERROR = 'fit_error_requires_review'


def full_scope(fits, audits, trees, input_ids):
    def index(rows):
        result = {}
        for row in rows:
            key = row['fit_input_id'], row['tree']
            if key in result or key[0] not in input_ids or key[1] not in trees:
                raise ValueError('Duplicate or unexpected full-grid identity')
            result[key] = row
        if len(result) != len(trees)*len(input_ids):
            raise ValueError('Incomplete full-grid scope')
        return result
    production, checked = index(fits), index(audits)
    flags, errors = [], []
    for key, row in production.items():
        audit = checked[key]
        status = row['status']
        if status not in (PASS, FLAG, ERROR) or audit['status'] != status:
            raise ValueError('Unknown or inconsistent audited status')
        if audit['numerical_fit_verified'] is not (status != ERROR):
            raise ValueError('Invalid numerical audit qualification')
        if status == FLAG:flags.append(row)
        if status == ERROR:errors.append(row)
    return production, flags, errors


def sources(config):
    bindings = dict(config['pins'])
    def bind(path, value=None):
        path = str(path); value = sha(path) if value is None else value
        if path in bindings and bindings[path] != value:
            raise ValueError('Conflicting source hash: '+path)
        bindings[path] = value
    fit_plan = json.loads(Path(config['fit_plan']).read_text())
    audit_plan = json.loads(Path(config['audit_plan']).read_text())
    assert audit_plan['fit_launch'] == config['fit_launch']
    for name in ('fit_launch', 'audit_launch'):
        launch = json.loads(Path(config[name]).read_text())
        expected_plan = config['fit_plan' if name == 'fit_launch' else 'audit_plan']
        assert launch['plan'] == expected_plan and launch['plan_sha256'] == sha(expected_plan)
        assert terminal_state(launch['unit']) == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    root = Path(fit_plan['output']); audit = Path(audit_plan['output'])
    fit_receipt = json.loads((root/'receipt.json').read_text())
    ar = json.loads((audit/'receipt.json').read_text())
    assert fit_receipt['status'] == 'complete_whole_protein_ml_dispositions_pending_full_audit'
    assert fit_receipt['plan_sha256'] == sha(config['fit_plan'])
    assert ar['status'] == 'complete_full_whole_protein_ml_output_audit'
    assert ar['source_receipt_sha256'] == sha(root/'receipt.json')
    assert ar['audit_plan_sha256'] == sha(config['audit_plan'])
    assert ar['status_counts'] == fit_receipt['status_counts']
    assert ar['tree_fit_dispositions'] == fit_receipt['tree_fit_dispositions'] == 375350
    assert fit_receipt['unique_inputs'] == 75070
    for folder, receipt in ((root, fit_receipt), (audit, ar)):
        bind(folder/'receipt.json')
        for name, h in receipt['artifacts'].items():bind(folder/name, h)
    for path, h in json.loads((root/'source_bindings.json').read_text()).items():bind(path, h)
    for plan in (fit_plan, audit_plan):
        for path, h in plan['pins'].items():bind(path, h)
    inputs = Path(fit_plan['inputs']); recipes = {}
    ir = json.loads((inputs/'receipt.json').read_text())
    bind(inputs/'receipt.json')
    for name, h in ir['artifacts'].items():bind(inputs/name, h)
    for path, h in bindings.items():assert sha(path) == h, path
    with (inputs/'input_manifest.jsonl').open() as f:
        for line in f:
            row = json.loads(line); key = row['fit_input_id']
            assert key not in recipes; recipes[key] = row
    assert len(recipes) == 75070
    factor_root = Path(fit_plan['factors'])
    trees = {p.stem for p in factor_root.glob('*.npz')}; assert len(trees) == 5
    with (root/'fit_manifest.jsonl').open() as f, (audit/'audit_manifest.jsonl').open() as a:
        production, flags, errors = full_scope((json.loads(l) for l in f),
            (json.loads(l) for l in a), trees, set(recipes))
    assert dict(Counter(r['status'] for r in production.values())) == fit_receipt['status_counts']
    for row in flags:bind(row['path'],row['sha256'])
    # Bind every flagged input's full audit record to its original five source fits.
    for identifier in {r['fit_input_id'] for r in flags}:
        bind(inputs/recipes[identifier]['path'],recipes[identifier]['sha256'])
        path = audit/'inputs'/identifier[:2]/(identifier+'.json'); bind(path)
        proof = json.loads(path.read_text()); identity = proof['identity']
        assert identity['fit_input_id'] == identifier and identity['input_sha256'] == recipes[identifier]['sha256']
        assert identity['fit_plan_sha256'] == sha(config['fit_plan']) and identity['audit_plan_sha256'] == sha(config['audit_plan'])
        assert identity['fit_hashes'] == {production[identifier,t]['path']:production[identifier,t]['sha256'] for t in trees}
        assert proof['results_sha256'] == digest(proof['results'])
        assert {r['tree'] for r in proof['results']} == trees and len(proof['results']) == 5
        for row in proof['results']:
            assert row['status'] == production[identifier,row['tree']]['status']
            assert row['numerical_fit_verified'] is (row['status'] != ERROR)
    for path,h in bindings.items():assert sha(path)==h,path
    return fit_plan, production, flags, errors, recipes, bindings


def arrays(item, recipe, fit_plan):
    path = Path(fit_plan['inputs'])/recipe['path']
    assert sha(path) == recipe['sha256']
    assert sha(item['path']) == item['sha256']
    original = json.loads(Path(item['path']).read_text())
    assert original['fit_input_id'] == item['fit_input_id'] and original['tree'] == item['tree']
    assert original['input_sha256'] == recipe['sha256']
    assert original['specification'] == recipe['recipe']['specification']
    assert original['payload']['status'] == FLAG and original['payload_sha256'] == digest(original['payload'])
    with np.load(path, allow_pickle=False) as a:
        matrix, bg, family, rows = a['matrix'], a['background'], a['family'], a['pattern_rows']
    with np.load(Path(fit_plan['factors'])/(item['tree']+'.npz'), allow_pickle=False) as a:
        assert np.all(rows >= 0) and np.all(rows < len(a['factor']))
        factor = a['factor'][rows]
    assert np.all(matrix[:,1] == 1)
    scales = np.std(matrix[:,2:], axis=0)
    assert np.all(scales > 0)
    np.testing.assert_array_equal(scales, original['payload']['covariate_scales'])
    x = np.column_stack([matrix[:,1], matrix[:,2:]/scales])
    return original, matrix, bg, family, factor, x, matrix[:,0], scales


def recover(refinement, bg, family, factor, x, y):
    cache = CachedMatchedML(bg, family, factor, x, y)
    upper = float(np.log1p(refinement['maximum_ratio']))
    def objective(theta):
        ratios = np.expm1(theta)
        value = profiled_ml(MatchedCovariance(bg, family, factor, 1., *ratios), x, y)['negative_profiled_ml']
        gradient = evaluate_gradient(cache, ratios)
        np.testing.assert_allclose(value, gradient['negative_profiled_ml'], rtol=1e-9, atol=1e-7)
        return value, gradient['log1p_ratio_gradient']
    def describe(theta):
        theta = np.asarray(theta); value, gradient = objective(theta)
        projected = np.where(theta <= 1e-7, np.minimum(gradient, 0),
            np.where(theta >= upper-1e-7, np.maximum(gradient, 0), gradient))
        return dict(theta=theta.tolist(), objective=float(value), gradient=gradient.tolist(),
            projected_gradient=projected.tolist(), maximum_absolute_projected_gradient=float(max(abs(projected))))
    starts = [(i,c) for i,c in enumerate(refinement['candidates'])
        if all(c['active']) and c['start'] != 'original_unoptimized_reference']
    assert len(starts) == 4
    starts.append((None, dict(start='1.0_original_initialization', theta=[float(np.log(2))]*3)))
    runs = []
    for index, candidate in starts:
        before = describe(candidate['theta'])
        fit = minimize(objective, np.asarray(candidate['theta']), jac=True, method='SLSQP',
            bounds=[(0., upper)]*3, options=dict(ftol=1e-12, maxiter=1000))
        runs.append(dict(source_candidate_index=index, start=candidate['start'], initial=before,
            endpoint=describe(fit.x), success=bool(fit.success), message=str(fit.message),
            iterations=int(fit.nit), evaluations=int(fit.nfev)))
    choices = [(c['objective'],c['success'],'original_candidate',i,c['theta']) for i,c in enumerate(refinement['candidates'])]
    choices += [(r['endpoint']['objective'],r['success'],'recovery_endpoint',i,r['endpoint']['theta']) for i,r in enumerate(runs)]
    selected = min(choices, key=lambda r:(r[0], not r[1])); summary = describe(selected[4])
    values = [r['endpoint']['objective'] for r in runs]
    checks = dict(all_recovery_optimizers_successful=all(r['success'] for r in runs),
        all_recovery_projected_gradients_pass=all(r['endpoint']['maximum_absolute_projected_gradient'] <= 1e-3 for r in runs),
        all_full_face_recovered_starts_agree=bool(np.ptp(values) <= 1e-5),
        selected_optimizer_success=bool(selected[1]), selected_projected_gradient_pass=summary['maximum_absolute_projected_gradient'] <= 1e-3,
        no_recovery_upper_bound_contact=all(max(r['endpoint']['theta']) < upper-1e-6 for r in runs),
        selected_no_upper_bound_contact=bool(max(selected[4]) < upper-1e-6),
        original_objective_not_worsened=bool(selected[0] <= refinement['negative_profiled_ml']),
        recovered_starts_agree_with_selected=bool(max(abs(v-summary['objective']) for v in values) <= 1e-5))
    fitted = profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(selected[4])), x, y)
    return dict(original_candidates=refinement['candidates'], recoveries=runs,
        selection=dict(kind=selected[2],index=selected[3],**summary), checks=checks,
        fitted={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in fitted.items()},
        optimizer=dict(method='SLSQP',ftol=1e-12,maxiter=1000,bounds=[0.,upper]))


def followup(original, bg, family, factor, x, y, scales):
    refined = refine(bg, family, factor, x, y, original['payload']['log1p_ratios'],
        maximum_ratio=10000., maxiter=1000)
    recovery = None
    if refined['status'] != 'ml_refinement_passed_numerical_checks':
        recovery = recover(refined, bg, family, factor, x, y)
    passed = refined['status'] == 'ml_refinement_passed_numerical_checks' if recovery is None else all(recovery['checks'].values())
    theta = refined['log1p_ratios'] if recovery is None else recovery['selection']['theta']
    fitted = profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(theta)), x, y)
    conversion = np.r_[1.,1/scales]
    return dict(status='numerical_followup_passed_pending_independent_readback' if passed else 'numerical_followup_requires_review',
        refinement=refined,recovery=recovery,selected_theta=theta,
        fitted={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in fitted.items()},
        raw_unit_beta=(fitted['beta']*conversion).tolist(),
        raw_unit_conditional_beta_covariance=(fitted['conditional_beta_covariance']*np.outer(conversion,conversion)).tolist())
