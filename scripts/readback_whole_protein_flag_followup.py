#!/usr/bin/env python3
"""Independently replay all saved full-grid numerical follow-up candidates."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import fcntl
import itertools
import json
import multiprocessing
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from cached_matched_ml import CachedMatchedML
from matched_ml_gradient import evaluate_gradient
from matched_mixed_covariance import MatchedCovariance
from readback_selected_refinement_candidates import independent_fit
from whole_protein_flag_followup import sources, arrays, FLAG, ERROR
from run_after_verified_dependencies import terminal_state
from run_whole_protein_ml import digest
from screen_duplication_alignment_reuse import sha


def check_result(result, original, bg, family, factor, x, y, scales, finite_difference_policy='legacy'):
    if result['status']=='numerical_followup_error_requires_review':
        assert result['error_type'] and result['traceback']
        return dict(status=result['status'],numerical_result_verified=False,candidates_checked=0,
            maximum_objective_error=0.,finite_difference_checks=[],
            independent_readback_status='followup_error_requires_review')
    p=result['refinement'];candidates=p['candidates'];upper=float(np.log1p(10000.))
    assert p['maximum_ratio']==original['payload']['maximum_ratio']==10000.
    grid=[([True]*3,'original_unoptimized_reference')]
    for face in itertools.product([False,True],repeat=3):
        grid.extend((list(face),start) for start in (['0.05','1.0','20.0'] if any(face) else ['exact_zero']))
        if all(face):grid.append((list(face),'original_parameters'))
    assert [(c['active'],c['start']) for c in candidates]==grid and len(candidates)==24
    np.testing.assert_array_equal(candidates[0]['theta'],original['payload']['log1p_ratios'])
    assert not candidates[0]['success']
    cache=CachedMatchedML(bg,family,factor,x,y);maximum=0.;checked=0
    def independent(theta):
        theta=np.asarray(theta);assert theta.shape==(3,) and np.isfinite(theta).all()
        assert (theta>=0).all() and (theta<=upper).all()
        return independent_fit(MatchedCovariance(bg,family,factor,1.,*np.expm1(theta)),x,y)
    def objective(theta,value):
        nonlocal maximum,checked
        fit=independent(theta);error=abs(fit['negative_profiled_ml']-value)
        np.testing.assert_allclose(fit['negative_profiled_ml'],value,rtol=1e-9,atol=1e-7)
        maximum=max(maximum,error);checked+=1;return fit
    def gradient(theta):
        theta=np.asarray(theta)
        grad=evaluate_gradient(cache,np.expm1(theta))['log1p_ratio_gradient']
        projected=np.where(theta<=1e-7,np.minimum(grad,0),np.where(theta>=upper-1e-7,np.maximum(grad,0),grad))
        return grad,projected
    for c in candidates:
        theta=np.asarray(c['theta']);assert (theta[~np.asarray(c['active'])]==0).all()
        objective(theta,c['objective'])
    best=min(candidates,key=lambda c:(c['objective'],not c['success']))
    theta=np.asarray(best['theta']);np.testing.assert_array_equal(p['log1p_ratios'],theta)
    np.testing.assert_allclose(p['ratios'],np.expm1(theta),rtol=1e-12,atol=1e-12)
    grad,projected=gradient(theta)
    np.testing.assert_allclose(p['analytic_gradient'],grad,rtol=1e-9,atol=1e-9)
    np.testing.assert_allclose(p['projected_gradient'],projected,rtol=1e-9,atol=1e-9)
    full=[c['objective'] for c in candidates if all(c['active']) and c['start']!='original_unoptimized_reference']
    checks=dict(best_optimizer_success=best['success'],projected_gradient_pass=bool(max(abs(projected))<=1e-3),
        upper_bound_contact=bool(max(theta)>=upper-1e-6),all_full_face_starts_agree=max(full)-min(full)<=1e-5,
        original_objective_not_worsened=best['objective']<=candidates[0]['objective'])
    assert checks==p['checks']
    passed=all(v for k,v in checks.items() if k!='upper_bound_contact') and not checks['upper_bound_contact']
    assert p['status']==('ml_refinement_passed_numerical_checks' if passed else 'ml_refinement_requires_review')
    assert p['original_objective']==candidates[0]['objective']
    np.testing.assert_allclose(p['original_objective'],original['payload']['negative_profiled_ml'],rtol=1e-9,atol=1e-7)
    assert p['objective_improvement']==p['original_objective']-best['objective']
    fit=independent(theta)
    for k,v in fit.items():np.testing.assert_allclose(p[k],v,rtol=1e-8,atol=1e-8)
    assert p['variance_profile_denominator']==len(y) and p['residual_degrees_of_freedom']==len(y)-x.shape[1]
    assert p['likelihood']=='ordinary_gaussian_ml'
    points=[('selected',theta,grad)]
    r=result['recovery']
    if passed:assert r is None
    else:
        assert r is not None and r['original_candidates']==candidates
        assert r['optimizer']==dict(method='SLSQP',ftol=1e-12,maxiter=1000,bounds=[0.,upper])
        expected=[(i,c['start'],c['theta']) for i,c in enumerate(candidates)
            if all(c['active']) and c['start']!='original_unoptimized_reference']
        expected.append((None,'1.0_original_initialization',[float(np.log(2))]*3))
        runs=r['recoveries'];assert len(runs)==len(expected)==5
        def replay(row):
            objective(row['theta'],row['objective']);g,pg=gradient(row['theta'])
            np.testing.assert_allclose(row['gradient'],g,rtol=1e-9,atol=1e-9)
            np.testing.assert_allclose(row['projected_gradient'],pg,rtol=1e-9,atol=1e-9)
            assert row['maximum_absolute_projected_gradient']==float(max(abs(pg)))
        for run,(index,label,start) in zip(runs,expected):
            assert run['source_candidate_index']==index and run['start']==label
            np.testing.assert_array_equal(run['initial']['theta'],start)
            assert isinstance(run['success'],bool) and run['message'] and run['evaluations']>0 and run['iterations']>=0
            replay(run['initial']);replay(run['endpoint'])
        choices=[(c['objective'],c['success'],'original_candidate',i,c['theta']) for i,c in enumerate(candidates)]
        choices += [(v['endpoint']['objective'],v['success'],'recovery_endpoint',i,v['endpoint']['theta']) for i,v in enumerate(runs)]
        choice=min(choices,key=lambda v:(v[0],not v[1]));selection=r['selection']
        assert selection['kind']==choice[2] and selection['index']==choice[3]
        np.testing.assert_array_equal(selection['theta'],choice[4]);replay(selection)
        values=[v['endpoint']['objective'] for v in runs]
        rc=dict(all_recovery_optimizers_successful=all(v['success'] for v in runs),
            all_recovery_projected_gradients_pass=all(v['endpoint']['maximum_absolute_projected_gradient']<=1e-3 for v in runs),
            all_full_face_recovered_starts_agree=max(values)-min(values)<=1e-5,
            selected_optimizer_success=choice[1],selected_projected_gradient_pass=selection['maximum_absolute_projected_gradient']<=1e-3,
            no_recovery_upper_bound_contact=all(max(v['endpoint']['theta'])<upper-1e-6 for v in runs),
            selected_no_upper_bound_contact=max(choice[4])<upper-1e-6,
            original_objective_not_worsened=choice[0]<=p['negative_profiled_ml'],
            recovered_starts_agree_with_selected=max(abs(v-selection['objective']) for v in values)<=1e-5)
        assert rc==r['checks'];passed=all(rc.values());theta=np.asarray(selection['theta'])
        fitted=independent(theta)
        for k,v in fitted.items():np.testing.assert_allclose(r['fitted'][k],v,rtol=1e-8,atol=1e-8)
        points=[('selected',theta,np.asarray(selection['gradient']))]+[(v['start'],np.asarray(v['initial']['theta']),np.asarray(v['initial']['gradient']))
            for v in runs if v['initial']['maximum_absolute_projected_gradient']>1e-3]
    np.testing.assert_array_equal(result['selected_theta'],theta)
    assert result['status']==('numerical_followup_passed_pending_independent_readback' if passed else 'numerical_followup_requires_review')
    fitted=independent(theta)
    for k,v in fitted.items():np.testing.assert_allclose(result['fitted'][k],v,rtol=1e-8,atol=1e-8)
    assert result['fitted']['variance_profile_denominator']==len(y) and result['fitted']['residual_degrees_of_freedom']==len(y)-x.shape[1]
    conversion=np.r_[1.,1/scales]
    np.testing.assert_allclose(result['raw_unit_beta'],fitted['beta']*conversion,rtol=1e-9,atol=1e-9)
    np.testing.assert_allclose(result['raw_unit_conditional_beta_covariance'],fitted['conditional_beta_covariance']*np.outer(conversion,conversion),rtol=1e-9,atol=1e-9)
    finite=[]
    if finite_difference_policy not in {'legacy', 'adaptive_v2'}:
        raise ValueError('unknown finite-difference policy')
    for label,point,g in points:
        for step in ((1e-5,1e-6) if finite_difference_policy=='legacy' else (1e-6,)):
            values=[]
            for axis in range(3):
                d=np.eye(3)[axis]*step
                value=lambda t:independent(t)['negative_profiled_ml']
                if point[axis]<step:fd=(-3*value(point)+4*value(point+d)-value(point+2*d))/(2*step)
                elif point[axis]>upper-step:fd=(3*value(point)-4*value(point-d)+value(point-2*d))/(2*step)
                else:fd=(value(point+d)-value(point-d))/(2*step)
                values.append(fd)
            maximum_error=float(max(abs(np.asarray(values)-g)))
            if finite_difference_policy=='legacy':
                np.testing.assert_allclose(values,g,rtol=1e-4,atol=1e-3)
            finite.append(dict(point=label,step=step,maximum_absolute_error=maximum_error,
                               within_legacy_tolerance=bool(np.allclose(values,g,rtol=1e-4,atol=1e-3))))
    finite_pass=all(row['within_legacy_tolerance'] for row in finite)
    return dict(status=result['status'],numerical_result_verified=finite_pass,candidates_checked=checked,
        maximum_objective_error=maximum,finite_difference_checks=finite,
        independent_readback_status=('independently_replayed_numerical_followup' if finite_pass else 'finite_difference_discrepancy_requires_review'))


def work(task):
    row,item,recipe,fit_plan,root,finite_difference_policy=task
    path=Path(root)/row['path'];assert sha(path)==row['sha256'];saved=json.loads(path.read_text())
    original,matrix,bg,family,factor,x,y,scales=arrays(item,recipe,fit_plan)
    identity=dict(fit_input_id=item['fit_input_id'],tree=item['tree'],original_fit=item['path'],
        original_fit_sha256=item['sha256'],input_sha256=recipe['sha256'],plan_sha256=fit_plan['followup_plan_sha256'])
    assert saved['identity']==identity and saved['result_sha256']==digest(saved['result'])
    assert saved['original_status']==FLAG and saved['original_checks']==original['payload']['checks']
    assert saved['columns']==original['specification']['columns'] and saved['records']==len(matrix)
    assert saved['scientific_eligibility'] is False
    np.testing.assert_array_equal(saved['covariate_scales'],scales)
    assert saved['result']['status']==row['status']
    checked = dict(fit_input_id=item['fit_input_id'],tree=item['tree'],source_sha256=row['sha256'],
        **check_result(saved['result'],original,bg,family,factor,x,y,scales,finite_difference_policy))
    assert sha(path)==row['sha256']
    return checked


def initialize():threadpool_limits(limits=1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',required=True,type=Path)
    parser.add_argument('--finite-difference-policy',choices=['legacy','adaptive_v2'],default='legacy')
    args=parser.parse_args();config=json.loads(args.plan.read_text())
    producer=json.loads(Path(config['source_plan']).read_text());ph=sha(config['source_plan'])
    launch=json.loads(Path(config['source_launch']).read_text())
    assert sha(launch['plan'])==launch['plan_sha256']
    state=terminal_state(launch['unit'])
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    historical_code_pins=config.get('historical_source_code_pins',{})
    fit_plan,production,flags,errors,recipes,bindings=sources(producer,historical_code_pins)
    for p,h in historical_code_pins.items():
        assert bindings[p]==h,p
    for p,h in config['pins'].items():
        assert p not in bindings or bindings[p]==h;bindings[p]=h
    # Historical producer code is retained in ``bindings``.  Record the
    # independently executed readback implementation separately, so a
    # superseding validator cannot silently masquerade as producer code.
    implementation_pins=config['implementation_pins']
    for p,h in implementation_pins.items():
        assert sha(p)==h,p
    bindings[str(args.plan)]=sha(args.plan)
    root=Path(producer['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_full_grid_optimization_followup_pending_independent_readback'
    assert receipt['plan_sha256']==ph and receipt['full_dispositions']==375350
    assert receipt['flagged_fits']==len(flags) and receipt['original_fit_errors_retained']==len(errors)
    assert receipt['scientific_eligibility'] is False
    bindings[str(root/'receipt.json')]=sha(root/'receipt.json')
    for n,h in receipt['artifacts'].items():bindings[str(root/n)]=h
    def verify():
        for p,h in bindings.items():
            if p not in historical_code_pins:
                assert sha(p)==h,p
    verify()
    seen=set()
    with (root/'scope_dispositions.jsonl').open() as f:
        for line in f:
            row=json.loads(line);key=row['fit_input_id'],row['tree'];assert key in production and key not in seen;seen.add(key)
            expected=production[key];scope=('optimization_flag_followup_required' if expected['status']==FLAG else
                'original_fit_error_retained_requires_review' if expected['status']==ERROR else 'original_numerical_pass_retained')
            assert row==dict(**expected,followup_scope_status=scope,scientific_eligibility=False)
    assert seen==set(production)
    expected={(r['fit_input_id'],r['tree']):r for r in flags};seen=set();rows=[]
    with (root/'followup_manifest.jsonl').open() as f:
        for line in f:
            r=json.loads(line);key=r['fit_input_id'],r['tree'];assert key in expected and key not in seen;seen.add(key);rows.append(r)
    assert seen==set(expected)
    fit_plan['followup_plan_sha256']=ph
    output=Path(config['output']);output.mkdir(parents=True,exist_ok=True)
    lock=(output/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    stored=output/'run_plan.json'
    if stored.exists():assert stored.read_bytes()==args.plan.read_bytes()
    else:stored.write_bytes(args.plan.read_bytes())
    counts=Counter();independent_counts=Counter();maximum=0.;total=candidates=0
    with (output/'readback_manifest.jsonl').open('w') as f,ProcessPoolExecutor(max_workers=config['workers'],
        mp_context=multiprocessing.get_context('spawn'),initializer=initialize) as pool:
        tasks=((r,expected[r['fit_input_id'],r['tree']],recipes[r['fit_input_id']],fit_plan,str(root),args.finite_difference_policy) for r in rows)
        pending=set()
        def collect():
            nonlocal pending,total,candidates,maximum
            done,pending=wait(pending,return_when=FIRST_COMPLETED)
            for future in done:
                result=future.result();counts[result['status']]+=1;independent_counts[result['independent_readback_status']]+=1;total+=1
                candidates+=result['candidates_checked'];maximum=max(maximum,result['maximum_objective_error'])
                f.write(json.dumps(result)+'\n')
            f.flush();print('Independent full-grid flag check',total,'/',len(flags),flush=True)
        for task in tasks:
            pending.add(pool.submit(work,task))
            if len(pending)>=2*config['workers']:collect()
        while pending:collect()
    assert total==len(flags) and dict(counts)==receipt['followup_status_counts']
    verify()
    result=dict(status='passed_full_grid_optimization_followup_readback_with_review_statuses_retained',
        full_dispositions=len(production),flagged_fits=total,original_fit_errors_retained=len(errors),counts=dict(counts),
        independent_readback_counts=dict(independent_counts),finite_difference_policy=args.finite_difference_policy,
        candidate_likelihoods_replayed=candidates,maximum_objective_error=maximum,
        source_receipt_sha256=sha(root/'receipt.json'),producer_terminal_state=state,source_hashes=bindings,
        implementation_hashes=implementation_pins,
        artifacts={'readback_manifest.jsonl':sha(output/'readback_manifest.jsonl')},scientific_eligibility=False,
        scope='Full 375350-disposition census and every flagged serialized result checked. All candidate LLs independently reconstructed by GLS normal equations; selected/raw coefficients and covariance transforms, gradients, recovery lineage and finite-difference scores checked. Original fit and follow-up errors retain unverified numerical status. Shared covariance operator/analytic-score library; no global optimum, inferential calibration or biological acceptance.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
