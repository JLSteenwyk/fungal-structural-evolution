#!/usr/bin/env python3
"""Validate candidate selection and recompute GLS/ML without profiled_ml."""
import json
import subprocess
import time
from collections import Counter
from pathlib import Path
import numpy as np
import psutil
from cached_matched_ml import CachedMatchedML
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from screen_duplication_alignment_reuse import sha


def independent_fit(cov,x,y):
    precision_x=cov.solve(x);precision_y=cov.solve(y)
    info=x.T@precision_x
    info=(info+info.T)/2
    beta=np.linalg.solve(info,x.T@precision_y)
    residual=y-x@beta
    quadratic=float(residual@cov.solve(residual));scale=quadratic/len(y)
    assert scale>0
    return dict(beta=beta,conditional_beta_covariance=scale*np.linalg.inv(info),residual_quadratic=quadratic,profiled_scale=scale,negative_profiled_ml=float((cov.logdet+len(y)*(np.log(2*np.pi*scale)+1))/2))


def main():
    pp=Path('metadata/selected_refinement_candidates_readback_plan_20260929.json')
    plan=json.loads(pp.read_text());bindings={str(pp):sha(pp),**plan['pins']}
    launch=json.loads(Path(plan['launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p=psutil.Process(launch['pid'])
            if abs(p.create_time()-launch['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    producer=json.loads(Path(launch['plan']).read_text());root=Path(producer['output'])
    receipt=json.loads((root/'receipt.json').read_text());assert receipt['status']=='complete_selected_refinement_candidates_pending_readback'
    bindings.update(receipt['source_hashes']);bindings[str(root/'receipt.json')]=sha(root/'receipt.json')
    bindings.update({str(root/n):h for n,h in receipt['artifacts'].items()})
    def verify():
        for path,digest in bindings.items():assert sha(path)==digest,path
    verify()
    source=Path(producer['refinements'])
    originals={(r['fit_input_id'],r['tree']):r for r in json.loads((source/'fit_manifest.json').read_text())}
    proposals=[]
    for directory in producer['proposals']:
        rows=[json.loads(line) for line in (Path(directory)/'dispositions.jsonl').open()]
        proposals.append({(r['fit_input_id'],r['tree']):r for r in rows})
    seen=set();counts=Counter();maximum_error=0.
    for item in json.loads((root/'manifest.json').read_text()):
        key=item['fit_input_id'],item['tree'];assert key not in seen and key in originals;seen.add(key)
        path=root/item['path'];assert sha(path)==item['sha256'];r=json.loads(path.read_text())
        oldpath=source/originals[key]['path'];saved=json.loads(oldpath.read_text());old=saved['payload']
        assert r['source_refinement']==str(oldpath) and r['source_refinement_sha256']==sha(oldpath)
        assert r['original_production_fit']==saved['source_fit_path'] and r['original_production_sha256']==saved['source_fit_sha256']
        assert r['source_refinement_status']==old['status'] and r['source_refinement_checks']==old['checks']
        assert (r['fit_input_id'],r['tree'])==key and r['polynomial_degree']==saved['polynomial_degree']
        if old['status']=='ml_refinement_passed_numerical_checks':
            theta=np.array(old['log1p_ratios']);kind='audited_refinement'
            assert r['proposal_index'] is None and r['hessian_step'] is None
        else:
            possible=[(i,p[key]) for i,p in enumerate(proposals) if p[key]['disposition'] in ['two_local_proposals_passed','two_lower_face_proposals_passed']]
            assert len(possible)==1
            i,entry=possible[0];trials=entry['proposal']['proposals'];assert len(trials)==2 and all(t['passed'] for t in trials)
            objectives=[t['objective'] for t in trials];assert max(objectives)-min(objectives)<=1e-5
            chosen=sorted(trials,key=lambda t:(t['objective'],t['hessian_step']))[0]
            theta=np.array(chosen['theta']);kind='interior_local_proposal' if i==0 else 'lower_face_local_proposal'
            assert r['proposal_index']==i and r['hessian_step']==chosen['hessian_step']
        assert r['selection_kind']==item['selection_kind']==kind
        np.testing.assert_array_equal(r['log1p_ratios'],theta)
        assert r['maximum_ratio']==old['maximum_ratio']
        assert (theta>=0).all() and (theta<np.log1p(r['maximum_ratio'])-1e-6).all()
        input_path=Path(producer['inputs'])/(key[0]+'.npz');assert sha(input_path)==r['input_sha256']==saved['input_sha256']
        with np.load(input_path) as a:values,bg,family,rows=a['matrix'],a['background'],a['family'],a['pattern_rows']
        with np.load(Path(producer['factors'])/(key[1]+'.npz')) as a:factor=a['factor'][rows]
        active=np.max(values[:,1:],axis=0)-np.min(values[:,1:],axis=0)>1e-12
        scales=np.std(values[:,1:][:,active],axis=0)
        np.testing.assert_array_equal(r['active_covariates'],active);np.testing.assert_array_equal(r['covariate_scales'],scales)
        x=np.column_stack([np.ones(len(values)),values[:,1:][:,active]/scales]);y=values[:,0]
        ratios=np.expm1(theta);np.testing.assert_array_equal(r['ratios'],ratios)
        fit=independent_fit(MatchedCovariance(bg,family,factor,1.,*ratios),x,y)
        for field,value in fit.items():np.testing.assert_allclose(r['fitted'][field],value,rtol=1e-8,atol=1e-8)
        maximum_error=max(maximum_error,abs(r['fitted']['negative_profiled_ml']-fit['negative_profiled_ml']))
        assert r['records']==r['fitted']['variance_profile_denominator']==len(y)
        assert r['fitted']['residual_degrees_of_freedom']==len(y)-x.shape[1]
        conversion=np.r_[1.,1/scales]
        np.testing.assert_allclose(r['raw_unit_beta'],fit['beta']*conversion,rtol=1e-8,atol=1e-8)
        np.testing.assert_allclose(r['raw_unit_conditional_beta_covariance'],fit['conditional_beta_covariance']*np.outer(conversion,conversion),rtol=1e-8,atol=1e-8)
        baseline=independent_fit(MatchedCovariance(bg,family,factor,1.,*np.expm1(old['log1p_ratios'])),x,y)
        delta=fit['negative_profiled_ml']-baseline['negative_profiled_ml']
        np.testing.assert_allclose(r['objective_change_vs_refinement'],delta,rtol=1e-8,atol=1e-8)
        assert r['objective_change_vs_refinement']<=0 and delta<=1e-8
        grad=evaluate_gradient(CachedMatchedML(bg,family,factor,x,y),ratios)['log1p_ratio_gradient']
        projected=np.where(theta<=1e-7,np.minimum(grad,0),grad)
        np.testing.assert_allclose(r['gradient'],grad,rtol=1e-9,atol=1e-9)
        np.testing.assert_allclose(r['projected_gradient'],projected,rtol=1e-9,atol=1e-9)
        assert max(abs(projected))<=1e-3
        if kind=='lower_face_local_proposal':
            fixed=np.asarray(old['log1p_ratios'])==0;assert np.all(theta[fixed]==0) and np.all(grad[fixed]>=0)
        assert r['status']=='selected_numerical_candidate_pending_readback'
        counts[kind]+=1
        if len(seen)%50==0:print('Checked selected candidates',len(seen),'/601',flush=True)
    assert seen==set(originals) and len(seen)==receipt['fits']==601 and dict(counts)==receipt['selection_counts']
    verify()
    proof=dict(status='passed_full_selected_refinement_candidate_readback',fits=len(seen),selection_counts=dict(counts),maximum_likelihood_error=maximum_error,source_receipt_sha256=sha(root/'receipt.json'),checker_sha256=sha(__file__),producer_terminal_state=state,scope='All selections and GLS/ML summaries reconstructed. Separate normal-equation solve/inverse implementation shares covariance solver and analytic-gradient library. No original fit overwrite, full-grid acceptance, calibrated uncertainty or biological inference.')
    with Path(plan['proof']).open('x') as h:json.dump(proof,h,indent=2);h.write('\n')
    print(json.dumps(proof),flush=True)


if __name__=='__main__':main()
