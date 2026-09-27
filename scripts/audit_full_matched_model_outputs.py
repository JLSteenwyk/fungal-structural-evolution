#!/usr/bin/env python3
"""Audit every completed model disposition and its parameter/unit bookkeeping."""
import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import psutil
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser()
    for name in ['plan','launch','output']:ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    launch=json.loads(args.launch.read_text());lh=sha(args.launch)
    assert launch['plan_sha256']==ph
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines())
    assert state=={'ActiveState':'inactive','ExecMainStatus':'0'},state
    assert sha(args.plan)==ph and sha(args.launch)==lh
    for path,h in plan['pins'].items():assert sha(path)==h,path
    root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['plan_sha256']==ph
    assert (root/'run_plan.json').read_bytes()==args.plan.read_bytes()
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h,name
    recipes={r['fit_input_id']:r for r in map(json.loads,(Path(plan['inventory'])/'unique_fit_recipes.jsonl').open())}
    trees={p.stem for p in Path(plan['factors']).glob('*.npz')};assert len(trees)==5
    expected={(identifier,tree) for identifier in recipes for tree in trees}
    seen=set();counts=Counter();candidate_count=0;maximum_direct_error=0.
    for line in (root/'fit_manifest.jsonl').open():
        entry=json.loads(line);key=entry['fit_input_id'],entry['tree'];assert key in expected and key not in seen;seen.add(key)
        path=Path(entry['path'])
        assert path==root/'fits'/key[0][:2]/key[0]/(key[1]+'.json')
        assert sha(path)==entry['sha256']
        result=json.loads(path.read_text());payload=result['payload']
        assert result['plan_sha256']==ph and result['fit_input_id']==key[0] and result['tree']==key[1]
        content=json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
        assert hashlib.sha256(content).hexdigest()==result['payload_sha256']
        assert entry['status']==payload['status'];counts[payload['status']]+=1
        if payload['status']=='fit_error_requires_review':
            assert payload['error_type'] and payload['traceback']
            continue
        assert payload['status'] in ['candidate_passed_numerical_optimization_checks','candidate_requires_optimization_review']
        assert payload['records']==recipes[key[0]]['records']
        candidates=payload['candidates'];assert len(candidates)==22
        expected_attempts=[]
        for face in itertools.product([False,True],repeat=3):
            expected_attempts.extend((list(face),start) for start in ([.05,1.,20.] if any(face) else [0.]))
        assert [(c['active'],c['start_ratio']) for c in candidates]==expected_attempts
        upper=np.log1p(payload['maximum_ratio'])
        for c in candidates:
            theta=np.array(c['theta']);active=np.array(c['active'])
            assert theta.shape==(3,) and np.isfinite(theta).all() and np.isfinite(c['objective'])
            assert (theta>=0).all() and (theta<=upper).all() and (theta[~active]==0).all()
        best=min(candidates,key=lambda c:c['objective'])
        np.testing.assert_allclose(payload['negative_profiled_reml'],best['objective'],rtol=1e-10,atol=1e-9)
        np.testing.assert_array_equal(payload['log1p_ratios'],best['theta'])
        np.testing.assert_allclose(payload['ratios'],np.expm1(best['theta']),rtol=1e-12,atol=1e-12)
        np.testing.assert_array_equal(payload['zero_boundary'],np.array(best['theta'])<=1e-7)
        np.testing.assert_array_equal(payload['upper_boundary'],np.array(best['theta'])>=upper-1e-6)
        flags=payload['checks'];gradient=np.array(payload['projected_gradient'])
        assert flags['best_optimizer_success']==best['success']
        assert flags['projected_gradient_pass']==bool(np.max(abs(gradient))<=1e-3)
        assert flags['upper_bound_contact']==any(payload['upper_boundary'])
        full=[c['objective'] for c in candidates if all(c['active'])]
        assert flags['all_full_face_starts_agree']==bool(np.ptp(full)<=1e-5)
        passed=flags['best_optimizer_success'] and flags['projected_gradient_pass'] and not flags['upper_bound_contact'] and flags['all_full_face_starts_agree']
        assert passed==(payload['status']=='candidate_passed_numerical_optimization_checks')
        active=np.array(payload['active_covariates']);assert active.shape==(4,)
        scales=np.array(payload['covariate_scales']);assert len(scales)==active.sum() and np.isfinite(scales).all() and (scales>0).all()
        coefficient_count=int(active.sum())+1
        beta=np.array(payload['beta']);cov=np.array(payload['conditional_beta_covariance'])
        assert beta.shape==(coefficient_count,) and cov.shape==(coefficient_count,coefficient_count)
        assert np.isfinite(beta).all() and np.isfinite(cov).all()
        np.testing.assert_allclose(cov,cov.T,rtol=1e-9,atol=1e-10)
        assert np.linalg.eigvalsh(cov).min()>0
        conversion=np.r_[1.,1./scales]
        np.testing.assert_allclose(payload['raw_unit_beta'],beta*conversion,rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(payload['raw_unit_conditional_beta_covariance'],cov*conversion[:,None]*conversion[None,:],rtol=1e-12,atol=1e-12)
        dof=payload['records']-coefficient_count
        assert payload['residual_degrees_of_freedom']==dof and dof>0
        np.testing.assert_allclose(payload['profiled_scale'],payload['residual_quadratic']/dof,rtol=1e-12,atol=1e-12)
        assert payload['profiled_scale']>0
        audit=payload['direct_candidate_readback'];assert audit['candidates_checked']==22
        assert 1<=audit['distinct_parameter_vectors']<=22
        candidate_count+=22;maximum_direct_error=max(maximum_direct_error,audit['maximum_objective_error'])
        if len(seen)%1000==0:print('Audited all output fields',len(seen),'/ 144040',flush=True)
    assert seen==expected and len(seen)==receipt['tree_fit_dispositions']==144040
    assert dict(counts)==receipt['status_counts']
    for path,h in plan['pins'].items():assert sha(path)==h,path
    result=dict(status='passed_full_working_model_output_integrity_audit',tree_fit_dispositions=len(seen),status_counts=dict(counts),
                reported_direct_candidate_checks=candidate_count,maximum_reported_direct_objective_error=maximum_direct_error,
                source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),
                scope='Every expected disposition, output checksum, source binding, face/start, bound/review classification and coefficient unit conversion checked. Producer direct numerical checks are counted, not independently repeated here. Error cases remain errors. This is output integrity, not calibration, biological model adequacy or final inference.')
    with args.output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')


if __name__=='__main__':main()
