#!/usr/bin/env python3
"""Validate polynomial worker serialization, direct ML, scaling and restart behavior."""
import json, tempfile
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from threadpoolctl import threadpool_limits
import run_full_polynomial_ml as runner
from screen_duplication_domain_alignment_coverage import sha


def main():
    threadpool_limits(1)
    rng=np.random.default_rng(9272026);n=80
    bg=np.repeat(np.arange(20),4);family=bg//4;factor=rng.normal(size=(n,3))*.1
    runner.FACTORS={f'tree{i}':factor*(i/4) for i in range(5)}
    u=rng.uniform(.1,.9,n);v=rng.uniform(.1,.9,n)
    cov=np.column_stack([u-v,rng.normal(size=n),rng.normal(size=n),np.zeros(n),u*u-v*v,u**3-v**3])
    y=.2+.3*cov[:,0]+rng.normal(size=n)*.4
    matrix=np.column_stack([y,cov]);checked=0;fit_counts={};maxerr=0.
    with tempfile.TemporaryDirectory(prefix='polynomial-ml-check-') as output:
        for degree in [1,2,3]:
            task=(str(degree)*64,degree,bg,family,np.arange(n),matrix[:,:degree+4],output,'fixture-plan')
            rows=runner.fit_case(task);assert len(rows)==5
            assert runner.fit_case(task)==rows
            for row in rows:
                saved=json.loads(Path(row['path']).read_text());p=saved['payload']
                assert p['status']!='fit_error_requires_review',p
                assert saved['polynomial_degree']==p['polynomial_degree']==degree
                assert p['likelihood']=='ordinary_gaussian_ml'
                assert p['variance_profile_denominator']==n and len(p['candidates'])==22
                active=np.array(p['active_covariates']);xraw=np.column_stack([np.ones(n),cov[:,:degree+3][:,active]])
                x=xraw.copy();x[:,1:]/=p['covariate_scales']
                np.testing.assert_allclose(x@p['beta'],xraw@p['raw_unit_beta'],atol=1e-12)
                f=runner.FACTORS[row['tree']]
                for c in p['candidates']:
                    b,g,s=np.expm1(c['theta']);r=np.eye(n)+b*(bg[:,None]==bg)+g*(family[:,None]==family)+s*(f@f.T)
                    ch=cho_factor(r);rx=cho_solve(ch,x);ry=cho_solve(ch,y)
                    beta=np.linalg.solve(x.T@rx,x.T@ry);res=y-x@beta;q=res@cho_solve(ch,res)
                    objective=.5*(2*np.log(np.diag(ch[0])).sum()+n*(1+np.log(2*np.pi*q/n)))
                    np.testing.assert_allclose(objective,c['objective'],rtol=1e-9,atol=1e-7)
                    maxerr=max(maxerr,abs(objective-c['objective']));checked+=1
                fit_counts[p['status']]=fit_counts.get(p['status'],0)+1
            try:runner.fit_case((*task[:-1],'changed-plan'))
            except AssertionError:pass
            else:raise AssertionError('Accepted changed checkpoint plan')
        bad=matrix[:,:5].copy();bad[:,4]=1.
        errors=runner.fit_case(('f'*64,1,bg,family,np.arange(n),bad,output,'fixture-plan'))
        assert len(errors)==5 and all(r['status']=='fit_error_requires_review' for r in errors)
        assert runner.fit_case(('f'*64,1,bg,family,np.arange(n),bad,output,'fixture-plan'))==errors
    result=dict(status='passed_polynomial_ml_worker_fixture',fits=15,direct_dense_candidates=checked,preserved_errors=5,successful_and_error_checkpoint_resume=True,changed_plan_rejected=True,maximum_dense_objective_error=maxerr,status_counts=fit_counts,script_sha256=sha(__file__),runner_sha256=sha('scripts/run_full_polynomial_ml.py'),scope='Synthetic worker and serialization validation only; no production fit acceptance or full-run ETA.')
    print(json.dumps(result,indent=2))
    with Path('metadata/polynomial_ml_runner_checks_20260927.json').open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')


if __name__=='__main__':main()
