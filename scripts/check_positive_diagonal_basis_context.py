#!/usr/bin/env python3
"""Qualify reused positive-diagonal products against complete separate audits."""
import argparse
from datetime import datetime, timezone
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha
from check_retained_shared_entity_candidates import inputs
from check_positive_diagonal_kernel_products import dense, rejected
from covariance_basis_context import ComponentKernelProducts
from independent_positive_diagonal_kernel_products import PositiveDiagonalKernelProducts
from positive_diagonal_basis_context import PositiveDiagonalBasisContext


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(exist_ok=False);assert not a.receipt.exists()
    exports=[];checks=reordered=altered_rejections=immutable=0;maximum_raw=maximum_projected=0.
    rng=np.random.default_rng(202610036)
    for exception,mode in itertools.product([False,True],['signed','unsigned']):
        src,rows,x,y,ops,audit,old,parent=inputs(exception,mode);n=len(rows)
        operators={'target_node':sparse.eye(n,format='csr'),'background_node':ops['background_node']}
        if exception:operators['model_pair']=ops['model_pair']
        operators['family_intercept']=ops['family_intercept']
        factor=src['factors'][audit['tree']].copy()
        context=PositiveDiagonalBasisContext(src['labels'],operators,factor).design(x)
        patterns=[np.ones(n),np.full(n,2.),np.where(np.arange(n)%2,.5,2.),np.linspace(.125,3.,n),
            1.+np.arange(n)*2.**-40,np.geomspace(1e-4,1e4,n)]
        first={}
        for j,d in enumerate(patterns):
            actual=context.audit(d)
            separate=ComponentKernelProducts(src['labels'],operators,d).tree(factor).audit(x)
            independent=PositiveDiagonalKernelProducts(src['labels'],operators,d).tree(factor).project(x)
            reference=dense(operators,factor,d,x)
            for name in ['raw_gram','projected_gram','raw_roundoff_envelope','projected_roundoff_envelope']:
                np.testing.assert_allclose(actual[name],separate[name],rtol=3e-9,atol=2e-8)
            for field,kind,expected in [('raw_gram','raw',reference[0]),('projected_gram','projected',reference[1])]:
                matrix=np.asarray(actual[field]);np.testing.assert_allclose(matrix,expected,rtol=3e-9,atol=2e-8)
                np.testing.assert_allclose(matrix,independent[kind],rtol=3e-9,atol=2e-8)
                if kind=='raw':maximum_raw=max(maximum_raw,float(np.max(abs(matrix-expected))))
                else:maximum_projected=max(maximum_projected,float(np.max(abs(matrix-expected))))
            for field in ['raw_diagnostics','reml_diagnostics']:
                for key in ['disposition','rank','rank_at_one_tenth_tolerance','rank_at_ten_times_tolerance','unresolved_kernel_names']:
                    assert actual[field][key]==separate[field][key],(j,field,key)
            if j in [0,1]:assert actual['raw_diagnostics']['disposition']=='dependent_covariance_bases_require_review'
            first[j]=actual;checks+=1
        # Changing call order never changes the cached inputs or results.
        for j in rng.permutation(len(patterns)):
            assert context.audit(patterns[j])==first[j];reordered+=1
        # The source factor is copied before read-only caching.
        factor[:]=np.nan
        assert context.audit(patterns[2])==first[2];immutable+=1
        # Returned nested lists are independent, not views into shared caches.
        returned=context.audit(patterns[2]);returned['kernel_names'][0]='foreign';returned['raw_gram'][0][0]=-1
        assert context.audit(patterns[2])==first[2];immutable+=1
        for d in [np.zeros(n),-np.ones(n),np.full(n,np.nan),np.full(n,np.inf),np.ones(n+1)]:
            rejected(lambda:context.audit(d));altered_rejections+=1
        for design in [x*0,np.column_stack([x,x[:,0]]),x[:-1],np.ones((n,n)),np.full_like(x,np.nan)]:
            rejected(lambda:context.parent.design(design));altered_rejections+=1
        exports.append(dict(pair_exception=exception,loading_mode=mode,cases=list(first.values()),scientific_eligibility=False))
    assert (checks,reordered,immutable,altered_rejections)==(24,24,8,40)
    artifact=a.output/'reused_positive_diagonal_gram_cases.json';artifact.write_text(json.dumps(exports,indent=2,allow_nan=False)+'\n')
    modules=['positive_diagonal_basis_context','check_positive_diagonal_basis_context',
        'independent_positive_diagonal_kernel_products','check_positive_diagonal_kernel_products',
        'covariance_basis_context','covariance_basis_audit','check_retained_shared_entity_candidates',
        'run_positive_diagonal_kernel_software_stage']
    result=dict(status='passed_reused_fresh_positive_diagonal_raw_reml_gram_contracts_v1',checked_utc=datetime.now(timezone.utc).isoformat(),
        complete_separate_diagonal_audits=checks,reordered_diagonal_replays=reordered,
        copied_source_and_detached_result_checks=immutable,invalid_diagonal_or_design_cases_rejected=altered_rejections,
        maximum_dense_raw_error=maximum_raw,maximum_dense_projected_error=maximum_projected,
        unchanged_comparison_rtol=3e-9,unchanged_comparison_atol=2e-8,
        source_hashes={'scripts/'+m+'.py':sha('scripts/'+m+'.py') for m in modules},artifacts={str(artifact):sha(artifact)},
        raw_reml_basis_qualification_complete=False,scientific_eligibility=False,
        scope='All24declared q5/q6 mode/diagonal cases compared against separate complete component, independent latent and dense audits, with retained constant-diagonal dependency reviews. Fresh constant-kernel reuse does not inherit saved uniform audits/envelopes or alter any model. Complete input/result isolation checked. No measured whole-data speedup, real weighted source qualification, variance fitting or scientific acceptance.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts','scope']}),flush=True)


if __name__=='__main__':main()
