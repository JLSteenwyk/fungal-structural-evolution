"""Synthetic end-to-end simulation/refit plumbing, not fitted-data calibration."""
from pathlib import Path
from unittest.mock import patch
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from refit_matched_simulation import refit_one,coverage_accounting,response_seed


def main():
    threadpool_limits(1)
    rng=np.random.default_rng(314159);n=48
    bg=np.repeat(np.arange(16),3);fam=bg//4;factor=rng.normal(size=(n,4))/2
    x=np.column_stack([np.ones(n),rng.normal(size=n)]);beta=np.array([.1,-.4])
    out=Path('results/model_validation/matched-simulation-refit-checks-20260928-v1');out.mkdir(parents=True,exist_ok=False)
    records=[]
    for i,ratios in enumerate([[0,0,0],[.4,.8,.2]]):
        batch=[]
        for replicate in range(3):
            r=refit_one(bg,fam,factor,x,beta,1.3,ratios,42,f'synthetic-{i}',replicate)
            assert r['status']!='refit_error_requires_review',r
            assert len(r['fit']['candidates'])==24
            write_json(out/f'case_{i}_replicate_{replicate}.json',r);batch.append(r);records.append(r)
        accounting=coverage_accounting(batch,2)
        assert accounting['attempted']==3 and accounting['qualified']+accounting['unresolved']==3
    # Failure accounting must retain failed attempts in the denominator.
    with patch('refit_matched_simulation.refine',side_effect=RuntimeError('fixture failure')):
        failed=refit_one(bg,fam,factor,x,beta,1.3,[0,0,0],42,'synthetic-0',99)
    assert failed['status']=='refit_error_requires_review'
    mixed=coverage_accounting(records[:3]+[failed],2)
    assert mixed['attempted']==4 and mixed['unresolved']>=1
    assert response_seed(42,'synthetic-0',0)==response_seed(42,'synthetic-0',0)
    assert response_seed(42,'synthetic-0',0)!=response_seed(42,'synthetic-0',1)
    assert response_seed(42,'synthetic-0',0)!=response_seed(42,'synthetic-1',0)
    try:coverage_accounting([records[0],records[0]],2)
    except ValueError:pass
    else:raise AssertionError('Duplicate replicates accepted')
    result=dict(status='passed_synthetic_end_to_end_refit_checks',simulated_refits=6,
        candidates_checked=144,qualified=sum(r['status']=='refit_numerically_checked' for r in records),
        review=sum(r['status']=='refit_requires_review' for r in records),failure_accounting=mixed,
        artifacts={p.name:sha(p) for p in out.iterdir()},
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/refit_matched_simulation.py'),Path('scripts/simulate_matched_working_model.py'),Path('scripts/refine_matched_reml_analytic.py')]},
        scope='Six synthetic implementation fixtures, not a biological pilot or coverage estimate. Full-grid production remains unlaunched.')
    write_json(out/'receipt.json',result)
    write_json(Path('metadata/matched_simulation_refit_checks_20260928.json'),result)
    print('checked',result['simulated_refits'],'refits;',result['qualified'],'qualified;',result['review'],'review')


if __name__=='__main__':main()
