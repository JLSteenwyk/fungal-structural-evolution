#!/usr/bin/env python3
"""Replay lower-face helper against separately checked real-case proposals."""
import json
from pathlib import Path
import numpy as np
from lower_face_gradient_polish import propose
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from screen_duplication_alignment_reuse import sha


def main():
    pp=Path('metadata/boundary_refinement_diagnostic_plan_20260929.json')
    plan=json.loads(pp.read_text())
    reference_path=Path('metadata/boundary_refinement_proposals_checked_20260929.json')
    reference=json.loads(reference_path.read_text())
    bindings={**plan['pins'],str(pp):sha(pp),str(reference_path):sha(reference_path),'scripts/lower_face_gradient_polish.py':sha('scripts/lower_face_gradient_polish.py')}
    for path,digest in bindings.items():assert sha(path)==digest,path
    original=json.loads(Path(plan['fit']).read_text())['payload']
    with np.load(plan['inputs']) as a:
        values,bg,family,rows=a['matrix'],a['background'],a['family'],a['pattern_rows']
    with np.load(plan['factor']) as a:factor=a['factor'][rows]
    active=np.ptp(values[:,1:],axis=0)>1e-12
    cov=values[:,1:][:,active]
    x=np.column_stack([np.ones(len(values)),cov/cov.std(axis=0)]);y=values[:,0]
    cache=CachedMatchedML(bg,family,factor,x,y)
    def objective(t):return profiled_ml(MatchedCovariance(bg,family,factor,1.,*np.expm1(t)),x,y)['negative_profiled_ml']
    def gradient(t):return evaluate_gradient(cache,np.expm1(t))['log1p_ratio_gradient']
    result=propose(objective,gradient,original['log1p_ratios'],np.log1p(original['maximum_ratio']))
    assert result['status']=='two_lower_face_proposals_passed'
    for trial,checked in zip(result['proposals'],reference['checks']):
        assert trial['hessian_step']==checked['hessian_step']
        np.testing.assert_allclose(trial['theta'],checked['theta'],rtol=0,atol=1e-14)
        np.testing.assert_allclose(trial['objective'],checked['objective'],rtol=0,atol=1e-9)
        for fd,old in zip(trial['finite_differences'],checked['checks']):
            assert fd['step']==old['step']
            np.testing.assert_allclose(fd['gradient'],old['direct_gradient'],rtol=1e-9,atol=1e-9)
    for path,digest in bindings.items():assert sha(path)==digest,path
    proof=dict(status='passed_lower_face_helper_real_case',source_hashes=bindings,checker_sha256=sha(__file__),result=result,scope='Helper reproduces separately checked proposals for one real case. No production fit replacement or full-grid acceptance.')
    with Path('metadata/lower_face_polish_real_case_20260929.json').open('x') as h:json.dump(proof,h,indent=2);h.write('\n')
    print(proof['status'])


if __name__=='__main__':main()
