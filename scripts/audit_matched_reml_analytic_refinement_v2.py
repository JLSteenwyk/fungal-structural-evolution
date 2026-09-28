#!/usr/bin/env python3
"""Wait for the full refinement, then independently read back every disposition."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
import psutil
from ancestral_chain_attempt import sha,write_json


def payload_digest(payload):
    return hashlib.sha256(json.dumps(payload,sort_keys=True,allow_nan=False).encode()).hexdigest()


def audit_payload(p,old):
    assert len(p['candidates'])==24
    candidates=p['candidates'];reference=candidates[0]
    assert reference['start']=='original_unoptimized_reference' and not reference['success']
    np.testing.assert_array_equal(reference['theta'],old['log1p_ratios'])
    np.testing.assert_allclose(reference['objective'],old['negative_profiled_reml'],rtol=1e-9,atol=1e-7)
    assert p['maximum_ratio']==old['maximum_ratio']
    upper=float(np.log1p(p['maximum_ratio']))
    faces=Counter(tuple(c['active']) for c in candidates[1:])
    assert len(faces)==8 and faces[(False,False,False)]==1 and faces[(True,True,True)]==4
    assert all(n==3 for face,n in faces.items() if any(face) and not all(face))
    for c in candidates:
        theta=np.asarray(c['theta']);assert theta.shape==(3,) and np.isfinite(theta).all() and np.all((theta>=0)&(theta<=upper))
        assert np.isfinite(c['objective'])
        if c['start']!='original_unoptimized_reference':assert np.all(theta[~np.asarray(c['active'],dtype=bool)]==0)
    best=min(candidates,key=lambda c:(c['objective'],not c['success']))
    np.testing.assert_array_equal(p['log1p_ratios'],best['theta'])
    np.testing.assert_allclose(p['negative_profiled_reml'],best['objective'],rtol=1e-9,atol=1e-7)
    theta=np.array(p['log1p_ratios']);gradient=np.array(p['analytic_gradient'])
    assert gradient.shape==(3,) and np.isfinite(gradient).all()
    projected=np.where(theta<=1e-7,np.minimum(gradient,0),np.where(theta>=upper-1e-7,np.maximum(gradient,0),gradient))
    np.testing.assert_array_equal(projected,p['projected_gradient'])
    np.testing.assert_allclose(np.expm1(theta),p['ratios'],rtol=1e-13,atol=1e-14)
    full=[c['objective'] for c in candidates if all(c['active']) and c['start']!='original_unoptimized_reference']
    checks=dict(best_optimizer_success=best['success'],projected_gradient_pass=bool(np.max(abs(projected))<=1e-3),
        upper_bound_contact=bool(np.any(theta>=upper-1e-6)),all_full_face_starts_agree=bool(np.ptp(full)<=1e-5),
        original_objective_not_worsened=bool(best['objective']<=reference['objective']))
    assert checks==p['checks']
    passed=all(value for key,value in checks.items() if key!='upper_bound_contact') and not checks['upper_bound_contact']
    assert p['status']==('refinement_passed_numerical_checks' if passed else 'refinement_requires_review')
    assert p['original_objective']==reference['objective'] and p['objective_improvement']==reference['objective']-best['objective']
    assert p['profiled_scale']>0 and p['residual_degrees_of_freedom']>0
    scales=np.asarray(p['covariate_scales']);assert (scales>0).all()
    conversion=np.r_[1.,1/scales]
    np.testing.assert_allclose(np.asarray(p['beta'])*conversion,p['raw_unit_beta'],rtol=1e-13,atol=1e-14)
    np.testing.assert_allclose(np.asarray(p['conditional_beta_covariance'])*conversion[:,None]*conversion[None,:],p['raw_unit_conditional_beta_covariance'],rtol=1e-13,atol=1e-14)
    return checks


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    config=json.loads(args.plan.read_text());ch=sha(args.plan)
    def verify():
        assert sha(args.plan)==ch
        for path,h in config['pins'].items():assert sha(path)==h,path
    verify();launch=json.loads(Path(config['producer_launch']).read_text())
    assert launch['plan_sha256']==sha(config['producer_plan'])
    while True:
        state=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','MainPID','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:break
        assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
        try:
            p=psutil.Process(launch['pid']);assert p.create_time()==launch['created'] and p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:time.sleep(1);continue
        time.sleep(30)
    assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
    verify();plan=json.loads(Path(config['producer_plan']).read_text());root=Path(plan['output'])
    receipt=json.loads((root/'receipt.json').read_text());assert receipt['plan_sha256']==sha(config['producer_plan'])
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h
    assert (root/'run_plan.json').read_bytes()==Path(config['producer_plan']).read_bytes()
    for path,h in plan['pins'].items():assert sha(path)==h,path
    production=json.loads(Path(plan['production_plan']).read_text())
    for path,h in production['pins'].items():assert sha(path)==h,path
    source=Path(plan['review_root'])/'remaining_review_cases.jsonl';cases=[json.loads(line) for line in source.read_text().splitlines()]
    expected={(r['fit_input_id'],r['tree']):r for r in cases};assert len(expected)==len(cases)==658
    seen=set();counts=Counter();remaining=[];rows=[]
    for line in (root/'manifest.jsonl').open():
        entry=json.loads(line);key=entry['fit_input_id'],entry['tree'];assert key in expected and key not in seen;seen.add(key)
        path=Path(entry['path']);assert path.resolve().is_relative_to(root.resolve());assert sha(path)==entry['sha256']
        r=json.loads(path.read_text());assert r['plan_sha256']==sha(config['producer_plan']) and r['source']==expected[key]
        p=r['payload'];assert r['payload_sha256']==payload_digest(p) and entry['status']==p['status']
        src=expected[key];assert sha(src['source_fit'])==src['source_fit_sha256']
        old=json.loads(Path(src['source_fit']).read_text())['payload'];counts[p['status']]+=1
        row=dict(fit_input_id=key[0],tree=key[1],status=p['status'],source_path=str(path),source_sha256=sha(path))
        if p['status']=='refinement_error_requires_review':remaining.append({**row,'reason':p['error_type']})
        else:
            checks=audit_payload(p,old)
            row.update(objective_improvement=p['objective_improvement'],maximum_projected_gradient=float(np.max(np.abs(np.asarray(p['projected_gradient'], dtype=float)))),
                original_raw_intercept=old['raw_unit_beta'][0],refined_raw_intercept=p['raw_unit_beta'][0],
                intercept_difference=p['raw_unit_beta'][0]-old['raw_unit_beta'][0],maximum_direct_candidate_objective_error=p['maximum_direct_candidate_objective_error'])
            if p['status']!='refinement_passed_numerical_checks':remaining.append({**row,'checks':checks})
        rows.append(row)
    assert seen==set(expected) and dict(counts)==receipt['status_counts'] and receipt['dispositions']==658
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    pd.DataFrame(rows).to_csv(out/'refinement_summary.tsv',sep='\t',index=False)
    write_json(out/'remaining_review.json',remaining);verify()
    write_json(out/'receipt.json',dict(status='complete_refinement_disposition_and_diagnostic_readback',plan_sha256=ch,
        source_receipt_sha256=sha(root/'receipt.json'),dispositions=len(rows),status_counts=dict(counts),remaining_review=len(remaining),
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='All 658 source-bound dispositions, 24-candidate enumeration, best selection, projected-gradient arithmetic, flags and coefficient units checked. Does not independently refit or recompute analytic derivatives from records; numerical fixtures and production direct-likelihood checks are separate. No calibrated confidence intervals, global-optimum proof or biological inference.'))


if __name__=='__main__':main()
