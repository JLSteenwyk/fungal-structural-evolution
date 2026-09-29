#!/usr/bin/env python3
"""Replay serialized whole-protein refinements, retaining unresolved flags."""
import argparse
import itertools
import json
import subprocess
import time
from collections import Counter
from pathlib import Path
import numpy as np
import psutil
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_ml_gradient import evaluate_gradient
from matched_mixed_covariance import MatchedCovariance
from screen_duplication_alignment_reuse import sha


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args()
    config=json.loads(args.plan.read_text())
    launch=json.loads(Path(config['launch']).read_text())
    bindings={str(args.plan):sha(args.plan),**config['pins']}
    def verify():
        for path,h in bindings.items():assert sha(path)==h,path
    verify()
    while psutil.pid_exists(launch['pid']):
        try:
            process=psutil.Process(launch['pid'])
            if process.create_time()!=launch['created'] or process.status()==psutil.STATUS_ZOMBIE:break
            assert process.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    verify()
    assert sha(launch['plan'])==launch['plan_sha256']
    plan=json.loads(Path(launch['plan']).read_text());root=Path(plan['output'])
    receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_frozen_whole_protein_refinement_pending_serialized_review'
    bindings.update(receipt['source_hashes']);bindings[str(root/'receipt.json')]=sha(root/'receipt.json')
    verify()
    proof=json.loads(Path(plan['original_replay']).read_text())
    expected={(r['fit_input_id'],r['tree']) for r in proof['results']}
    original_rows={}
    for line in Path(plan['frozen_manifest']).read_text().splitlines():
        r=json.loads(line)
        if (r['fit_input_id'],r['tree']) in expected:original_rows[r['fit_input_id'],r['tree']]=r
    assert set(original_rows)==expected
    inputs=Path(plan['inputs']);recipes={}
    for line in (inputs/'input_manifest.jsonl').open():
        row=json.loads(line)
        if row['fit_input_id'] in {k[0] for k in expected}:recipes[row['fit_input_id']]=row
    seen=set();counts=Counter();candidate_count=0;maximum_error=0.
    for item in receipt['dispositions']:
        key=item['fit_input_id'],item['tree'];assert key in expected and key not in seen;seen.add(key)
        path=root/item['path'];assert sha(path)==item['sha256'];bindings[str(path)]=item['sha256']
        saved=json.loads(path.read_text());source=original_rows[key]
        assert (saved['fit_input_id'],saved['tree'])==key
        assert saved['source_fit']==source['path'] and saved['source_fit_sha256']==source['sha256']==sha(source['path'])
        original=json.loads(Path(source['path']).read_text());old=original['payload']
        recipe=recipes[key[0]];ip=inputs/recipe['path'];assert sha(ip)==recipe['sha256']==saved['input_sha256']==original['input_sha256']
        assert saved['columns']==original['specification']['columns']==recipe['recipe']['specification']['columns']
        with np.load(ip,allow_pickle=False) as a:matrix,bg,family,index=a['matrix'],a['background'],a['family'],a['pattern_rows']
        fp=Path(plan['factors'])/(key[1]+'.npz');assert sha(fp)==bindings[str(fp)]
        with np.load(fp,allow_pickle=False) as a:factor=a['factor'][index]
        assert np.all(matrix[:,1]==1)
        scales=np.std(matrix[:,2:],axis=0);np.testing.assert_array_equal(saved['covariate_scales'],scales)
        x=np.column_stack([matrix[:,1],matrix[:,2:]/scales]);y=matrix[:,0];p=saved['payload']
        candidates = p['candidates']
        grid = [([True] * 3, 'original_unoptimized_reference')]
        for face in itertools.product([False, True], repeat=3):
            for start in (['0.05', '1.0', '20.0'] if any(face) else ['exact_zero']):
                grid.append((list(face), start))
            if all(face):
                grid.append((list(face), 'original_parameters'))
        assert [(c['active'], c['start']) for c in candidates] == grid and len(candidates) == 24
        np.testing.assert_array_equal(candidates[0]['theta'], old['log1p_ratios'])
        assert not candidates[0]['success'] and p['maximum_ratio'] == old['maximum_ratio']
        upper = np.log1p(p['maximum_ratio'])
        for c in candidates:
            theta = np.asarray(c['theta'])
            assert theta.shape == (3,) and np.isfinite(theta).all() and (theta >= 0).all() and (theta <= upper).all()
            assert (theta[~np.asarray(c['active'])] == 0).all()
            direct = profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(theta)), x, y)
            error = abs(direct['negative_profiled_ml'] - c['objective'])
            np.testing.assert_allclose(direct['negative_profiled_ml'], c['objective'], rtol=1e-9, atol=1e-7)
            maximum_error = max(maximum_error, error)
            candidate_count += 1
        best = min(candidates, key=lambda c: (c['objective'], not c['success']))
        theta = np.asarray(best['theta'])
        np.testing.assert_array_equal(p['log1p_ratios'], theta)
        np.testing.assert_allclose(p['ratios'], np.expm1(theta), rtol=1e-12, atol=1e-12)
        grad = evaluate_gradient(CachedMatchedML(bg, family, factor, x, y), np.expm1(theta))['log1p_ratio_gradient']
        np.testing.assert_allclose(p['analytic_gradient'], grad, rtol=1e-9, atol=1e-9)
        projected = np.where(theta <= 1e-7, np.minimum(grad, 0), np.where(theta >= upper - 1e-7, np.maximum(grad, 0), grad))
        np.testing.assert_allclose(p['projected_gradient'], projected, rtol=1e-9, atol=1e-9)
        full = [c['objective'] for c in candidates if all(c['active']) and c['start'] != 'original_unoptimized_reference']
        checks = dict(best_optimizer_success=best['success'], projected_gradient_pass=bool(max(abs(projected)) <= 1e-3), upper_bound_contact=bool(any(theta >= upper - 1e-6)), all_full_face_starts_agree=bool(max(full) - min(full) <= 1e-5), original_objective_not_worsened=bool(best['objective'] <= candidates[0]['objective']))
        assert p['checks'] == checks
        passed = all(v for k, v in checks.items() if k != 'upper_bound_contact') and not checks['upper_bound_contact']
        assert p['status'] == item['status'] == ('ml_refinement_passed_numerical_checks' if passed else 'ml_refinement_requires_review')
        fitted = profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(theta)), x, y)
        for name in ['negative_profiled_ml', 'beta', 'profiled_scale', 'residual_quadratic', 'conditional_beta_covariance']:
            np.testing.assert_allclose(p[name], fitted[name], rtol=1e-8, atol=1e-8)
        conversion = np.r_[1., 1. / scales]
        np.testing.assert_allclose(p['raw_unit_beta'], fitted['beta'] * conversion, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(p['raw_unit_conditional_beta_covariance'], fitted['conditional_beta_covariance'] * conversion[:, None] * conversion[None, :], rtol=1e-10, atol=1e-10)
        assert p['variance_profile_denominator'] == len(y) and p['likelihood'] == 'ordinary_gaussian_ml'
        assert p['original_objective'] == candidates[0]['objective'] and p['objective_improvement'] == candidates[0]['objective'] - best['objective']
        counts[p['status']] += 1
        assert item['checks']==p['checks'] and item['objective_improvement']==p['objective_improvement']
        assert saved['candidate_likelihoods_replayed']==24
        print('Replayed whole-protein refinement',len(seen),'/',len(expected),flush=True)
    assert seen==expected and len(seen)==receipt['fits']==proof['flagged_dispositions']
    assert candidate_count==receipt['candidate_likelihoods_replayed']==24*len(seen)
    verify()
    result=dict(status='passed_serialized_frozen_whole_protein_refinement_readback',fits=len(seen),counts=dict(counts),candidate_likelihoods_replayed=candidate_count,maximum_objective_error=maximum_error,source_receipt_sha256=sha(root/'receipt.json'),source_hashes=bindings,producer_terminal_state=state,scope='All candidates, coefficients, covariance transformations, gradients, numerical decisions and source identities replayed from saved files in a separate process. Mathematical libraries are shared. Unresolved statuses remain; no production overlay, global-optimum proof, later-flag resolution or calibrated inference.')
    with Path(config['output']).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)


if __name__=='__main__':main()
