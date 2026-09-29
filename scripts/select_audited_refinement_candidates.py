#!/usr/bin/env python3
"""Build a separate numerical-candidate snapshot from audited fits/proposals."""
import json
import subprocess
from collections import Counter
from pathlib import Path
import numpy as np
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from screen_duplication_alignment_reuse import sha


def select_trial(proposal):
    trials=proposal['proposals']
    assert len(trials)==2 and all(t['passed'] for t in trials)
    assert {t['hessian_step'] for t in trials}=={1e-5,1e-6}
    assert max(t['objective'] for t in trials)-min(t['objective'] for t in trials)<=1e-5
    return min(trials,key=lambda t:(t['objective'],t['hessian_step']))


def main():
    pp=Path('metadata/selected_refinement_candidates_plan_20260929.json')
    plan=json.loads(pp.read_text());bindings={str(pp):sha(pp),**plan['pins']}
    for entry in plan['audits']:
        state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',entry['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
        proof=json.loads(Path(entry['proof']).read_text());root=Path(entry['root'])
        receipt=json.loads((root/'receipt.json').read_text())
        assert proof['status']==entry['status'] and proof['fits']==601 and proof['source_receipt_sha256']==sha(root/'receipt.json')
        bindings[str(root/'receipt.json')]=sha(root/'receipt.json')
        bindings[entry['proof']]=sha(entry['proof'])
        bindings.update(receipt['source_hashes'])
        bindings.update({str(root/n):h for n,h in receipt['artifacts'].items()})
    def verify():
        for path,digest in bindings.items():assert sha(path)==digest,path
    verify()
    refinements=Path(plan['refinements']);manifest=json.loads((refinements/'fit_manifest.json').read_text())
    expected={(r['fit_input_id'],r['tree']) for r in manifest};assert len(expected)==len(manifest)==601
    proposals=[]
    for root in plan['proposals']:
        rows={}
        for line in (Path(root)/'dispositions.jsonl').open():
            row=json.loads(line);key=row['fit_input_id'],row['tree'];assert key not in rows;rows[key]=row
        assert set(rows)==expected;proposals.append(rows)
    factors={}
    for path in Path(plan['factors']).glob('*.npz'):
        with np.load(path) as a:factors[path.stem]=a['factor']
    out=Path(plan['output']);out.mkdir(exist_ok=False)
    counts=Counter();selected=[]
    for item in manifest:
        key=item['fit_input_id'],item['tree'];source=refinements/item['path']
        assert sha(source)==item['sha256'];saved=json.loads(source.read_text());p=saved['payload']
        proposal_index=None;step=None
        if p['status']=='ml_refinement_passed_numerical_checks':
            assert all(rows[key]['disposition']=='already_passed_no_proposal' for rows in proposals)
            theta=np.asarray(p['log1p_ratios']);selected_objective=p['negative_profiled_ml'];kind='audited_refinement'
        else:
            assert p['checks']==dict(best_optimizer_success=True,projected_gradient_pass=False,upper_bound_contact=False,all_full_face_starts_agree=True,original_objective_not_worsened=True)
            choices=[(i,rows[key]) for i,rows in enumerate(proposals) if rows[key]['disposition'] in ['two_local_proposals_passed','two_lower_face_proposals_passed']]
            assert len(choices)==1
            proposal_index,row=choices[0];trial=select_trial(row['proposal'])
            theta=np.asarray(trial['theta']);selected_objective=trial['objective'];step=trial['hessian_step']
            kind='interior_local_proposal' if proposal_index==0 else 'lower_face_local_proposal'
        input_path=Path(plan['inputs'])/(key[0]+'.npz');assert sha(input_path)==saved['input_sha256']
        with np.load(input_path) as a:
            values,bg,family,indices=a['matrix'],a['background'],a['family'],a['pattern_rows']
        active=np.ptp(values[:,1:],axis=0)>1e-12;scales=values[:,1:][:,active].std(axis=0)
        np.testing.assert_array_equal(active,p['active_covariates']);np.testing.assert_array_equal(scales,p['covariate_scales'])
        x=np.column_stack([np.ones(len(values)),values[:,1:][:,active]/scales]);y=values[:,0];factor=factors[key[1]][indices]
        ratios=np.expm1(theta);upper=np.log1p(p['maximum_ratio'])
        assert np.isfinite(theta).all() and np.all(theta>=0) and np.all(theta<upper-1e-6)
        fitted=profiled_ml(MatchedCovariance(bg,family,factor,1.,*ratios),x,y)
        np.testing.assert_allclose(fitted['negative_profiled_ml'],selected_objective,rtol=1e-10,atol=1e-8)
        baseline=profiled_ml(MatchedCovariance(bg,family,factor,1.,*np.expm1(p['log1p_ratios'])),x,y)
        assert fitted['negative_profiled_ml']<=baseline['negative_profiled_ml']
        grad=evaluate_gradient(CachedMatchedML(bg,family,factor,x,y),ratios)['log1p_ratio_gradient']
        projected=np.where(theta<=1e-7,np.minimum(grad,0),grad)
        assert max(abs(projected))<=1e-3
        if kind=='lower_face_local_proposal':
            fixed=np.asarray(p['log1p_ratios'])==0
            assert np.all(theta[fixed]==0) and np.all(grad[fixed]>=0)
        conversion=np.r_[1.,1./scales]
        result=dict(fit_input_id=key[0],tree=key[1],polynomial_degree=saved['polynomial_degree'],source_refinement=str(source),source_refinement_sha256=item['sha256'],original_production_fit=saved['source_fit_path'],original_production_sha256=saved['source_fit_sha256'],input_sha256=saved['input_sha256'],selection_kind=kind,proposal_index=proposal_index,hessian_step=step,source_refinement_status=p['status'],source_refinement_checks=p['checks'],records=len(y),active_covariates=active.tolist(),covariate_scales=scales.tolist(),log1p_ratios=theta.tolist(),ratios=ratios.tolist(),maximum_ratio=p['maximum_ratio'],gradient=grad.tolist(),projected_gradient=projected.tolist(),objective_change_vs_refinement=fitted['negative_profiled_ml']-baseline['negative_profiled_ml'],fitted={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in fitted.items()},raw_unit_beta=(fitted['beta']*conversion).tolist(),raw_unit_conditional_beta_covariance=(fitted['conditional_beta_covariance']*conversion[:,None]*conversion[None,:]).tolist(),status='selected_numerical_candidate_pending_readback')
        path=out/(key[0]+'__'+key[1]+'.json')
        with path.open('x') as h:json.dump(result,h,indent=2,allow_nan=False);h.write('\n')
        assert json.loads(path.read_text())==result
        selected.append(dict(fit_input_id=key[0],tree=key[1],path=path.name,sha256=sha(path),selection_kind=kind));counts[kind]+=1
        if len(selected)%50==0:print('Selected/recomputed frozen candidates',len(selected),'/601',flush=True)
    verify()
    (out/'manifest.json').write_text(json.dumps(selected,indent=2)+'\n')
    receipt=dict(status='complete_selected_refinement_candidates_pending_readback',fits=len(selected),selection_counts=dict(counts),source_hashes=bindings,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Separate frozen numerical-candidate snapshot. Lowest direct proposal objective wins; smaller Hessian step breaks exact ties only. Original fits and flags retained. No production integration, full-grid acceptance, calibrated intervals or biological inference.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_hashes','artifacts']}),flush=True)


if __name__=='__main__':main()
