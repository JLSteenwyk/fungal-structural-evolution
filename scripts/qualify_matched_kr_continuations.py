"""Separate numerical qualification overlay; selected parameters are unchanged."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import t
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from calibrate_matched_kr import evaluate_refit,summarize
from matched_calibration_intervals import summarize_replicates
from simulate_matched_working_model import simulate
from matched_mixed_covariance import MatchedCovariance,profiled_reml
from matched_reml_gradient import evaluate_gradient
from cached_matched_likelihood import CachedMatchedLikelihood


def qualification(fit, retries):
    assert not fit['checks']['all_full_face_starts_agree']
    checks={k:v for k,v in fit['checks'].items() if k!='all_full_face_starts_agree'}
    full=[c for c in fit['candidates'] if all(c['active']) and c['start']!='original_unoptimized_reference']
    assert len(full)==4 and set(retries).issubset({c['start'] for c in full})
    objectives=[];stationary=[];success=[];strict=[]
    for candidate in full:
        retry=retries.get(candidate['start'])
        if retry is None:objectives.append(candidate['objective']);continue
        assert retry['original_candidate']==candidate
        end=retry['attempts'][-1];objectives.append(end['objective'])
        stationary.append(end['projected_gradient_maximum']<=1e-3)
        strict.append(end['projected_gradient_maximum']<=1e-6)
        success.append(end['optimizer_success'])
    checks.update(completed_full_face_starts_agree=bool(np.ptp(objectives)<=1e-5),
                  completed_endpoints_agree_with_selected=bool(np.max(abs(np.array(objectives)-fit['negative_profiled_reml']))<=1e-5),
                  continued_endpoints_stationary=bool(stationary and all(stationary)),
                  continued_endpoints_optimizer_success=bool(success and all(success)))
    qualified=all(v for k,v in checks.items() if k!='upper_bound_contact') and not checks['upper_bound_contact']
    return dict(qualified=bool(qualified),checks=checks,full_face_objective_spread=float(np.ptp(objectives)),
                stricter_retry_threshold_all_met=bool(strict and all(strict)))


def main():
    threadpool_limits(1)
    root=Path('results/model_validation/matched-kr-simulations-20260928-v1')
    retryroot=Path('results/model_validation/matched-kr-start-retries-20260928-v1')
    plan=json.loads((root/'run_plan.json').read_text())
    original=json.loads((root/'receipt.json').read_text())
    retryreceipt=json.loads((retryroot/'receipt.json').read_text())
    audit=Path('results/model_validation/matched-kr-simulation-audit-20260928-v1/receipt.json')
    assert json.loads(audit.read_text())['source_receipt_sha256']==sha(root/'receipt.json')
    for name,digest in original['artifacts'].items():assert sha(root/name)==digest
    for name,digest in plan['pins'].items():assert sha(name)==digest
    retries={};retry_paths={}
    for row in retryreceipt['rows']:
        path=retryroot/row['path'];assert sha(path)==row['sha256'];saved=json.loads(path.read_text())
        key=(saved['case'],saved['replicate']);start=saved['start']
        assert start not in retries.setdefault(key,{})
        retries[key][start]=saved;retry_paths.setdefault(key,{})[str(path)]=row['sha256']
    assert len(retries)==184 and sum(map(len,retries.values()))==185
    out=Path('results/model_validation/matched-kr-continuation-overlay-20260928-v1');out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((root/'manifest.json').read_text());assert len(manifest)==15984
    expected={(c['id'],r) for c in plan['cases'] for r in range(plan['replicates'])}
    assert {(r['case'],r['replicate']) for r in manifest}==expected
    kr={};baseline={};overlays=[];unchanged=0;qualified=0
    for case in plan['cases']:
        with np.load(case['design'],allow_pickle=False) as h:a={k:h[k] for k in h.files}
        interval_records=[];baseline_records=[]
        for entry in [r for r in manifest if r['case']==case['id']]:
            assert sha(entry['path'])==entry['sha256']
            saved=json.loads(Path(entry['path']).read_text());r=saved['refit'];key=(case['id'],r['replicate'])
            assert saved['case']==case and saved['plan_sha256']==sha(root/'run_plan.json')
            if key not in retries:
                assert r['status']=='refit_numerically_checked'
                interval_records.append(saved['interval']);baseline_records.append(r);unchanged+=1;continue
            assert r['status']=='refit_requires_review'
            for retry in retries[key].values():assert retry['source_sha256']==entry['sha256'] and retry['source_path']==entry['path']
            q=qualification(r['fit'],retries[key]);assert q['qualified'],q
            y=simulate(a['background'],a['family'],a['factor'],a['design'],a['beta'],case['scale'],case['ratios'],1,
                       np.random.default_rng(np.random.SeedSequence(r['seed_entropy'])))[:,0]
            assert hashlib.sha256(np.asarray(y,dtype='<f8').tobytes()).hexdigest()==r['response_sha256']
            f=r['fit'];check=profiled_reml(MatchedCovariance(a['background'],a['family'],a['factor'],1.,*f['ratios']),a['design'],y)
            for name in ['beta','profiled_scale','conditional_beta_covariance','negative_profiled_reml']:
                np.testing.assert_allclose(check[name],f[name],rtol=1e-7,atol=1e-8)
            cache=CachedMatchedLikelihood(a['background'],a['family'],a['factor'],a['design'],y)
            g=evaluate_gradient(cache,f['ratios'])['log1p_ratio_gradient'];theta=np.array(f['log1p_ratios']);upper=np.log1p(f['maximum_ratio'])
            pg=np.where(theta<=1e-7,np.minimum(g,0),np.where(theta>=upper-1e-7,np.maximum(g,0),g))
            assert np.max(abs(pg))<=1e-3
            # Adapter is local only; original serialized fit/status remain unchanged.
            effective=dict(r,status='refit_numerically_checked')
            interval=evaluate_refit(effective,a['background'],a['family'],a['factor'],a['design'],a['beta'])
            assert interval['status']=='interval_dispositions_recorded'
            z=(np.array(f['beta'])-a['beta'])/np.sqrt(np.diag(f['conditional_beta_covariance']))
            effective['nominal_095_t_interval_covers_truth']=(abs(z)<=t.ppf(.975,f['residual_degrees_of_freedom'])).tolist()
            record=dict(case=case['id'],replicate=r['replicate'],source_path=entry['path'],source_sha256=entry['sha256'],
                        original_status=r['status'],qualification=q,retry_bindings=retry_paths[key],
                        interval=interval,conditional_t_covers_truth=effective['nominal_095_t_interval_covers_truth'],
                        selected_parameters_unchanged=True,selected_projected_gradient_maximum=float(np.max(abs(pg))))
            path=out/f"{case['id']}-{r['replicate']:05d}.json";write_json(path,record)
            overlays.append(dict(case=case['id'],replicate=r['replicate'],path=path.name,sha256=sha(path)))
            interval_records.append(interval);baseline_records.append(effective);qualified+=1
        assert len(interval_records)==len(baseline_records)==999
        kr[case['id']]=summarize(interval_records)
        baseline[case['id']]=summarize_replicates(baseline_records,case['coefficients'])
    assert qualified==184 and unchanged==15800
    write_json(out/'overlay_manifest.json',overlays);write_json(out/'kr_coverage.json',kr);write_json(out/'conditional_t_coverage.json',baseline)
    result=dict(status='completed_continuation_qualification_overlay_pending_independent_readback',
        total_attempts=15984,unchanged_original_qualified=unchanged,newly_qualified_by_continuation=qualified,
        source_receipt_sha256=sha(root/'receipt.json'),source_audit_receipt_sha256=sha(audit),
        retry_receipt_sha256=sha(retryroot/'receipt.json'),script_sha256=sha(__file__),
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='Same frozen responses, selected parameters and15984-denominator. Qualification '
              'uses original1e-3 gradient/1e-5 start-agreement criteria plus successful endpoint '
              'continuations. Stricter retry diagnostic preserved. Numerical protocol update '
              'after initial study, not new replicates or biological/coverage certification.')
    write_json(out/'receipt.json',result);write_json(Path('metadata/matched_kr_continuation_overlay_20260928.json'),result)
    print({k:v for k,v in result.items() if k!='artifacts'})


if __name__=='__main__':main()
