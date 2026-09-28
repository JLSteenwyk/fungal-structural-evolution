"""Check selected-estimate mapping, original units, omissions and review flags."""
import copy
import json
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from evaluate_selected_matched_intervals import evaluate,COVARIATES
from integrate_matched_refinements import overlay
from matched_kr_scalar import scalar_interval
from matched_mixed_covariance import MatchedCovariance,profiled_reml


def main():
    threadpool_limits(1)
    root=Path('results/structural_comparisons/full-working-model-grid-export-20260927-v1')
    receipt=json.loads((root/'receipt.json').read_text())
    assert sha(root/'unique_fits.parquet')==receipt['artifacts']['unique_fits.parquet']
    frame,_=overlay(pd.read_parquet(root/'unique_fits.parquet'),[])
    frame=frame.set_index(['fit_input_id','tree'],drop=False)
    reference=Path('results/model_validation/matched-information-timing-20260928-v1')
    r=json.loads((reference/'receipt.json').read_text())
    for name,digest in r['artifacts'].items():assert sha(reference/name)==digest
    results=[];sources={}
    for job in [j for j in r['jobs'] if j['records']==148]:
        spec=job['job'];row=frame.loc[(spec['fit_input_id'],spec['tree'])].to_dict()
        with np.load(spec['cache'],allow_pickle=False) as h:a={k:h[k] for k in h.files}
        with np.load(spec['factor'],allow_pickle=False) as h:f=h['factor']
        old=json.loads(Path(spec['original_fit']).read_text())['payload']
        result=evaluate(row,a,f)
        with np.load(reference/job['artifact'],allow_pickle=False) as h:c={k:h[k] for k in h.files}
        conversion=np.r_[1.,1/a['covariate_scales']]
        names=['intercept']+[n for n,keep in zip(COVARIATES,a['active_covariates']) if keep]
        for i,name in enumerate(names):
            interval=scalar_interval(c,old['beta'],np.eye(len(names))[i])
            checked=result['intervals'][name]
            for key in ['estimate','lower','upper']:
                np.testing.assert_allclose(checked[key],interval[key]*conversion[i],rtol=1e-8,atol=1e-9)
            np.testing.assert_allclose(checked['denominator_df'],interval['denominator_df'],rtol=1e-8)
        results.append(result)
        for key in ['cache','factor','original_fit']:sources[spec[key]]=sha(spec[key])
    assert len(results)==5
    flagged=copy.deepcopy(row);flagged['selected_review_required']=True
    assert evaluate(flagged,a,f)['status']=='selected_fit_requires_review'
    bad=copy.deepcopy(row);bad['selected_intercept']+=1
    try:evaluate(bad,a,f)
    except AssertionError:pass
    else:raise AssertionError('Changed selected coefficient accepted')
    # A separate fixture has two omitted, identically zero covariates.
    rng=np.random.default_rng(20260928);n=48;bg=np.repeat(np.arange(16),3);fam=bg//4
    f=rng.normal(size=(n,4))/2;covariates=np.column_stack([rng.normal(size=(n,2)),np.zeros((n,2))])
    scales=np.std(covariates[:,:2],axis=0);x=np.column_stack([np.ones(n),covariates[:,:2]/scales]);y=rng.normal(size=n)
    fit=profiled_reml(MatchedCovariance(bg,fam,f,1.,.3,.2,.4),x,y)
    raw=fit['beta']*np.r_[1.,1/scales]
    synthetic=dict(fit_input_id='synthetic',tree='fixture',selected_source_sha256='fixture',selection='fixture',
        selected_review_required=False,records=n,selected_intercept=raw[0],
        selected_profiled_scale=fit['profiled_scale'],selected_negative_profiled_reml=fit['negative_profiled_reml'],
        selected_conditional_intercept_variance=fit['conditional_beta_covariance'][0,0],
        selected_variance_ratio_background=.3,selected_variance_ratio_family_component=.2,selected_variance_ratio_species=.4)
    for i,name in enumerate(COVARIATES):synthetic['selected_coefficient_'+name]=raw[i+1] if i<2 else None
    arrays=dict(matrix=np.column_stack([y,covariates]),active_covariates=np.array([True,True,False,False]),
                covariate_scales=scales,background=bg,family=fam,pattern_rows=np.arange(n))
    omitted=evaluate(synthetic,arrays,f)
    assert all(omitted['intervals'][name]==dict(status='omitted_constant_covariate',estimate=None) for name in COVARIATES[2:])
    out=Path('results/model_validation/selected-matched-interval-checks-20260928-v1');out.mkdir(parents=True,exist_ok=False)
    write_json(out/'real_input_results.json',results);write_json(out/'omission_fixture.json',omitted)
    result=dict(status='passed_selected_estimate_interval_mapping_checks',real_tree_cases=5,
        original_unit_checks=True,review_flag_retained=True,changed_selected_coefficient_rejected=True,
        omitted_covariates_retained_as_null=2,source_files=sources,
        source_export_sha256=sha(root/'unique_fits.parquet'),reference_receipt_sha256=sha(reference/'receipt.json'),
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/evaluate_selected_matched_intervals.py'),
              Path('scripts/integrate_matched_refinements.py'),Path('scripts/matched_kr_scalar.py'),
              Path('scripts/matched_covariance_information_fast.py'),Path('scripts/matched_kr_covariance.py')]},
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='Selected-row interface checked on original numerical fits plus omitted-column fixture. '
              'Refined full-grid overlay not yet complete; no full-grid interval launch or calibrated inference.')
    write_json(out/'receipt.json',result);write_json(Path('metadata/selected_matched_interval_checks_20260928.json'),result)
    print(result['status'])


if __name__=='__main__':main()
