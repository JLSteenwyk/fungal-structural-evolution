"""Candidate scalar intervals for one selected full-grid fit in original units."""
import numpy as np
from matched_mixed_covariance import MatchedCovariance,profiled_reml
from matched_covariance_information_fast import covariance_information
from matched_kr_scalar import scalar_interval

COVARIATES=('identity_difference','original_coverage_difference',
            'log_aligned_length_ratio','confidence_fraction_difference')


def evaluate(row, arrays, factor):
    result=dict(fit_input_id=row['fit_input_id'],tree=row['tree'],
                selected_source_sha256=row['selected_source_sha256'],selection=row['selection'])
    review=row['selected_review_required']
    if not isinstance(review,(bool,np.bool_)):
        raise ValueError('Explicit selected numerical-review flag required')
    if review:
        return dict(**result,status='selected_fit_requires_review',intervals={})
    matrix=arrays['matrix'];active=np.asarray(arrays['active_covariates'],bool)
    scales=np.asarray(arrays['covariate_scales'],float)
    assert matrix.shape==(row['records'],5) and active.shape==(4,)
    assert np.isfinite(matrix).all()
    np.testing.assert_array_equal(active,np.ptp(matrix[:,1:],axis=0)>1e-12)
    np.testing.assert_allclose(scales,np.std(matrix[:,1:][:,active],axis=0),rtol=1e-13,atol=1e-14)
    assert np.all(scales>0) and np.all(abs(matrix[:,1:][:,~active])<=1e-12)
    names=['intercept']+[name for name,keep in zip(COVARIATES,active) if keep]
    for name,keep in zip(COVARIATES,active):
        if not keep:
            value=row['selected_coefficient_'+name]
            assert value is None or np.isnan(value), 'Omitted covariate must remain null'
    x=np.column_stack([np.ones(len(matrix)),matrix[:,1:][:,active]/scales])
    f=factor[arrays['pattern_rows']]
    ratios=np.array([row['selected_variance_ratio_'+name] for name in ['background','family_component','species']])
    fit=profiled_reml(MatchedCovariance(arrays['background'],arrays['family'],f,1.,*ratios),x,matrix[:,0])
    conversion=np.r_[1.,1/scales]
    raw=np.array([row['selected_intercept']]+[row['selected_coefficient_'+name] for name in names[1:]])
    np.testing.assert_allclose(fit['beta']*conversion,raw,rtol=1e-7,atol=1e-8)
    for key in ['profiled_scale','negative_profiled_reml']:
        np.testing.assert_allclose(fit[key],row['selected_'+key],rtol=1e-7,atol=1e-8)
    np.testing.assert_allclose(fit['conditional_beta_covariance'][0,0],row['selected_conditional_intercept_variance'],rtol=1e-7,atol=1e-8)
    c=covariance_information(arrays['background'],arrays['family'],f,x,
                              row['selected_profiled_scale']*np.r_[1.,ratios])
    np.testing.assert_allclose(c['conditional_beta_covariance'],fit['conditional_beta_covariance'],rtol=1e-7,atol=1e-8)
    intervals={name:dict(status='omitted_constant_covariate',estimate=None) for name,keep in zip(COVARIATES,active) if not keep}
    for i,name in enumerate(names):
        contrast=np.eye(len(names))[i]*conversion[i]
        interval=scalar_interval(c,fit['beta'],contrast)
        if interval['status']=='candidate_interval_pending_coverage_validation':
            np.testing.assert_allclose(interval['estimate'],raw[i],rtol=1e-7,atol=1e-8)
        intervals[name]=interval
    return dict(**result,status='candidate_interval_dispositions_pending_validation',intervals=intervals,
                information_rank=c['numerical_information_rank'],
                information_eigenvalues=c['normalized_information_eigenvalues'].tolist(),
                information_method=c['evaluation_method'],
                scope='Selected numerical estimates, original coefficient units, marginal candidate '
                      'intervals only. Coverage, multiple testing and model adequacy remain unqualified.')
