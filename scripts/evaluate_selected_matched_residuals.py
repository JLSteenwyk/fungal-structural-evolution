"""Bind descriptive residuals to the full selected-estimate grid."""
import numpy as np
from evaluate_selected_matched_intervals import COVARIATES
from matched_marginal_residual_diagnostics import diagnostics


def evaluate(row,arrays,factor):
    identity={k:row[k] for k in ['fit_input_id','tree','selected_source_sha256','selection']}
    if row['selected_review_required']:
        return dict(**identity,status='selected_fit_requires_review')
    matrix=arrays['matrix'];active=arrays['active_covariates'];scales=arrays['covariate_scales']
    assert matrix.shape==(row['records'],5) and active.shape==(4,)
    np.testing.assert_array_equal(active,np.ptp(matrix[:,1:],axis=0)>1e-12)
    np.testing.assert_allclose(scales,np.std(matrix[:,1:][:,active],axis=0),rtol=1e-13,atol=1e-14)
    assert np.all(scales>0) and np.all(abs(matrix[:,1:][:,~active])<=1e-12)
    names=['intercept']+[name for name,keep in zip(COVARIATES,active) if keep]
    x=np.column_stack([np.ones(len(matrix)),matrix[:,1:][:,active]/scales])
    ratios=np.array([row['selected_variance_ratio_'+name] for name in ['background','family_component','species']])
    result=diagnostics(arrays['background'],arrays['family'],factor[arrays['pattern_rows']],x,matrix[:,0],row['selected_profiled_scale']*np.r_[1.,ratios])
    if result['status']=='descriptive_marginal_residual_diagnostics':
        raw=np.asarray(result.pop('beta'))*np.r_[1.,1/scales]
        expected=[row['selected_intercept']]+[row['selected_coefficient_'+name] for name in names[1:]]
        np.testing.assert_allclose(raw,expected,rtol=1e-7,atol=1e-8)
        result['raw_unit_beta']=raw.tolist();result['coefficient_names']=names
        for item in result['covariate_bins']:item['covariate']=names[item['design_column']]
    return dict(**identity,**result)
