#!/usr/bin/env python3
"""Fit conditional sequence/structure site associations with marker fixed effects."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from audit_genus_codon_trees import sha


def checked(path):
    receipt=json.loads((path/'receipt.json').read_text())
    for name,digest in receipt['artifacts'].items():
        if sha(path/name)!=digest:raise ValueError('Changed artifact: '+name)
    return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['frame','accessibility','plan','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a new immutable output directory')
    plan=json.loads(args.plan.read_text());fr=checked(args.frame);ar=checked(args.accessibility)
    if sha(args.frame/'receipt.json')!=plan['source_frame_sha256'] or sha(args.accessibility/'receipt.json')!=plan['source_accessibility_sha256']:raise ValueError('Source provenance differs from plan')
    sites=pd.read_csv(args.frame/'site_rate_exposure.tsv',sep='\t')
    keys=['marker','paired_column_1based']
    if sites.duplicated(keys).any() or len(sites)!=plan['sites'] or sites.marker.nunique()!=plan['markers']:raise ValueError('Unexpected site universe')
    cols=keys+['matrix_column_1based','taxon_id','ca_plddt','normalization_status','rsa_tien2013_theoretical','rsa_miller1987']
    raw=pd.read_csv(args.accessibility/'normalized_paired_sites.tsv.gz',sep='\t',usecols=cols)
    if raw.duplicated(keys+['taxon_id']).any() or len(raw)!=fr['taxon_site_observations']:raise ValueError('Taxon-site grid differs')
    if not raw.normalization_status.eq('internal_residue_normalized').all():raise ValueError('Unnormalized exposure')
    if not np.isfinite(raw[['ca_plddt','rsa_tien2013_theoretical','rsa_miller1987']]).all().all():raise ValueError('Nonfinite covariate')
    group=raw.groupby(keys,sort=True)
    confidence=group.ca_plddt.quantile(.25).rename('ca_plddt_q25').reset_index()
    source_counts=group.size().rename('source_observed_taxa').reset_index()
    source_matrix=group.matrix_column_1based.agg(['min','max']).reset_index()
    joined=sites.merge(confidence,on=keys,validate='one_to_one').merge(source_counts,on=keys,validate='one_to_one').merge(source_matrix,on=keys,validate='one_to_one')
    if len(joined)!=len(sites) or not joined.observed_taxa.eq(joined.source_observed_taxa).all() or not joined.matrix_column_1based.eq(joined['min']).all() or not joined['min'].eq(joined['max']).all():raise ValueError('Confidence/exposure join differs')
    for scale,column in [('tien2013','rsa_tien2013_theoretical'),('miller1987','rsa_miller1987')]:
        med=group[column].median().rename('recomputed_rsa').reset_index()
        check=joined.merge(med,on=keys,validate='one_to_one')
        if not np.allclose(check[scale+'_rsa_median'],check.recomputed_rsa,atol=1e-12,rtol=0):raise ValueError('RSA medians differ')
    joined=joined.drop(columns=['source_observed_taxa','min','max']).sort_values(keys).reset_index(drop=True)
    policy=plan['marker_review_policy']
    review=joined.loc[joined.marker_review_status!='no_current_copy_review_flag'].groupby(['marker','marker_review_status']).size()
    observed=[{'marker':str(marker),'status':str(status),'sites':int(count)} for (marker,status),count in review.items()]
    if observed!=policy['expected_flagged_markers']:raise ValueError('Copy-review universe differs from explicitly reviewed plan')
    if policy['mode']=='exclude_reviewed_markers':
        joined=joined[~joined.marker.isin([r['marker'] for r in observed])].copy().reset_index(drop=True)
    elif policy['mode']!='retain_flagged_exploratory':raise ValueError('Unknown copy-review policy')
    if len(joined)!=plan['analysis_sites'] or joined.marker.nunique()!=plan['analysis_markers']:raise ValueError('Reviewed analysis universe differs')
    markers=sorted(joined.marker.unique());groups=pd.Categorical(joined.marker,categories=markers).codes
    dummy=pd.get_dummies(joined.marker,prefix='marker',dtype=float)
    focal=['aa_log1p_rate','rsa_centered','aa_by_rsa']
    coefficients=[];summaries=[];allcoeff=[]
    for alphabet,rate,scale,controls in itertools.product(plan['grid']['alphabet_model'],plan['grid']['rate_model'],plan['grid']['rsa_scale'],plan['grid']['controls']):
        model_id='__'.join([alphabet,rate,scale,controls])
        aa=np.log1p(joined['aa_'+rate+'_rate'].to_numpy())
        response=np.log1p(joined[alphabet+'_'+rate+'_rate'].to_numpy())
        rsa=joined[scale+'_rsa_median'].to_numpy()-.25
        x=pd.DataFrame({'aa_log1p_rate':aa,'rsa_centered':rsa,'aa_by_rsa':aa*rsa,'observed_fraction':joined.observed_fraction,'ca_plddt_q25_fraction':joined.ca_plddt_q25/100})
        if controls=='composition_adjusted':
            x['aa_entropy_nats']=joined.aa_entropy_nats
            for amino in 'CDEFGHIKLMNPQRSTVWY':x['aa_fraction_'+amino]=joined['aa_count_'+amino]/joined.observed_taxa
        x=pd.concat([x,dummy],axis=1)
        values=x.to_numpy(dtype=float);n,k=values.shape
        if not np.isfinite(values).all() or not np.isfinite(response).all():raise ValueError('Nonfinite design')
        rank=np.linalg.matrix_rank(values)
        if rank!=k:raise ValueError('Rank-deficient model: '+model_id)
        fitted=sm.OLS(response,values).fit(cov_type='cluster',cov_kwds={'groups':groups,'use_correction':True,'df_correction':True},use_t=True)
        # Independent least-squares solution and sandwich assembly verify the library result.
        beta=np.linalg.lstsq(values,response,rcond=None)[0]
        residual=response-values@beta
        if not np.allclose(beta,fitted.params,rtol=1e-7,atol=1e-9):raise ValueError('Independent coefficients differ')
        pinv=np.linalg.pinv(values);bread=pinv@pinv.T
        scores=np.stack([values[groups==g].T@residual[groups==g] for g in range(len(markers))])
        correction=len(markers)/(len(markers)-1)*(n-1)/(n-k)
        covariance=correction*bread@(scores.T@scores)@bread
        if not np.allclose(covariance,fitted.cov_params(),rtol=1e-7,atol=1e-10):raise ValueError('Independent clustered covariance differs')
        intervals=fitted.conf_int();within_y=response-pd.Series(response).groupby(groups).transform('mean').to_numpy()
        within_r2=1-float(residual@residual)/float(within_y@within_y)
        summary={'model_id':model_id,'alphabet_model':alphabet,'rate_model':rate,'rsa_scale':scale,'controls':controls,'sites':n,'markers':len(markers),'parameters_including_marker_intercepts':k,'design_rank':rank,'design_condition_number':float(np.linalg.cond(values)),'cluster_reference_df':float(fitted.df_resid_inference),'within_marker_r_squared':within_r2,'residual_sum_squares':float(residual@residual),'independent_coefficient_max_difference':float(np.max(np.abs(beta-fitted.params))),'independent_covariance_max_difference':float(np.max(np.abs(covariance-fitted.cov_params())))}
        summaries.append(summary)
        for j,name in enumerate(x.columns):
            row={'model_id':model_id,'alphabet_model':alphabet,'rate_model':rate,'rsa_scale':scale,'controls':controls,'term':name,'coefficient':float(fitted.params[j]),'marker_cluster_standard_error':float(fitted.bse[j]),'t_statistic':float(fitted.tvalues[j]),'p_value_conditional_cluster_t':float(fitted.pvalues[j]),'ci95_lower':float(intervals[j,0]),'ci95_upper':float(intervals[j,1])}
            allcoeff.append(row)
            if name in focal:coefficients.append(dict(row))
        print('Fitted',model_id,flush=True)
    q=multipletests([r['p_value_conditional_cluster_t'] for r in coefficients],method='fdr_bh')[1]
    for row,value in zip(coefficients,q):row['bh_q_across_all_72_focal_tests']=float(value)
    if len(coefficients)!=72 or len(summaries)!=plan['fits']:raise ValueError('Incomplete planned fit grid')
    args.output.mkdir(parents=True)
    joined.to_csv(args.output/'site_covariates.tsv',sep='\t',index=False)
    for name,rows in [('focal_coefficients.tsv',coefficients),('all_coefficients.tsv',allcoeff),('fit_summary.tsv',summaries)]:pd.DataFrame(rows).to_csv(args.output/name,sep='\t',index=False)
    result={'status':'complete_exploratory_conditional_site_coupling_fits','sites':len(joined),'markers':len(markers),'fits':len(summaries),'focal_tests':len(coefficients),'confidence_taxon_site_rows':len(raw),'plan_sha256':sha(args.plan),'source_frame_sha256':sha(args.frame/'receipt.json'),'source_accessibility_sha256':sha(args.accessibility/'receipt.json'),'script_sha256':sha(Path(__file__)),'versions':{'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__},'artifacts':{p.name:sha(p) for p in args.output.iterdir()},'interpretation':plan['interpretation'],'marker_review_policy':policy,'source_frame_sites':len(sites),'verification':'Full taxon/site confidence join and both RSA medians checked; every fit cross-checked by independent NumPy least squares and manually assembled finite-sample-corrected marker-cluster covariance. These numerical checks do not validate the statistical assumptions.'}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))


if __name__=='__main__':main()
