#!/usr/bin/env python3
"""Paired conditional bootstrap differences at common baseline references."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

TERMS = ['aa_log1p_rate','rsa_centered','aa_by_rsa']
CONTRASTS = ['aa_at_RSA_025','rsa_at_aa_log_rate_0','interaction',
             'aa_at_baseline_mean_RSA','rsa_at_baseline_mean_aa_log_rate']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked(folder):
    r=json.loads((folder/'receipt.json').read_text())
    for name,digest in r['artifacts'].items():
        if sha(folder/name)!=digest:raise ValueError('Changed source artifact: '+name)
    return r


def pair_counts(base,sensitivity):
    if (not np.array_equal(base['markers'],sensitivity['markers'])
            or not np.array_equal(base['counts'],sensitivity['counts'])):
        raise ValueError('Bootstrap marker identities or paired draws differ')


def transform(mean_rsa,mean_aa):
    return np.array([[1,0,0],[0,1,0],[0,0,1],[1,0,mean_rsa-.25],[0,1,mean_aa]],dtype=float)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['baseline','sensitivity','baseline-resampling','sensitivity-resampling','output']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    folders=[a.baseline,a.sensitivity,a.baseline_resampling,a.sensitivity_resampling]
    receipts=[checked(x) for x in folders]
    b,s,br,sr=receipts
    if (any(x['status']!='complete_exploratory_conditional_site_coupling_fits' for x in [b,s])
            or any(x['status']!='complete_conditional_marker_resampling' for x in [br,sr])
            or br['source_fit_receipt_sha256']!=sha(a.baseline/'receipt.json')
            or sr['source_fit_receipt_sha256']!=sha(a.sensitivity/'receipt.json')
            or b['markers']!=s['markers'] or b['sites']!=s['sites']
            or b['marker_review_policy']['mode']!=s['marker_review_policy']['mode']
            or b['marker_review_policy']['expected_flagged_markers']!=s['marker_review_policy']['expected_flagged_markers']):
        raise ValueError('Mismatched analysis provenance or policy')
    pins={str(x/'receipt.json'):sha(x/'receipt.json') for x in folders}
    with np.load(a.baseline_resampling/'marker_resampling_counts.npz',allow_pickle=False) as bc, np.load(a.sensitivity_resampling/'marker_resampling_counts.npz',allow_pickle=False) as sc:
        pair_counts(bc,sc);draws=bc['counts'].shape[0];markers=bc['markers'].tolist()
        if len(markers)!=b['markers'] or draws*24!=br['bootstrap_fits'] or br['bootstrap_fits']!=sr['bootstrap_fits']:
            raise ValueError('Incomplete bootstrap dimensions')
    read=lambda p:pd.read_csv(p,sep='\t',float_precision='round_trip')
    bf,sf=(read(x/'site_covariates.tsv') for x in [a.baseline,a.sensitivity])
    keys=['marker','paired_column_1based','matrix_column_1based']
    pd.testing.assert_frame_equal(bf[keys],sf[keys],check_exact=True)
    if sorted(bf.marker.unique())!=markers or len(bf)!=b['sites']:raise ValueError('Changed site universe')
    specs=read(a.baseline/'fit_summary.tsv')
    other=read(a.sensitivity/'fit_summary.tsv')
    columns=['model_id','alphabet_model','rate_model','rsa_scale','controls']
    pd.testing.assert_frame_equal(specs[columns],other[columns],check_exact=True)
    if len(specs)!=24 or specs.model_id.duplicated().any():raise ValueError('Wrong specification grid')
    points=[read(x/'focal_coefficients.tsv').set_index(['model_id','term']) for x in [a.baseline,a.sensitivity]]
    rows=[];arrays={}
    for spec in specs.itertuples():
        mean_rsa=float(bf[spec.rsa_scale+'_rsa_median'].mean())
        mean_aa=float(np.log1p(bf['aa_'+spec.rate_model+'_rate']).mean())
        t=transform(mean_rsa,mean_aa)
        with np.load(a.baseline_resampling/(spec.model_id+'.npz'),allow_pickle=False) as bn, np.load(a.sensitivity_resampling/(spec.model_id+'.npz'),allow_pickle=False) as sn:
            if not np.array_equal(bn['coefficient_names'],sn['coefficient_names']):raise ValueError('Different coefficient order')
            names=bn['coefficient_names'].tolist();ii=[names.index(term) for term in TERMS]
            bb,ss=bn['coefficients'][:,ii],sn['coefficients'][:,ii]
            if bb.shape!=ss.shape or bb.shape!=(draws,3):raise ValueError('Wrong draw grid')
            eligible=~(bn['singular']|sn['singular'])
            if not eligible.any() or not np.isfinite(bb[eligible]).all() or not np.isfinite(ss[eligible]).all():raise ValueError('Invalid paired coefficients')
            delta=np.full((draws,5),np.nan);delta[eligible]=(ss[eligible]-bb[eligible])@t.T
        point=[np.array([x.loc[(spec.model_id,term),'coefficient'] for term in TERMS])@t.T for x in points]
        arrays[spec.model_id]=(delta,eligible)
        for j,contrast in enumerate(CONTRASTS):
            lo,hi=np.quantile(delta[eligible,j],[.025,.975])
            rows.append(dict(model_id=spec.model_id,alphabet_model=spec.alphabet_model,rate_model=spec.rate_model,
                             rsa_scale=spec.rsa_scale,controls=spec.controls,contrast=contrast,
                             baseline_mean_rsa=mean_rsa,baseline_mean_aa_log1p_rate=mean_aa,
                             baseline_estimate=point[0][j],sensitivity_estimate=point[1][j],
                             sensitivity_minus_baseline=point[1][j]-point[0][j],paired_draws=int(eligible.sum()),
                             excluded_singular_draws=int((~eligible).sum()),percentile_lower=lo,percentile_upper=hi))
    a.output.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(a.output/'paired_difference_summary.tsv',sep='\t',index=False)
    for model,(delta,eligible) in arrays.items():
        np.savez_compressed(a.output/(model+'.npz'),contrasts=np.array(CONTRASTS),difference=delta,eligible=eligible)
    for folder in folders:
        checked(folder)
        if sha(folder/'receipt.json')!=pins[str(folder/'receipt.json')]:raise ValueError('Source changed during comparison')
    result=dict(status='complete_paired_conditional_fcs_bootstrap_comparison',models=24,markers=b['markers'],sites=b['sites'],
                draws_per_model=draws,contrasts=len(rows),source_receipts=pins,script_sha256=sha(__file__),
                artifacts={p.name:sha(p) for p in a.output.iterdir()},
                scope='Identical marker identities and multinomial counts paired across baseline/FCS draws. Both contrasts use fixed baseline site-weighted reference means. Unadjusted percentile difference intervals; no multiple-testing significance claim, equivalence test or causal/lineage-wide inference. Conditional on estimated rates/topologies and predicted structures; cross-marker phylogenetic dependence remains.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['artifacts','source_receipts']},indent=2))


if __name__=='__main__':main()
