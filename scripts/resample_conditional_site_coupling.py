#!/usr/bin/env python3
"""Quantify marker influence and paired marker-bootstrap sensitivity of coupling fits."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from audit_genus_codon_trees import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['fits','plan','output']:parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a new immutable output directory')
    plan=json.loads(args.plan.read_text());receipt=json.loads((args.fits/'receipt.json').read_text())
    if sha(args.fits/'receipt.json')!=plan['source_fit_receipt_sha256']:raise ValueError('Source plan differs')
    for name,digest in receipt['artifacts'].items():
        if sha(args.fits/name)!=digest:raise ValueError('Changed fit artifact')
    frame=pd.read_csv(args.fits/'site_covariates.tsv',sep='\t')
    fitrows=pd.read_csv(args.fits/'fit_summary.tsv',sep='\t')
    original=pd.read_csv(args.fits/'all_coefficients.tsv',sep='\t')
    markers=sorted(frame.marker.unique());g=len(markers)
    groups=[np.flatnonzero(frame.marker.to_numpy()==marker) for marker in markers]
    if g!=plan['markers'] or len(fitrows)!=plan['models']:raise ValueError('Analysis grid differs')
    rng=np.random.default_rng(plan['seed']);weights=rng.multinomial(g,np.full(g,1/g),size=plan['marker_bootstrap_draws'])
    if not np.all(weights.sum(axis=1)==g):raise ValueError('Invalid bootstrap counts')
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'marker_resampling_counts.npz',markers=np.asarray(markers),counts=weights)
    summaries=[];influence=[];numeric=[];total_singular=0
    for _,spec in fitrows.iterrows():
        model_id=spec.model_id;rate=spec.rate_model;scale=spec.rsa_scale
        aa=np.log1p(frame['aa_'+rate+'_rate'].to_numpy());rsa=frame[scale+'_rsa_median'].to_numpy()-.25
        y=np.log1p(frame[spec.alphabet_model+'_'+rate+'_rate'].to_numpy())
        columns=['aa_log1p_rate','rsa_centered','aa_by_rsa','observed_fraction','ca_plddt_q25_fraction']
        arrays=[aa,rsa,aa*rsa,frame.observed_fraction.to_numpy(),frame.ca_plddt_q25.to_numpy()/100]
        if spec.controls=='composition_adjusted':
            columns.append('aa_entropy_nats');arrays.append(frame.aa_entropy_nats.to_numpy())
            for amino in 'CDEFGHIKLMNPQRSTVWY':columns.append('aa_fraction_'+amino);arrays.append((frame['aa_count_'+amino]/frame.observed_taxa).to_numpy())
        x=np.column_stack(arrays);xd=np.empty_like(x);yd=np.empty_like(y)
        cross=[];rhs=[]
        for idx in groups:
            xd[idx]=x[idx]-x[idx].mean(axis=0);yd[idx]=y[idx]-y[idx].mean()
            cross.append(xd[idx].T@xd[idx]);rhs.append(xd[idx].T@yd[idx])
        cross=np.stack(cross);rhs=np.stack(rhs);xx=cross.sum(axis=0);xy=rhs.sum(axis=0)
        beta=np.linalg.solve(xx,xy)
        source=original[original.model_id==model_id].set_index('term').loc[columns,'coefficient'].to_numpy()
        if not np.allclose(beta,source,rtol=1e-7,atol=1e-9):raise ValueError('Absorbed coefficients differ from full marker-intercept model')
        # All reference values remain at the full observed distribution.
        transform=np.zeros((5,len(columns)))
        transform[0,0]=transform[1,1]=transform[2,2]=1
        transform[3,0]=1;transform[3,2]=rsa.mean()
        transform[4,1]=1;transform[4,2]=aa.mean()
        point=transform@beta
        omitted=[]
        for j,marker in enumerate(markers):
            b=np.linalg.solve(xx-cross[j],xy-rhs[j]);omitted.append(b)
            effects=transform@b
            for k,name in enumerate(plan['contrasts']):
                influence.append({'model_id':model_id,'omitted_marker':marker,'contrast':name,'estimate_without_marker':effects[k],'full_estimate':point[k],'difference':effects[k]-point[k]})
        omitted=np.stack(omitted)
        xx_boot=(weights@cross.reshape(g,-1)).reshape(len(weights),len(columns),len(columns))
        xy_boot=weights@rhs
        singular=np.linalg.matrix_rank(xx_boot)!=len(columns);total_singular+=int(singular.sum())
        draws=np.full((len(weights),len(columns)),np.nan)
        draws[~singular]=np.linalg.solve(xx_boot[~singular],xy_boot[~singular,:,None])[:,:,0]
        if not np.isfinite(draws[~singular]).all():raise ValueError('Nonfinite estimable bootstrap fit')
        # Expanded-row least squares provides a separate check of one predetermined
        # omission and one predetermined bootstrap sample for every specification.
        chosen=int(hashlib.sha256(model_id.encode()).hexdigest()[:8],16)
        omission=chosen%g;keep=np.concatenate([idx for j,idx in enumerate(groups) if j!=omission])
        direct_loo=np.linalg.lstsq(xd[keep],yd[keep],rcond=None)[0]
        if not np.allclose(direct_loo,omitted[omission],rtol=1e-7,atol=1e-9):raise ValueError('Expanded leave-one-out fit differs')
        candidates=np.flatnonzero(~singular)
        if not len(candidates):raise ValueError('No estimable bootstrap samples')
        sample=int(candidates[chosen%len(candidates)])
        idx=np.concatenate([np.tile(group,int(n)) for group,n in zip(groups,weights[sample]) if n])
        direct_boot=np.linalg.lstsq(xd[idx],yd[idx],rcond=None)[0]
        if not np.allclose(direct_boot,draws[sample],rtol=1e-7,atol=1e-9):raise ValueError('Expanded bootstrap fit differs')
        effects=draws@transform.T;loo_effects=omitted@transform.T
        np.savez_compressed(args.output/(model_id+'.npz'),coefficient_names=np.asarray(columns),coefficients=draws,contrasts=np.asarray(plan['contrasts']),contrast_values=effects,singular=singular,leave_one_out_coefficients=omitted)
        for k,name in enumerate(plan['contrasts']):
            values=effects[~singular,k];low,high=np.quantile(values,[.025,.975])
            deltas=loo_effects[:,k]-point[k];most=int(np.argmax(np.abs(deltas)))
            summaries.append({'model_id':model_id,'alphabet_model':spec.alphabet_model,'rate_model':rate,'rsa_scale':scale,'controls':spec.controls,'contrast':name,'point_estimate':point[k],'bootstrap_percentile_lower':low,'bootstrap_percentile_upper':high,'bootstrap_estimable_draws':len(values),'bootstrap_singular_draws':int(singular.sum()),'loo_minimum':loo_effects[:,k].min(),'loo_maximum':loo_effects[:,k].max(),'loo_sign_reversals':int(np.sum(loo_effects[:,k]*point[k]<0)),'most_influential_marker':markers[most],'largest_absolute_loo_change':abs(deltas[most]),'observed_mean_log1p_aa_rate':aa.mean(),'observed_mean_rsa':rsa.mean()+.25})
        numeric.append({'model_id':model_id,'point_coefficient_max_difference':float(np.max(abs(beta-source))),'verified_omission':markers[omission],'omission_max_difference':float(np.max(abs(direct_loo-omitted[omission]))),'verified_bootstrap_index_zero_based':sample,'bootstrap_max_difference':float(np.max(abs(direct_boot-draws[sample]))),'maximum_bootstrap_normal_equation_residual':float(np.max(abs(np.einsum('bij,bj->bi',xx_boot[~singular],draws[~singular])-xy_boot[~singular])))})
        print('Resampled',model_id,'singular',int(singular.sum()),flush=True)
    for name,rows in [('contrast_summary.tsv',summaries),('marker_influence.tsv',influence),('numerical_checks.tsv',numeric)]:pd.DataFrame(rows).to_csv(args.output/name,sep='\t',index=False)
    result={'status':'complete_conditional_marker_resampling','models':len(fitrows),'markers':g,'bootstrap_fits':len(fitrows)*len(weights),'singular_bootstrap_fits':total_singular,'leave_one_marker_out_fits':len(fitrows)*g,'contrasts':len(summaries),'plan_sha256':sha(args.plan),'source_fit_receipt_sha256':sha(args.fits/'receipt.json'),'script_sha256':sha(Path(__file__)),'artifacts':{p.name:sha(p) for p in args.output.iterdir()},'interpretation':plan['interpretation'],'interval_interpretation':plan['intervals'],'verification':'All point coefficients match full marker-intercept fits. One deterministic omission and one deterministic resample per model independently match expanded-row least squares. Counts and every bootstrap solve retained, with singularity flags.'}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))


if __name__=='__main__':main()
