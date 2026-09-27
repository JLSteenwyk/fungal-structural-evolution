"""Descriptive covariate balance with explicit zero/insufficient-variance outcomes."""
import numpy as np
FEATURES=['sequence_distance','log_positive_sequence_distance','mean_log_length','log_length_asymmetry','mean_plddt','minimum_plddt','mean_lowconf_fraction','maximum_lowconf_fraction']
def feature_vector(n):
    lengths=np.array([n['length_a'],n['length_b']],float);p=np.array([n['mean_ca_plddt_a'],n['mean_ca_plddt_b']],float);f=np.array([n['fraction_ca_plddt_below50_a'],n['fraction_ca_plddt_below50_b']],float);d=n['sequence_distance']
    assert np.isfinite(lengths).all() and (lengths>0).all() and d>=0
    return np.array([d,np.log(d) if d>0 else np.nan,np.log(lengths).mean(),abs(np.log(lengths[0]/lengths[1])),p.mean(),p.min(),f.mean(),f.max()])
def balance(x,y,baseline):
    keep=np.isfinite(x)&np.isfinite(y);x=x[keep];y=y[keep];baseline=baseline[np.isfinite(baseline)];n=len(x);nb=len(baseline)
    out=dict(pairs=n,baseline_targets=nb,target_mean='',control_mean='',target_sd='',control_sd='',mean_difference='',standardized_mean_difference='',smd_status='no_pairs',mean_absolute_difference='',p95_absolute_difference='',maximum_absolute_difference='',baseline_target_mean=float(baseline.mean()) if nb else '',selection_mean_shift='',selection_shift_in_baseline_sd='',selection_shift_status='no_pairs')
    if not n:return out
    xm=float(x.mean());ym=float(y.mean());diff=x-y;absolute=np.abs(diff)
    out.update(target_mean=xm,control_mean=ym,mean_difference=xm-ym,mean_absolute_difference=float(absolute.mean()),p95_absolute_difference=float(np.quantile(absolute,.95)),maximum_absolute_difference=float(absolute.max()))
    if n<2:out['smd_status']='insufficient_pairs'
    else:
        vx=float(x.var(ddof=1));vy=float(y.var(ddof=1));scale=((vx+vy)/2)**.5;out.update(target_sd=vx**.5,control_sd=vy**.5)
        if scale==0:out['smd_status']='zero_pooled_variance'
        else:out.update(standardized_mean_difference=(xm-ym)/scale,smd_status='estimable')
    if nb:
        shift=xm-float(baseline.mean());out['selection_mean_shift']=shift
        if nb<2:out['selection_shift_status']='insufficient_baseline'
        else:
            sd=float(baseline.std(ddof=1))
            if sd==0:out['selection_shift_status']='zero_baseline_variance'
            else:out.update(selection_shift_in_baseline_sd=shift/sd,selection_shift_status='estimable')
    return out
