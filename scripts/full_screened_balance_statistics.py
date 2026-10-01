"""Finite-pair descriptive balance, with exact constant-value variance handling."""
import numpy as np
from background_control_balance import FEATURES, feature_vector


def balance(x, y, baseline):
    keep = np.isfinite(x) & np.isfinite(y); x = x[keep]; y = y[keep]; baseline = baseline[np.isfinite(baseline)]
    n, nb = len(x), len(baseline)
    mean = lambda a: float(a[0]) if len(a) and np.all(a == a[0]) else float(a.mean())
    sd = lambda a: 0.0 if np.all(a == a[0]) else float(a.std(ddof=1))
    bm = mean(baseline) if nb else ''
    out = dict(pairs=n,baseline_targets=nb,target_mean='',control_mean='',target_sd='',control_sd='',mean_difference='',standardized_mean_difference='',smd_status='no_pairs',mean_absolute_difference='',p95_absolute_difference='',maximum_absolute_difference='',baseline_target_mean=bm,selection_mean_shift='',selection_shift_in_baseline_sd='',selection_shift_status='no_pairs')
    if not n: return out
    xm, ym = mean(x), mean(y); absolute = np.abs(x-y)
    out.update(target_mean=xm,control_mean=ym,mean_difference=xm-ym,mean_absolute_difference=float(absolute.mean()),p95_absolute_difference=float(np.quantile(absolute,.95)),maximum_absolute_difference=float(absolute.max()))
    if n < 2: out['smd_status'] = 'insufficient_pairs'
    else:
        sx, sy = sd(x), sd(y); scale = np.hypot(sx,sy)/np.sqrt(2); out.update(target_sd=sx,control_sd=sy)
        if scale == 0: out['smd_status'] = 'zero_pooled_variance'
        else: out.update(standardized_mean_difference=(xm-ym)/scale,smd_status='estimable')
    if not nb: out['selection_shift_status'] = 'no_baseline'
    else:
        shift = xm-bm; out['selection_mean_shift'] = shift
        if nb < 2: out['selection_shift_status'] = 'insufficient_baseline'
        elif sd(baseline) == 0: out['selection_shift_status'] = 'zero_baseline_variance'
        else: out.update(selection_shift_in_baseline_sd=shift/sd(baseline),selection_shift_status='estimable')
    return out
