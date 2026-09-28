"""Attach candidate KR intervals to numerical refits and retain unresolved draws."""
import numpy as np
from matched_covariance_information_fast import covariance_information
from matched_kr_scalar import scalar_interval
from matched_calibration_intervals import exact_coverage_bounds


def evaluate_refit(record, background, family, factor, design, true_beta, alpha=.05):
    x, truth = np.asarray(design, float), np.asarray(true_beta, float)
    if x.ndim != 2 or truth.shape != (x.shape[1],) or not np.isfinite(truth).all():
        raise ValueError('Generating coefficients must match fixed design')
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError('Invalid interval alpha')
    result = dict(fit_id=record['fit_id'], replicate=record['replicate'],
                  response_sha256=record['response_sha256'], true_beta=truth.tolist(),
                  alpha=alpha, status='unresolved', intervals=[])
    if record['status'] != 'refit_numerically_checked':
        if record['status'] not in ['refit_requires_review', 'refit_error_requires_review']:
            raise ValueError('Unknown refit disposition')
        result['reason'] = record['status']
        return result
    fit = record['fit']
    try:
        contractions = covariance_information(background, family, factor, x,
            fit['profiled_scale']*np.array([1., *fit['ratios']]))
        np.testing.assert_allclose(contractions['conditional_beta_covariance'],
                                   fit['conditional_beta_covariance'], rtol=1e-7, atol=1e-8)
        for i in range(len(truth)):
            interval = scalar_interval(contractions, fit['beta'], np.eye(len(truth))[i], alpha)
            if interval['status'] == 'candidate_interval_pending_coverage_validation':
                interval['covers_generating_value'] = bool(interval['lower'] <= truth[i] <= interval['upper'])
            result['intervals'].append(interval)
        result.update(status='interval_dispositions_recorded',
                      information_method=contractions['evaluation_method'],
                      fitted_variances=contractions['absolute_variances'].tolist())
    except Exception as error:
        # A partially completed coefficient vector is not counted as qualified.
        result.update(status='unresolved', reason='interval_evaluation_error',
                      error_type=type(error).__name__, error=str(error), intervals=[])
    return result


def summarize(records, monte_carlo_alpha=.05):
    if not records:
        raise ValueError('Nonempty fixed-size replicate accounting required')
    keys = [(r['fit_id'], r['replicate']) for r in records]
    if len(set(keys)) != len(keys) or len({k[0] for k in keys}) != 1:
        raise ValueError('Distinct replicates from one generating fit required')
    truth, alpha = records[0]['true_beta'], records[0]['alpha']
    covered = np.zeros(len(truth), dtype=int)
    unknown = np.zeros(len(truth), dtype=int)
    for r in records:
        if r['true_beta'] != truth or r['alpha'] != alpha:
            raise ValueError('Cannot pool different generating coefficients or interval levels')
        if r['status'] == 'unresolved':
            unknown += 1
            continue
        if r['status'] != 'interval_dispositions_recorded' or len(r['intervals']) != len(truth):
            raise ValueError('Invalid interval disposition vector')
        for i, interval in enumerate(r['intervals']):
            if interval['status'] != 'candidate_interval_pending_coverage_validation':
                unknown[i] += 1
            else:
                flag = interval['covers_generating_value']
                if type(flag) is not bool:
                    raise ValueError('Explicit boolean coverage outcome required')
                covered[i] += flag
    return dict(attempted=len(records), marginal_coverage_intervals=[
        exact_coverage_bounds(int(k), int(u), len(records), monte_carlo_alpha)
        for k,u in zip(covered, unknown)],
        scope='All attempted draws retained. Fixed-size independent replicates of one generating '
              'model required; no optional-stopping validity or simultaneous coverage claimed.')
