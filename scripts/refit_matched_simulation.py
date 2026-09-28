"""One reproducible full refit per simulated response, retaining review outcomes."""
import hashlib
import numpy as np
from scipy.stats import t
from simulate_matched_working_model import simulate
from refine_matched_reml_analytic import refine


def response_seed(master_seed, fit_id, replicate):
    if type(master_seed) is not int or master_seed < 0 or type(replicate) is not int or replicate < 0:
        raise ValueError('Nonnegative integer seed and replicate index required')
    if not isinstance(fit_id, str) or not fit_id:
        raise ValueError('Nonempty source fit identifier required')
    words = np.frombuffer(hashlib.sha256(fit_id.encode()).digest(), dtype='<u4').tolist()
    return [master_seed, replicate, *words]


def refit_one(background, family, factor, design, beta, scale, ratios,
              master_seed, fit_id, replicate, maximum_ratio=10000., maxiter=1000):
    entropy = response_seed(master_seed, fit_id, replicate)
    ratios = np.asarray(ratios, dtype=float)
    if ratios.shape != (3,) or not np.isfinite(ratios).all() or (ratios < 0).any() or (ratios > maximum_ratio).any():
        raise ValueError('Generating ratios outside declared fit bounds')
    response = simulate(background, family, factor, design, beta, scale, ratios,
                        1, np.random.default_rng(np.random.SeedSequence(entropy)))[:, 0]
    record = dict(fit_id=fit_id, replicate=replicate, seed_entropy=entropy,
                  response_sha256=hashlib.sha256(np.asarray(response,dtype='<f8').tobytes()).hexdigest())
    try:
        fitted = refine(background, family, factor, design, response, np.log1p(ratios), maximum_ratio, maxiter)
        record['fit'] = fitted
        if fitted['status'] != 'refinement_passed_numerical_checks':
            record['status'] = 'refit_requires_review'
            return record
        covariance = np.asarray(fitted['conditional_beta_covariance'])
        variance = np.diag(covariance)
        if not np.isfinite(variance).all() or (variance <= 0).any():
            raise ValueError('Invalid fitted coefficient variance')
        se = np.sqrt(variance)
        standardized = (np.asarray(fitted['beta'])-np.asarray(beta))/se
        critical = float(t.ppf(.975, fitted['residual_degrees_of_freedom']))
        record.update(status='refit_numerically_checked', coefficient_errors=(np.asarray(fitted['beta'])-beta).tolist(),
                      standard_errors=se.tolist(), studentized_errors=standardized.tolist(),
                      nominal_t_critical=critical, nominal_095_t_interval_covers_truth=(abs(standardized)<=critical).tolist())
    except Exception as error:
        record.update(status='refit_error_requires_review', error_type=type(error).__name__, error=str(error))
    return record


def coverage_accounting(records, coefficients):
    if not records or type(coefficients) is not int or coefficients < 1:
        raise ValueError('Nonempty replicate accounting required')
    keys = [(r['fit_id'],r['replicate']) for r in records]
    if len(set(keys)) != len(keys) or len({key[0] for key in keys}) != 1:
        raise ValueError('Unique replicates of one source fit required')
    good = [r for r in records if r['status']=='refit_numerically_checked']
    assert all(r['status'] in ['refit_numerically_checked','refit_requires_review','refit_error_requires_review'] for r in records)
    covers = np.zeros(coefficients,dtype=int)
    for r in good:
        values = r['nominal_095_t_interval_covers_truth']
        assert len(values)==coefficients and all(type(v) is bool for v in values)
        covers += values
    unknown = len(records)-len(good)
    return dict(attempted=len(records),qualified=len(good),unresolved=unknown,
                covered_qualified=covers.tolist(),
                coverage_lower_if_all_unresolved_fail=(covers/len(records)).tolist(),
                coverage_upper_if_all_unresolved_cover=((covers+unknown)/len(records)).tolist(),
                scope='Bounds from unresolved refits only, not Monte Carlo confidence intervals. Nominal t intervals are being evaluated, not asserted to have 95% coverage.')
