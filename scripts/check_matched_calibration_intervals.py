"""Enumerate binomial outcomes and partial-observation completions."""
from pathlib import Path
import numpy as np
from scipy.stats import binom
from ancestral_chain_attempt import sha, write_json
from matched_calibration_intervals import exact_coverage_bounds, summarize_replicates


def main():
    minimum_coverage = 1.; cases = 0
    for n in [1, 5, 20, 100]:
        intervals = [exact_coverage_bounds(k,0,n) for k in range(n+1)]
        for probability in np.r_[0.,.001,np.linspace(.01,.99,99),.999,1.]:
            covered = np.array([r['lower']-1e-14<=probability<=r['upper']+1e-14 for r in intervals])
            coverage = float(binom.pmf(np.arange(n+1),n,probability)[covered].sum())
            assert coverage >= .95-1e-12
            minimum_coverage = min(minimum_coverage,coverage); cases += 1
        # All completed-outcome intervals must lie inside the unresolved envelope.
        for known in range(n+1):
            for unknown in range(n-known+1):
                r = exact_coverage_bounds(known,unknown,n)
                assert r['lower'] <= intervals[known]['lower']+1e-14
                assert r['upper'] >= intervals[known+unknown]['upper']-1e-14
                if unknown == n:assert r['lower']==0 and r['upper']==1
    endpoints = [exact_coverage_bounds(k,0,999,alpha=1e-12) for k in [0,1,500,998,999]]
    assert all(np.isfinite([r['lower'],r['upper']]).all() for r in endpoints)
    invalid=[(-1,0,5,.05),(2,4,5,.05),(0,0,0,.05),(1.5,0,5,.05),(1,0,5,0),(1,0,5,float('nan'))]
    for args in invalid:
        try:exact_coverage_bounds(*args)
        except ValueError:pass
        else:raise AssertionError('Invalid interval inputs accepted')
    example=exact_coverage_bounds(940,10,999)
    records=[dict(fit_id='fixture',replicate=0,status='refit_numerically_checked',nominal_095_t_interval_covers_truth=[True,False]),
             dict(fit_id='fixture',replicate=1,status='refit_error_requires_review')]
    summary=summarize_replicates(records,2)
    assert summary['accounting']['attempted']==2 and summary['accounting']['unresolved']==1
    assert summary['marginal_monte_carlo_intervals']==[exact_coverage_bounds(1,1,2),exact_coverage_bounds(0,1,2)]
    result=dict(status='passed_fixed_sample_exact_binomial_interval_checks',
        enumerated_probability_cases=cases,minimum_enumerated_coverage=minimum_coverage,
        unresolved_envelopes_checked=True,replicate_summary_checked=True,invalid_cases_rejected=len(invalid),
        illustrative_counts_only=example,
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/matched_calibration_intervals.py')]},
        scope='Finite enumeration checks of interval implementation, not project coverage. Requires fixed sample size and iid simulation under a single model; no optional stopping or joint multiplicity claim.')
    write_json(Path('metadata/matched_calibration_interval_checks_20260928.json'),result)
    print(result)


if __name__=='__main__':main()
