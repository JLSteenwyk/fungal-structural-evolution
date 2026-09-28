"""Replay six full-refit fixtures, evaluate candidate intervals and failure counts."""
import copy
import hashlib
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json
from simulate_matched_working_model import simulate
from matched_mixed_covariance import MatchedCovariance, profiled_reml
from calibrate_matched_kr import evaluate_refit, summarize


def main():
    threadpool_limits(1)
    root = Path('results/model_validation/matched-simulation-refit-checks-20260928-v1')
    receipt = json.loads((root/'receipt.json').read_text())
    assert receipt['status'] == 'passed_synthetic_end_to_end_refit_checks'
    for name, digest in receipt['artifacts'].items():
        assert sha(root/name) == digest
    rng = np.random.default_rng(314159); n = 48
    bg = np.repeat(np.arange(16),3); fam = bg//4; factor = rng.normal(size=(n,4))/2
    x = np.column_stack([np.ones(n),rng.normal(size=n)]); beta = np.array([.1,-.4])
    out = Path('results/model_validation/matched-kr-refit-accounting-20260928-v1')
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for case,ratios in enumerate([[0,0,0],[.4,.8,.2]]):
        batch = []
        for replicate in range(3):
            record = json.loads((root/f'case_{case}_replicate_{replicate}.json').read_text())
            y = simulate(bg,fam,factor,x,beta,1.3,ratios,1,
                np.random.default_rng(np.random.SeedSequence(record['seed_entropy'])))[:,0]
            assert hashlib.sha256(np.asarray(y,dtype='<f8').tobytes()).hexdigest() == record['response_sha256']
            fitted = record['fit']
            check = profiled_reml(MatchedCovariance(bg,fam,factor,1.,*fitted['ratios']),x,y)
            for key in ['beta','profiled_scale','conditional_beta_covariance','negative_profiled_reml']:
                np.testing.assert_allclose(check[key],fitted[key],rtol=1e-7,atol=1e-8)
            result = evaluate_refit(record,bg,fam,factor,x,beta)
            assert result['status'] == 'interval_dispositions_recorded', result
            assert all(i['status'] == 'candidate_interval_pending_coverage_validation' for i in result['intervals'])
            np.testing.assert_array_equal(result['fitted_variances'],fitted['profiled_scale']*np.array([1.,*fitted['ratios']]))
            write_json(out/f'case_{case}_replicate_{replicate}.json',result)
            rows.append(result);batch.append(result)
        summary = summarize(batch)
        assert summary['attempted'] == 3
        assert all(r['unresolved'] == 0 for r in summary['marginal_coverage_intervals'])
    record['replicate'] = 99
    with patch('calibrate_matched_kr.covariance_information',side_effect=RuntimeError('fixture failure')):
        failed = evaluate_refit(record,bg,fam,factor,x,beta)
    assert failed['status'] == 'unresolved'
    mixed = summarize(rows[3:]+[failed])
    assert mixed['attempted'] == 4
    assert all(r['unresolved'] == 1 for r in mixed['marginal_coverage_intervals'])
    for status in ['refit_requires_review','refit_error_requires_review']:
        record['status'] = status
        retained = evaluate_refit(record,bg,fam,factor,x,beta)
        assert retained['status'] == 'unresolved' and retained['reason'] == status
    partial = copy.deepcopy(rows[3]); partial['replicate'] = 98
    partial['intervals'][0] = dict(status='information_requires_review')
    mixed_partial = summarize(rows[3:]+[partial])
    assert [a['unresolved'] for a in mixed_partial['marginal_coverage_intervals']] == [1,0]
    rejected = 0
    different = copy.deepcopy(rows[0]);different['replicate']=97;different['alpha']=.01
    for batch in [[rows[0],rows[0]], [rows[0],rows[3]], [rows[0],different]]:
        try: summarize(batch)
        except ValueError: rejected += 1
        else: raise AssertionError('Invalid pooling accepted')
    result = dict(status='passed_six_refit_kr_interval_and_failure_accounting_checks',
        replayed_responses=6, candidate_intervals=12, failure_summary=mixed,
        partial_coefficient_failure_summary=mixed_partial, invalid_pooling_rejected=rejected,
        source_receipt_sha256=sha(root/'receipt.json'),
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/calibrate_matched_kr.py'),
            Path('scripts/matched_kr_scalar.py'),Path('scripts/matched_kr_covariance.py'),
            Path('scripts/matched_covariance_information_fast.py'),Path('scripts/matched_covariance_information.py'),
            Path('scripts/matched_mixed_covariance.py'),Path('scripts/matched_calibration_intervals.py'),
            Path('scripts/simulate_matched_working_model.py')]},
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='Six saved synthetic refits replayed. Intervals use each refitted covariance, '
              'not generating variances. Plumbing check only; no fitted-data coverage claim '
              'or full-grid calibration launch. Unresolved attempts remain in denominators.')
    write_json(out/'receipt.json',result)
    write_json(Path('metadata/matched_kr_refit_accounting_checks_20260928.json'),result)
    print(result['status'])


if __name__ == '__main__':
    main()
