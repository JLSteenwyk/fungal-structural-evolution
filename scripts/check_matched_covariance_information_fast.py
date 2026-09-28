"""Compare residual-trace shortcut with dense and original blocked contractions."""
import itertools
import time
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json
from matched_covariance_information_fast import covariance_information
from matched_covariance_information import covariance_information as blocked


def main():
    threadpool_limits(1)
    rng = np.random.default_rng(20260928)
    n = 48
    bg = np.repeat(np.arange(12), 4)
    fam = bg // 3
    x = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
    raw = rng.normal(size=(n, 4)) / 2
    cases = 0
    maximum = 0.
    maximum_relative_scale_error = 0.
    maximum_blocked_error = 0.
    factors = [raw, np.column_stack([raw, raw[:, 0]]), np.zeros((n, 0))]
    for f in factors:
        kernels = [np.eye(n), (bg[:, None] == bg).astype(float),
                   (fam[:, None] == fam).astype(float), f @ f.T]
        for ratios in itertools.product([0., .7, 10000.], repeat=3):
            for scale in [.1, 3.]:
                v = scale * np.array([1., *ratios])
                actual = covariance_information(bg, fam, f, x, v, block_size=7)
                old = blocked(bg, fam, f, x, v, block_size=7)
                np.testing.assert_allclose(actual['expected_reml_information'],
                                           old['expected_reml_information'], rtol=1e-7, atol=1e-7)
                maximum_blocked_error = max(maximum_blocked_error, float(np.max(abs(actual['expected_reml_information'] - old['expected_reml_information']))))
                dense = sum(a * k for a, k in zip(v, kernels))
                inverse = np.linalg.inv(dense)
                t = inverse @ x
                phi = np.linalg.inv(x.T @ t)
                p = inverse - t @ phi @ t.T
                expected = np.array([[.5 * np.trace(p @ a @ p @ b) for b in kernels] for a in kernels])
                first = np.array([-t.T @ a @ t for a in kernels])
                second = np.array([[t.T @ a @ inverse @ b @ t for b in kernels] for a in kernels])
                for key, reference in [('conditional_beta_covariance', phi),
                                       ('first_contractions', first),
                                       ('second_contractions', second),
                                       ('expected_reml_information', expected)]:
                    np.testing.assert_allclose(actual[key], reference, rtol=1e-6, atol=1e-6)
                    maximum = max(maximum, float(np.max(abs(actual[key] - reference))))
                    maximum_relative_scale_error = max(maximum_relative_scale_error, float(np.max(abs(actual[key] - reference)/(1 + abs(reference)))))
                np.testing.assert_array_equal(actual['exact_zero_components'], v == 0)
                if not f.shape[1]:
                    assert actual['numerical_information_rank'] < 4
                cases += 1
    # Exact confounding: species kernel equals the family kernel.
    f = (fam[:, None] == np.unique(fam)).astype(float)
    result = covariance_information(bg, fam, f, x, [1., .2, .3, .4], block_size=5)
    assert result['numerical_information_rank'] == 3
    np.testing.assert_allclose(result['expected_reml_information'][2], result['expected_reml_information'][3])
    # Blocking changes storage only, including blocks larger than n.
    other = covariance_information(bg, fam, f, x, [1., .2, .3, .4], block_size=100)
    np.testing.assert_allclose(result['expected_reml_information'], other['expected_reml_information'], rtol=1e-12, atol=1e-12)
    fallback = covariance_information(bg, fam, np.eye(n), x, [1., 0., 0., 1e8])
    assert fallback['evaluation_method'] == 'blocked_reference_cancellation_fallback'
    old = blocked(bg, fam, np.eye(n), x, [1., 0., 0., 1e8])
    np.testing.assert_array_equal(fallback['expected_reml_information'], old['expected_reml_information'])
    invalid = 0
    for changes in [dict(block_size=0), dict(block_size=True), dict(variances=[0, 1, 1, 1]),
                    dict(variances=[1, -1, 1, 1]), dict(variances=[1, 2, 3]),
                    dict(design=np.ones((n, 2))), dict(design=np.full((n, 1), np.nan))]:
        args = dict(background=bg, family=fam, factor=f, design=x, variances=[1, .2, .3, .4])
        args.update(changes)
        try:
            covariance_information(**args)
        except ValueError:
            invalid += 1
        else:
            raise AssertionError('Invalid input accepted')
    n = 600
    bg = np.repeat(np.arange(100), 6); fam = bg // 5
    f = rng.normal(size=(n, 242))/16
    x = np.column_stack([np.ones(n), rng.normal(size=(n, 4))])
    timings = {}; results = {}
    for name, fn in [('blocked', blocked), ('shortcut', covariance_information)]:
        start = time.perf_counter()
        results[name] = fn(bg, fam, f, x, [.7, .2, .4, .3])
        timings[name] = time.perf_counter()-start
    np.testing.assert_allclose(results['shortcut']['expected_reml_information'], results['blocked']['expected_reml_information'], rtol=1e-9, atol=1e-8)
    record = dict(status='passed_residual_trace_shortcut_checks', cancellation_fallback_checked=True, dense_cases=cases,
        maximum_absolute_disagreement=maximum, maximum_scaled_disagreement=maximum_relative_scale_error,
        maximum_blocked_information_disagreement=maximum_blocked_error,
        synthetic_600_record_timing_seconds=timings, exact_confounding_rank=result['numerical_information_rank'],
        invalid_inputs_rejected=invalid, block_partition_equivalence=True,
        pins={str(p): sha(p) for p in [Path(__file__), Path('scripts/matched_covariance_information.py'), Path('scripts/matched_covariance_information_fast.py'),
                                     Path('scripts/matched_mixed_covariance.py')]},
        scope='Synthetic shortcut, blocked and dense comparisons including large and zero variances. No information '
              'inverse, KR intervals or coverage claims. Blocked exact reference has quadratic '
              'work and is not approved for full-grid execution.')
    write_json(Path('metadata/matched_covariance_information_fast_checks_20260928.json'), record)
    print(record)


if __name__ == '__main__':
    main()
