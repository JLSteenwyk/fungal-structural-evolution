#!/usr/bin/env python3
"""Compare a separate direct-lag implementation with the locked FFT oracle."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
import warnings

import arviz
import numpy as np
from ancestral_chain_diagnostics import diagnose as oracle
from independent_ancestral_scalar_diagnostics_v2 import diagnose, compare


def fixtures():
    rng = np.random.default_rng(20261003)
    for draws in [20, 21, 25, 50, 75, 500, 750]:
        for correlation in [-.95, -.5, 0., .5, .95]:
            values = rng.normal(size=(4, draws))
            for i in range(1, draws):
                values[:, i] += correlation * values[:, i - 1]
            yield f'ar-{draws}-{correlation}', values
        values = rng.normal(size=(4, draws))
        yield f'iid-{draws}', values
        changed = values.copy(); changed[0] += 4
        yield f'location-mismatch-{draws}', changed
        changed = values.copy(); changed[0] *= 15
        yield f'scale-mismatch-{draws}', changed
        yield f'heavy-tail-{draws}', rng.standard_cauchy(size=(4, draws))
        yield f'categorical-ties-{draws}', rng.integers(0, 5, size=(4, draws))
        yield f'binary-ties-{draws}', rng.integers(0, 2, size=(4, draws))
        changed = values.copy(); changed[0] = 1
        yield f'one-constant-chain-{draws}', changed
        yield f'constant-disagreeing-chains-{draws}', np.repeat(np.arange(4.)[:, None], draws, axis=1)
        changed = values.copy(); changed[1, 3] = float('inf')
        yield f'nonfinite-{draws}', changed
        yield f'constant-halves-{draws}', np.tile(np.r_[np.zeros(draws // 2), np.ones(draws - draws // 2)], (4, 1))
    yield 'insufficient', np.zeros((4, 19))
    # Real tied-draw regression retains exact raw input/report provenance.
    regression = Path('results/ancestral/independent-scalar-boundary-regression-20261002-v1/trace.npy')
    yield 'source-tied-quantile-boundary', np.load(regression)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert arviz.__version__ == '0.22.0' and np.__version__ == '2.2.6'
    started = time.perf_counter(); errors = Counter(); states = Counter(); checked = []; numerical_reviews = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        for name, values in fixtures():
            reference, separate = oracle(values), diagnose(values)
            differences, unresolved = compare(reference, separate, values=values)
            if unresolved:
                numerical_reviews.append(dict(fixture=name, unresolved=unresolved,
                    original_rhat=reference['rhat'], independent_rhat=separate['rhat']))
            for key, value in differences.items():
                errors[key] = max(errors[key], value)
            states[reference['status']] += 1
            checked.append(name)
    invalid = []
    for shape in [(3, 30), (5, 30), (120,)]:
        try:
            diagnose(np.ones(shape))
        except ValueError:
            invalid.append(list(shape))
        else:
            raise AssertionError('Invalid chain shape accepted')
    values = np.random.default_rng(77).normal(size=(4, 750))
    reference = oracle(values); altered = dict(reference); altered['bulk_ess'] += 1
    try:
        compare(reference, altered)
    except AssertionError:
        pass
    else:
        raise AssertionError('Changed numeric export accepted')
    paths = [Path(__file__), Path('scripts/independent_ancestral_scalar_diagnostics_v2.py'), Path('scripts/ancestral_chain_diagnostics.py')]
    result = dict(status='passed_independent_ancestral_scalar_direct_lag_contracts',
        fixtures=len(checked), fixture_names=checked, fixture_status_counts=dict(states),
        maximum_absolute_errors=dict(errors), invalid_chain_shapes_rejected=invalid,
        singular_split_variance_fixtures_retained_as_numeric_review=numerical_reviews,
        altered_numeric_export_rejected=True, arviz_oracle_version=arviz.__version__,
        quantile_boundary_regression_sha256=sha('results/ancestral/independent-scalar-boundary-regression-20261002-v1/trace.npy'),
        numpy_version=np.__version__, elapsed_seconds=time.perf_counter() - started,
        source_hashes={str(p): sha(p) for p in paths}, scientific_eligibility=False,
        scope='Software direct-lag versus locked FFT numerical comparison, including tied/binary/odd/antithetical/autocorrelated/shifted/heavy-tailed/constant/nonfinite traces. No production posterior, native-parser validation or biological pilot. The independent implementation imports no ArviZ/project diagnostic estimator.')
    with args.output.open('x') as f:
        f.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
