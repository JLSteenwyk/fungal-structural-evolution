#!/usr/bin/env python3
"""Time complete synthetic ML fits at observed record/group/factor dimensions."""
import json
import time
from pathlib import Path
import numpy as np
from cached_matched_ml import profiled_ml
from matched_mixed_covariance import MatchedCovariance
from fit_matched_ml import fit_variance_ratios
from screen_duplication_domain_alignment_coverage import sha


def main():
    paths = [__file__, 'scripts/cached_matched_ml.py', 'scripts/matched_ml_gradient.py',
             'scripts/matched_mixed_covariance.py', 'scripts/fit_matched_ml.py',
             'metadata/matched_ml_optimizer_checks_20260927.json']
    pins = {p:sha(p) for p in paths}
    out = Path('results/model_validation/matched-ml-size-benchmark-20260927-v1')
    out.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(913264)
    rows = []
    for n, groups, families in [(148, 120, 44), (2320, 1500, 400), (10963, 6478, 1634)]:
        bg = np.concatenate([np.arange(groups), rng.integers(0, groups, n-groups)])
        membership = np.concatenate([np.arange(families), rng.integers(0, families, groups-families)])
        family = membership[bg]
        factor = rng.normal(size=(n, 242))/np.sqrt(242)
        a, b = rng.uniform(size=(2, n))
        covariates = np.column_stack([a-b, rng.normal(size=(n, 3)), a*a-b*b, a**3-b**3])
        covariates /= covariates.std(axis=0)
        design = np.column_stack([np.ones(n), covariates])
        y = design@np.array([.1, .3, -.1, .2, .1, -.2, .1]) + rng.normal(size=n)
        y += .5*rng.normal(size=groups)[bg] + .7*rng.normal(size=families)[family]
        y += factor@rng.normal(size=242)
        print(f'Starting complete 22-attempt synthetic fit: {n} records', flush=True)
        start = time.perf_counter()
        fit = fit_variance_ratios(bg, family, factor, design, y)
        fitting_seconds = time.perf_counter()-start
        assert len(fit['candidates']) == 22
        before = time.perf_counter()
        maximum = 0.
        for candidate in fit['candidates']:
            ratios = np.expm1(candidate['theta'])
            direct = profiled_ml(MatchedCovariance(bg, family, factor, 1., *ratios), design, y)
            np.testing.assert_allclose(candidate['objective'], direct['negative_profiled_ml'], rtol=1e-8, atol=1e-7)
            maximum = max(maximum, abs(candidate['objective']-direct['negative_profiled_ml']))
        validation_seconds = time.perf_counter()-before
        path = out/f'fit-{n}.json'
        path.write_text(json.dumps(fit, indent=2)+'\n')
        row = dict(records=n, backgrounds=groups, family_components=families, species_factor_columns=242,
                   fixed_columns=7, fit_seconds=fitting_seconds, candidate_readback_seconds=validation_seconds,
                   candidate_checks=22, maximum_direct_objective_error=maximum,
                   status=fit['status'], evaluations=sum(c['evaluations'] for c in fit['candidates']),
                   maximum_projected_gradient=max(map(abs, fit['projected_gradient'])))
        rows.append(row)
        print(json.dumps(row), flush=True)
    for path, digest in pins.items():
        assert sha(path) == digest, path
    receipt = dict(status='complete_synthetic_ml_size_benchmark', runs=rows, pins=pins,
                   artifacts={p.name:sha(p) for p in out.iterdir()},
                   scope='Three synthetic timing cases at observed dimensions; all 66 candidate objectives checked using direct residual evaluation. Optimization flags preserved. Timings are not biological results, a representative convergence distribution, a full-grid ETA or a substitute for full-data fitting. Existing production fits unchanged.')
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    main()
