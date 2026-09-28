"""Evaluate exact information at existing fitted values; no interval claims."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json
from audit_matched_simulation_cache import replay
from matched_covariance_information import covariance_information


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    ph = sha(args.plan)
    threadpool_limits(1)

    def verify():
        assert sha(args.plan) == ph
        for path, digest in plan['pins'].items():
            assert sha(path) == digest, path

    verify()
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for job in plan['jobs']:
        with np.load(job['cache'], allow_pickle=False) as h:
            a = {k: h[k] for k in h.files}
        with np.load(job['factor'], allow_pickle=False) as h:
            factor = h['factor']
        old = json.loads(Path(job['original_fit']).read_text())['payload']
        readback = replay(a, old, factor)
        x = np.column_stack([np.ones(len(a['matrix'])),
                             a['matrix'][:, 1:][:, a['active_covariates']] / a['covariate_scales']])
        factor = factor[a['pattern_rows']]
        started = time.perf_counter()
        result = covariance_information(a['background'], a['family'], factor, x,
                                       old['profiled_scale'] * np.array([1., *old['ratios']]),
                                       block_size=plan['block_size'])
        elapsed = time.perf_counter() - started
        np.testing.assert_allclose(result['conditional_beta_covariance'],
                                   old['conditional_beta_covariance'], rtol=1e-7, atol=1e-8)
        arrays = {k: v for k, v in result.items() if isinstance(v, np.ndarray)}
        path = out / (job['fit_input_id'] + '-' + job['tree'] + '.npz')
        np.savez_compressed(path, **arrays)
        with np.load(path, allow_pickle=False) as saved:
            assert set(saved.files) == set(arrays)
            for key, value in arrays.items():
                np.testing.assert_array_equal(saved[key], value)
        row = dict(job=job, records=len(x), elapsed_seconds=elapsed,
                   numerical_information_rank=result['numerical_information_rank'],
                   rank_tolerance=result['rank_tolerance'],
                   normalized_information_eigenvalues=result['normalized_information_eigenvalues'].tolist(),
                   exact_zero_components=result['exact_zero_components'].tolist(),
                   original_replay=readback, artifact=path.name, sha256=sha(path))
        rows.append(row)
        write_json(out / (path.stem + '.json'), dict(plan_sha256=ph, **row))
        print(len(rows), '/', len(plan['jobs']), len(x), 'records', round(elapsed, 3),
              'seconds; information rank', row['numerical_information_rank'], flush=True)
    verify()
    write_json(out / 'receipt.json', dict(status='completed_selected_design_information_evaluations',
        plan_sha256=ph, jobs=rows, artifacts={p.name: sha(p) for p in out.iterdir()},
        scope='All 15 previously selected size/tree timing designs at original fitted values. '
              'Numerical rank is a diagnostic, not a regularity or interval-coverage guarantee. '
              'Availability/size-selected cases cannot estimate full-grid prevalence or runtime.'))


if __name__ == '__main__':
    main()
