"""Read back completed reference and compare all 15 pinned real-input cases."""
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json
from matched_covariance_information_fast import covariance_information


def main():
    threadpool_limits(1)
    plan_path = Path('metadata/matched_information_timing_plan_20260928.json')
    plan = json.loads(plan_path.read_text())
    launch = json.loads(Path('metadata/matched_information_timing_launch_20260928.json').read_text())
    state = dict(line.split('=', 1) for line in subprocess.check_output([
        'systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'Result',
        '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0'), state
    root = Path(plan['output'])
    source = root / 'receipt.json'
    receipt = json.loads(source.read_text())
    assert receipt['status'] == 'completed_selected_design_information_evaluations'
    assert receipt['plan_sha256'] == sha(plan_path) == launch['plan_sha256']
    pins = dict(plan['pins'])
    pins.update({str(plan_path): sha(plan_path), str(source): sha(source),
                 'scripts/matched_covariance_information_fast.py': sha('scripts/matched_covariance_information_fast.py'),
                 str(Path(__file__)): sha(__file__)})

    def verify():
        for path, digest in pins.items():
            assert sha(path) == digest, path
        for name, digest in receipt['artifacts'].items():
            assert sha(root / name) == digest, name

    verify()
    assert len(receipt['jobs']) == len(plan['jobs']) == 15
    assert [r['job'] for r in receipt['jobs']] == plan['jobs']
    out = Path('results/model_validation/matched-information-shortcut-comparison-20260928-v1')
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for row in receipt['jobs']:
        job = row['job']
        path = root / row['artifact']
        assert sha(path) == row['sha256']
        with np.load(path, allow_pickle=False) as h:
            original = {key: h[key] for key in h.files}
        saved = json.loads(path.with_suffix('.json').read_text())
        assert saved == dict(plan_sha256=sha(plan_path), **row)
        with np.load(job['cache'], allow_pickle=False) as h:
            a = {key: h[key] for key in h.files}
        with np.load(job['factor'], allow_pickle=False) as h:
            factor = h['factor'][a['pattern_rows']]
        old = json.loads(Path(job['original_fit']).read_text())['payload']
        x = np.column_stack([np.ones(len(a['matrix'])),
                             a['matrix'][:, 1:][:, a['active_covariates']] / a['covariate_scales']])
        assert row['records'] == len(x)
        v = old['profiled_scale'] * np.array([1., *old['ratios']])
        np.testing.assert_array_equal(original['absolute_variances'], v)
        np.testing.assert_array_equal(original['exact_zero_components'], v == 0)
        info = original['expected_reml_information']
        assert info.shape == (4, 4) and np.isfinite(info).all()
        np.testing.assert_allclose(info, info.T, rtol=1e-12, atol=1e-12)
        scale = np.sqrt(np.diag(info))
        eig = np.linalg.eigvalsh(info / np.outer(scale, scale))
        np.testing.assert_allclose(eig, row['normalized_information_eigenvalues'], rtol=1e-10, atol=1e-12)
        assert int(np.sum(eig > row['rank_tolerance'])) == row['numerical_information_rank']
        started = time.perf_counter()
        fast = covariance_information(a['background'], a['family'], factor, x, v, plan['block_size'])
        elapsed = time.perf_counter() - started
        errors = {}
        for key, expected in original.items():
            np.testing.assert_allclose(fast[key], expected, rtol=1e-7, atol=1e-8, err_msg=key)
            if expected.dtype.kind != 'b':
                errors[key] = float(np.max(abs(fast[key]-expected)))
        assert fast['numerical_information_rank'] == row['numerical_information_rank']
        target = out / path.name
        np.savez_compressed(target, **{key: fast[key] for key in original})
        rows.append(dict(job=job, records=len(x), reference_seconds=row['elapsed_seconds'],
                         shortcut_seconds=elapsed, evaluation_method=fast['evaluation_method'],
                         numerical_information_rank=fast['numerical_information_rank'],
                         maximum_absolute_errors=errors, artifact=target.name, sha256=sha(target)))
        print(len(rows), '/15', len(x), 'records', round(elapsed, 3), 'seconds; agreement passed', flush=True)
    verify()
    result = dict(status='passed_all_15_real_input_information_shortcut_comparisons',
        reference_terminal_state=state, source_pins=pins, jobs=rows,
        reference_seconds_sum=sum(r['reference_seconds'] for r in rows),
        shortcut_seconds_sum=sum(r['shortcut_seconds'] for r in rows),
        artifacts={p.name: sha(p) for p in out.iterdir()},
        scope='Full readback of selected reference artifacts and real-input shortcut comparison. '
              'Single timings on availability/size-selected inputs, not full-grid runtime estimates. '
              'No calibrated intervals, boundary regularity or scientific inference.')
    write_json(out / 'receipt.json', result)
    write_json(Path('metadata/matched_information_shortcut_comparison_20260928.json'), result)


if __name__ == '__main__':
    main()
