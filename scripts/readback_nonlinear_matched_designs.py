#!/usr/bin/env python3
"""Reconstruct quadratic and cubic designs independently with SQL moments."""
import argparse
import json
import subprocess
import time
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
import psutil
from screen_duplication_domain_alignment_coverage import sha


def numerical_fixtures():
    """Exercise SQL covariance versus direct spectra, including polynomial terms."""
    rng = np.random.default_rng(291704)
    a, b = rng.uniform(size=(2, 500))
    matrices = [np.column_stack([a-b, rng.normal(size=(500, 3)), a*a-b*b, a**3-b**3])]
    duplicate = matrices[0].copy()
    duplicate[:, -1] = duplicate[:, -2]
    matrices.append(duplicate)
    db = duckdb.connect()
    db.execute('SET threads=1')
    for x, expected_rank in zip(matrices, [6, 5]):
        db.register('fixture', pd.DataFrame(x, columns=[f'x{i}' for i in range(6)]))
        values = db.execute('SELECT ' + ','.join(f'covar_pop(x{i},x{j})' for i in range(6) for j in range(6)) + ' FROM fixture').fetchone()
        covariance = np.array(values).reshape(6, 6)
        sd = np.sqrt(np.diag(covariance))
        eigen = np.linalg.eigvalsh(covariance / np.outer(sd, sd))[::-1]
        direct = np.linalg.svd((x-x.mean(axis=0))/x.std(axis=0)/np.sqrt(len(x)), compute_uv=False)
        assert int((eigen > 1e-14).sum()) == expected_rank
        np.testing.assert_allclose(np.sqrt(np.maximum(eigen, 0)), direct, rtol=1e-7, atol=3e-8)
    db.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--launch', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    numerical_fixtures()
    plan = json.loads(args.plan.read_text())
    launch = json.loads(args.launch.read_text())
    ph, lh = sha(args.plan), sha(args.launch)
    assert launch['plan_sha256'] == ph
    while True:
        try:
            process = psutil.Process(launch['pid'])
            if process.create_time() != launch['create_time'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            assert process.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(line.split('=', 1) for line in subprocess.check_output(['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'ExecMainStatus', '-p', 'Result'], text=True).splitlines())
    assert state == {'ActiveState': 'inactive', 'ExecMainStatus': '0', 'Result': 'success'}, state
    assert sha(args.plan) == ph and sha(args.launch) == lh
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['plan_sha256'] == ph
    assert receipt['status'] == 'complete_nonlinear_design_diagnostics_pending_independent_readback'
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
    actual = pd.read_csv(root / 'design_diagnostics.tsv', sep='\t', dtype={'target_order': str, 'background_order': str})
    for suffix in ['min', 'max', 'mean', 'constant']:
        assert actual.loc[actual.polynomial_degree.eq(2), 'identity_power_3_difference_' + suffix].isna().all()
    covs = ['identity_difference', 'original_coverage_difference', 'log_aligned_length_ratio', 'confidence_fraction_difference', 'identity_power_2_difference', 'identity_power_3_difference']
    setting_keys = ['boundary', 'mask', 'cohort', 'screen', 'target_order', 'background_order']
    keys = ['guide', 'policy', 'scenario_id']
    source = Path(plan['summaries'])
    nonlinear = Path(plan['nonlinear'])
    nr = json.loads((nonlinear / 'receipt.json').read_text())
    assert nr['status'] == 'complete_polynomial_identity_contrasts_full_dual_implementation_check'
    for name, digest in nr['artifacts'].items():
        assert sha(nonlinear / name) == digest
    parts = json.loads((source / 'partition_manifest.json').read_text())
    db = duckdb.connect()
    db.execute('SET threads=1')
    db.execute("SET memory_limit='8GB'")
    db.execute("CREATE TABLE nodes AS SELECT * FROM read_csv(?,delim='\t',header=true,all_varchar=true) WHERE role='target'", [plan['nodes']])
    db.execute("CREATE TABLE selections AS SELECT s.*,n.guide,n.family,n.family_component,n.taxon_a FROM read_csv(?,delim='\t',header=true,all_varchar=true) s JOIN nodes n ON s.target_id=n.node_id", [plan['selections']])
    assert db.execute('SELECT count(*) FROM selections').fetchone()[0] == 2786912
    counts = {'targets': 'target_id', 'backgrounds': 'background_id', 'families': 'family', 'family_components': 'family_component', 'taxa': 'taxon_a'}
    fields = ['count(*) AS records'] + ['count(DISTINCT ' + column + ') AS ' + label for label, column in counts.items()]
    for cov in covs:
        fields += [fun + '(' + cov + ') AS ' + cov + '_' + label for fun, label in [('min', 'min'), ('max', 'max'), ('avg', 'mean')]]
    for i in range(6):
        for j in range(i, 6):
            fields.append('covar_pop(' + covs[i] + ',' + covs[j] + ') AS c' + str(i) + str(j))
    checked = 0
    max_singular_error = 0.
    for part in parts:
        path = source / part['path']
        assert sha(path) == part['sha256']
        db.execute('CREATE OR REPLACE TABLE records AS SELECT s.*,v.* EXCLUDE(domain_config_id),n.identity_power_2_difference,n.identity_power_3_difference FROM selections s JOIN read_parquet(?) v USING(domain_config_id) JOIN read_parquet(?) n USING(domain_config_id)', [str(path), str(nonlinear / f"{part['index']:03d}.parquet")])
        stats = db.execute('SELECT ' + ','.join(keys + fields) + ' FROM records GROUP BY ' + ','.join(keys)).df().set_index(keys)
        for role in ['target', 'background']:
            reuse = db.execute('SELECT ' + ','.join(keys) + ',max(uses) AS maximum_' + role + '_reuse FROM (SELECT ' + ','.join(keys) + ',' + role + '_id,count(*) uses FROM records GROUP BY ' + ','.join(keys) + ',' + role + '_id) GROUP BY ' + ','.join(keys)).df().set_index(keys)
            stats = stats.join(reuse)
        subset = actual
        for name in setting_keys:
            subset = subset[subset[name].eq(str(part[name]))]
        for degree in [2, 3]:
            degree_subset = subset[subset.polynomial_degree.eq(degree)].set_index(keys)
            assert len(degree_subset) == len(stats) == 432 and degree_subset.index.is_unique
            assert set(degree_subset.index) == set(stats.index)
            for key, row in stats.iterrows():
                observed = degree_subset.loc[key]
                for name in ['records'] + list(counts) + ['maximum_target_reuse', 'maximum_background_reuse']:
                    assert int(row[name]) == int(observed[name]), (key, name)
                low = np.array([row[c + '_min'] for c in covs[:degree+3]])
                high = np.array([row[c + '_max'] for c in covs[:degree+3]])
                active = high - low > 1e-12
                for i, cov in enumerate(covs[:degree+3]):
                    for suffix in ['min', 'max', 'mean']:
                        np.testing.assert_allclose(observed[cov + '_' + suffix], row[cov + '_' + suffix], rtol=1e-10, atol=1e-12)
                    assert int(observed[cov + '_constant']) == int(not active[i])
                matrix = np.array([[row['c' + str(min(i,j)) + str(max(i,j))] for j in range(degree+3)] for i in range(degree+3)])[np.ix_(active, active)]
                sd = np.sqrt(np.diag(matrix))
                correlation = matrix / np.outer(sd, sd)
                eigen = np.linalg.eigvalsh(correlation)[::-1]
                assert (eigen >= -1e-10).all()
                singular = np.sqrt(np.maximum(eigen, 0))
                saved = np.array(json.loads(observed.scaled_singular_values))
                np.testing.assert_allclose(singular, saved, rtol=1e-7, atol=3e-8)
                if len(singular):
                    max_singular_error = max(max_singular_error, float(np.max(abs(singular - saved))))
                rank = int((eigen > 1e-14).sum())
                assert rank == int(observed.centered_rank)
                assert rank + 1 == int(observed.intercept_design_rank)
                assert int(active.sum()) == int(observed.active_covariates)
                assert int(rank == active.sum()) == int(observed.full_rank_after_constant_removal)
                if rank == active.sum() and len(singular):
                    np.testing.assert_allclose(observed.scaled_condition, singular[0] / singular[-1], rtol=1e-6, atol=1e-8)
                else:
                    assert pd.isna(observed.scaled_condition)
                assert int(np.all((low <= 0) & (high >= 0))) == int(observed.zero_in_each_marginal_range)
                checked += 1
        print('Independently checked design setting', part['index'] + 1, '/ 192', flush=True)
    assert checked == len(actual) == receipt['strata'] == 165888
    assert int(actual.full_rank_after_constant_removal.eq(0).sum()) == receipt['rank_deficient_after_constant_removal']
    for cov in covs:
        assert int(actual[cov + '_constant'].fillna(0).sum()) == receipt['constant_covariate_cells'][cov]
    assert int(actual.zero_in_each_marginal_range.eq(0).sum()) == receipt['strata_with_zero_outside_some_marginal_range']
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    result = dict(status='passed_full_nonlinear_matched_design_diagnostics_readback', strata=checked,
                  maximum_singular_value_error=max_singular_error,
                  source_receipt_sha256=sha(root / 'receipt.json'), script_sha256=sha(__file__),
                  scope='Both polynomial degrees, all 192 settings and 432 strata. Every count, covariate range/mean/constant flag, reuse maximum, scaled singular spectrum/rank/condition and marginal-zero flag independently reconstructed using SQL aggregation and covariance eigenvalues. No inferential model fitted.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
