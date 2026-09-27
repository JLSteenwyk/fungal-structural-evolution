#!/usr/bin/env python3
"""Reconstruct all design diagnostics with SQL moments and eigendecomposition."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--launch', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    launch = json.loads(args.launch.read_text())
    ph, lh = sha(args.plan), sha(args.launch)
    assert launch['plan_sha256'] == ph
    while True:
        try:
            process = psutil.Process(launch['pid'])
            if process.create_time() != launch['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            assert process.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(line.split('=', 1) for line in subprocess.check_output(['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == {'ActiveState': 'inactive', 'ExecMainStatus': '0'}, state
    assert sha(args.plan) == ph and sha(args.launch) == lh
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['plan_sha256'] == ph
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
    actual = pd.read_csv(root / 'design_diagnostics.tsv', sep='\t', dtype={'target_order': str, 'background_order': str})
    covs = ['identity_difference', 'original_coverage_difference', 'log_aligned_length_ratio', 'confidence_fraction_difference']
    setting_keys = ['boundary', 'mask', 'cohort', 'screen', 'target_order', 'background_order']
    keys = ['guide', 'policy', 'scenario_id']
    source = Path(plan['summaries'])
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
    for i in range(4):
        for j in range(i, 4):
            fields.append('covar_pop(' + covs[i] + ',' + covs[j] + ') AS c' + str(i) + str(j))
    checked = 0
    max_singular_error = 0.
    for part in parts:
        path = source / part['path']
        assert sha(path) == part['sha256']
        db.execute('CREATE OR REPLACE TABLE records AS SELECT s.*,v.* EXCLUDE(domain_config_id) FROM selections s JOIN read_parquet(?) v USING(domain_config_id)', [str(path)])
        stats = db.execute('SELECT ' + ','.join(keys + fields) + ' FROM records GROUP BY ' + ','.join(keys)).df().set_index(keys)
        for role in ['target', 'background']:
            reuse = db.execute('SELECT ' + ','.join(keys) + ',max(uses) AS maximum_' + role + '_reuse FROM (SELECT ' + ','.join(keys) + ',' + role + '_id,count(*) uses FROM records GROUP BY ' + ','.join(keys) + ',' + role + '_id) GROUP BY ' + ','.join(keys)).df().set_index(keys)
            stats = stats.join(reuse)
        subset = actual
        for name in setting_keys:
            subset = subset[subset[name].eq(str(part[name]))]
        subset = subset.set_index(keys)
        assert len(subset) == len(stats) == 432 and subset.index.is_unique
        assert set(subset.index) == set(stats.index)
        for key, row in stats.iterrows():
            observed = subset.loc[key]
            for name in ['records'] + list(counts) + ['maximum_target_reuse', 'maximum_background_reuse']:
                assert int(row[name]) == int(observed[name]), (key, name)
            low = np.array([row[c + '_min'] for c in covs])
            high = np.array([row[c + '_max'] for c in covs])
            active = high - low > 1e-12
            for i, cov in enumerate(covs):
                for suffix in ['min', 'max', 'mean']:
                    np.testing.assert_allclose(observed[cov + '_' + suffix], row[cov + '_' + suffix], rtol=1e-10, atol=1e-12)
                assert int(observed[cov + '_constant']) == int(not active[i])
            matrix = np.array([[row['c' + str(min(i,j)) + str(max(i,j))] for j in range(4)] for i in range(4)])[np.ix_(active, active)]
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
    assert checked == len(actual) == receipt['strata'] == 82944
    assert int(actual.full_rank_after_constant_removal.eq(0).sum()) == receipt['rank_deficient_after_constant_removal']
    for cov in covs:
        assert int(actual[cov + '_constant'].sum()) == receipt['constant_covariate_cells'][cov]
    assert int(actual.zero_in_each_marginal_range.eq(0).sum()) == receipt['strata_with_zero_outside_some_marginal_range']
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    result = dict(status='passed_full_matched_design_diagnostics_readback', strata=checked,
                  maximum_singular_value_error=max_singular_error,
                  source_receipt_sha256=sha(root / 'receipt.json'), script_sha256=sha(__file__),
                  scope='Every count, covariate range/mean/constant flag, reuse maximum, scaled singular spectrum/rank/condition and marginal-zero flag independently reconstructed using SQL aggregation and covariance eigenvalues. No inferential model fitted.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
