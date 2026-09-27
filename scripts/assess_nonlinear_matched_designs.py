#!/usr/bin/env python3
"""Assess quadratic and cubic matched-record designs; preserve all settings."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import qr
from screen_duplication_domain_alignment_coverage import sha

COVARIATES = ['identity_difference', 'original_coverage_difference',
              'log_aligned_length_ratio', 'confidence_fraction_difference',
              'identity_power_2_difference', 'identity_power_3_difference']
KEYS = ['guide', 'policy', 'scenario_id']


def assess(x):
    """Rank after explicit constant removal and population-SD scaling.

    Intercept adds one to centered rank. Constants remain reported, never silently
    fit; absolute range <=1e-12 is treated as numerically constant.
    """
    assert len(x) and np.isfinite(x).all()
    low, high = x.min(axis=0), x.max(axis=0)
    active = high - low > 1e-12
    centered = x[:, active] - x[:, active].mean(axis=0)
    scales = np.sqrt(np.mean(centered ** 2, axis=0))
    singular = np.linalg.svd(centered / scales / np.sqrt(len(x)), compute_uv=False)
    if active.any():
        _, triangular, _ = qr(centered / scales / np.sqrt(len(x)), mode='economic', pivoting=True)
        independent = np.linalg.svd(triangular, compute_uv=False)
        np.testing.assert_allclose(singular, independent, atol=1e-12, rtol=1e-10)
        assert int(np.sum(independent > 1e-7)) == int(np.sum(singular > 1e-7))
    rank = int(np.sum(singular > 1e-7))
    result = dict(records=len(x), active_covariates=int(active.sum()),
                  centered_rank=rank, intercept_design_rank=rank + 1,
                  full_rank_after_constant_removal=int(rank == active.sum()),
                  scaled_singular_values=json.dumps(singular.tolist(), separators=(',', ':')),
                  scaled_condition=(float(singular[0] / singular[-1]) if rank == active.sum() and len(singular) else ''),
                  zero_in_each_marginal_range=int(np.all((low <= 0) & (high >= 0))))
    for i, name in enumerate(COVARIATES[:x.shape[1]]):
        result[name + '_min'] = low[i]
        result[name + '_max'] = high[i]
        result[name + '_mean'] = x[:, i].mean()
        result[name + '_constant'] = int(not active[i])
    return result


def fixtures():
    a = np.arange(20.) - 9.5
    x = np.column_stack([a, a * 2, np.ones(20), a ** 2])
    r = assess(x)
    assert r['active_covariates'] == 3 and r['centered_rank'] == 2
    assert r['log_aligned_length_ratio_constant'] == 1
    z = assess(np.zeros((10, 4)))
    assert z['intercept_design_rank'] == 1 and z['active_covariates'] == 0
    rng = np.random.default_rng(9371)
    assert assess(rng.normal(size=(100, 4)))['centered_rank'] == 4
    expanded = rng.normal(size=(100, 6))
    assert assess(expanded)['centered_rank'] == 6
    expanded[:, 5] = 2 * expanded[:, 4]
    assert assess(expanded)['centered_rank'] == 5


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)
    def verify():
        assert sha(args.plan) == plan_hash
        for path, digest in plan['pins'].items():
            assert sha(path) == digest, path
    verify()
    fixtures()
    source = Path(plan['summaries'])
    nonlinear = Path(plan['nonlinear'])
    nr = json.loads((nonlinear / 'receipt.json').read_text())
    assert nr['status'] == 'complete_polynomial_identity_contrasts_full_dual_implementation_check'
    for name, digest in nr['artifacts'].items():
        assert sha(nonlinear / name) == digest
    partitions = json.loads((source / 'partition_manifest.json').read_text())
    assert len(partitions) == 192
    selections = pd.read_csv(plan['selections'], sep='\t', usecols=['target_id', 'background_id', 'policy', 'scenario_id', 'domain_config_id'])
    assert len(selections) == 2786912
    nodes = pd.read_csv(plan['nodes'], sep='\t')
    target = nodes[nodes.role.eq('target')][['node_id', 'guide', 'family', 'family_component', 'taxon_a']].rename(columns={'node_id': 'target_id', 'taxon_a': 'taxon_id'})
    selected = selections.merge(target, on='target_id', validate='many_to_one')
    assert len(selected) == len(selections)
    expected = pd.read_csv(source / 'record_summary.tsv', sep='\t', dtype={'target_order': str, 'background_order': str})
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    rows = []
    setting_keys = ['boundary', 'mask', 'cohort', 'screen', 'target_order', 'background_order']
    for part in partitions:
        path = source / part['path']
        assert sha(path) == part['sha256']
        values = pd.read_parquet(path, columns=['domain_config_id'] + COVARIATES[:4])
        extra = pd.read_parquet(nonlinear / f"{part['index']:03d}.parquet")
        assert set(values.domain_config_id) == set(extra.domain_config_id)
        values = values.merge(extra[['domain_config_id'] + COVARIATES[4:]], on='domain_config_id', validate='one_to_one')
        records = selected.merge(values, on='domain_config_id', validate='many_to_one')
        grid = expected
        for key in setting_keys:
            grid = grid[grid[key].eq(str(part[key]))]
        assert len(grid) == 432
        counts = grid.set_index(KEYS).matched_records.to_dict()
        seen = set()
        for key, frame in records.groupby(KEYS, sort=True):
            assert len(frame) == counts[key]
            seen.add(key)
            row = dict(zip(KEYS, key))
            row.update({k: part[k] for k in setting_keys})
            # Reuse counts only; assess each nested fixed-effect space separately.
            for name, column in [('targets', 'target_id'), ('backgrounds', 'background_id'),
                                 ('families', 'family'), ('family_components', 'family_component'), ('taxa', 'taxon_id')]:
                row[name] = int(frame[column].nunique())
            row['maximum_target_reuse'] = int(frame.target_id.value_counts().max())
            row['maximum_background_reuse'] = int(frame.background_id.value_counts().max())
            for degree in [2, 3]:
                model_row = dict(row, polynomial_degree=degree)
                model_row.update(assess(frame[COVARIATES[:degree+3]].to_numpy()))
                rows.append(model_row)
        # The existing audited grid has no empty groups; fail explicitly if that changes.
        assert seen == set(counts) and all(n > 0 for n in counts.values())
        print('Assessed full covariate design', part['index'] + 1, '/ 192', flush=True)
    assert len(rows) == 165888
    table = pd.DataFrame(rows)
    table.to_csv(out / 'design_diagnostics.tsv', sep='\t', index=False)
    verify()
    result = dict(status='complete_nonlinear_design_diagnostics_pending_independent_readback',
                  plan_sha256=plan_hash, settings=len(partitions), strata=len(table),
                  rank_deficient_after_constant_removal=int((table.full_rank_after_constant_removal == 0).sum()),
                  constant_covariate_cells={name: int(table[name + '_constant'].fillna(0).sum()) for name in COVARIATES},
                  strata_with_zero_outside_some_marginal_range=int((table.zero_in_each_marginal_range == 0).sum()),
                  rank_absolute_tolerance=1e-7, constant_absolute_range_tolerance=1e-12,
                  artifacts={'design_diagnostics.tsv': sha(out / 'design_diagnostics.tsv')},
                  scope='All 192 settings and 432 strata each, for both quadratic and cubic specifications; absent cubic statistics in quadratic rows remain blank. Singular values checked by pivoted QR plus SVD. Direct SVD of centered population-SD-scaled equal-record covariates. Intercept rank adds one. Positive alternative weights preserve exact rank but can alter conditioning. Marginal zero coverage is not joint convex-hull support. No response fit, uncertainty estimate, effective sample size, phylogenetic adjustment or causal claim.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
