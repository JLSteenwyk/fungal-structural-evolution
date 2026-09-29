"""Describe the complete audited residual grid without adequacy thresholds.

Each fit has equal weight within its tree. Fits overlap in observations and
settings; empirical ranges across fits are neither confidence intervals nor
sampling distributions. Failed dispositions remain in denominator tables.
"""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha, write_json

GOOD = 'descriptive_marginal_residual_diagnostics'
METRICS = ['mean', 'second_raw_moment', 'third_raw_moment',
           'fourth_raw_moment', 'fraction_absolute_above_2',
           'fraction_absolute_above_3', 'maximum_fixed_effect_variance_fraction']


def summarize(frame):
    assert not frame.duplicated(['fit_input_id', 'tree']).any()
    rows = frame.copy()
    valid = rows.status.eq(GOOD)
    n = rows.loc[valid, 'records'].to_numpy(float)
    assert np.isfinite(n).all() and (n > 0).all() and (n == np.floor(n)).all()
    for threshold in [2, 3]:
        counts = rows.loc[valid, f'absolute_above_{threshold}'].to_numpy(float)
        assert np.isfinite(counts).all() and (counts >= 0).all()
        assert (counts <= n).all() and (counts == np.floor(counts)).all()
        rows.loc[valid, f'fraction_absolute_above_{threshold}'] = counts / n
    assert (rows.loc[valid, 'absolute_above_3'] <= rows.loc[valid, 'absolute_above_2']).all()
    counts = rows.groupby(['tree', 'status'], dropna=False).size().rename('fits').reset_index()
    counts['tree_total_fits'] = counts.groupby('tree')['fits'].transform('sum')
    summaries = []
    for tree, group in rows.groupby('tree', sort=True):
        accepted = group[group.status.eq(GOOD)]
        for metric in METRICS:
            values = accepted[metric].to_numpy(float)
            assert np.isfinite(values).all(), (tree, metric)
            q = np.quantile(values, [0, .025, .25, .5, .75, .975, 1]) if len(values) else [None]*7
            summaries.append(dict(tree=tree, metric=metric, total_fits=len(group),
                                  available_fits=len(values), unavailable_fits=len(group)-len(values),
                                  **dict(zip(['minimum','q025','q25','median','q75','q975','maximum'], q))))
    return rows, counts, pd.DataFrame(summaries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-plan', type=Path, required=True)
    parser.add_argument('--audit-launch', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.audit_plan.read_text())
    launch = json.loads(args.audit_launch.read_text())
    assert sha(args.audit_plan) == launch['plan_sha256']
    state = dict(line.split('=', 1) for line in subprocess.check_output(
        ['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState',
         '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0'), state
    root = Path(plan['output'])
    receipt_path = root/'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'complete_all_selected_residual_diagnostic_replays'
    assert receipt['plan_sha256'] == sha(args.audit_plan)
    assert receipt['inputs'] == 28808 and receipt['fits'] == 144040
    bindings = {str(args.audit_plan): sha(args.audit_plan),
                str(args.audit_launch): sha(args.audit_launch),
                str(receipt_path): sha(receipt_path), **plan['pins']}
    bindings.update({str(root/name): h for name, h in receipt['artifacts'].items()})
    def verify():
        for name, digest in bindings.items():
            assert sha(name) == digest, name
    verify()
    frame = pd.read_parquet(root/'residual_fit_summaries.parquet')
    assert len(frame) == 144040 and frame.fit_input_id.nunique() == 28808
    assert frame.tree.nunique() == 5 and frame.groupby('tree').size().eq(28808).all()
    assert frame.groupby('fit_input_id').tree.nunique().eq(5).all()
    assert frame.status.value_counts().to_dict() == receipt['status_counts']
    rows, counts, summary = summarize(frame)
    args.output.mkdir(parents=True, exist_ok=False)
    rows.to_parquet(args.output/'all_fit_diagnostics.parquet', index=False)
    counts.to_csv(args.output/'status_counts.tsv', sep='\t', index=False)
    summary.to_csv(args.output/'descriptive_ranges.tsv', sep='\t', index=False)
    # Check the complete serialized tables, not just the in-memory summaries.
    pd.testing.assert_frame_equal(pd.read_parquet(args.output/'all_fit_diagnostics.parquet'), rows)
    pd.testing.assert_frame_equal(pd.read_csv(args.output/'status_counts.tsv', sep='\t'), counts)
    pd.testing.assert_frame_equal(pd.read_csv(args.output/'descriptive_ranges.tsv', sep='\t'), summary,
                                  check_dtype=False, rtol=1e-12, atol=1e-12)
    verify()
    write_json(args.output/'receipt.json', dict(
        status='complete_descriptive_audited_residual_grid_summary', fits=len(rows),
        inputs=28808, trees=5, summary_rows=len(summary), source_bindings=bindings,
        script_sha256=sha(__file__), audit_terminal_state=state,
        artifacts={p.name: sha(p) for p in args.output.iterdir()},
        scope='Equal-weight per-fit descriptive ranges, with all dispositions retained. '
              'Overlapping fits are not independent replicates. Ranges are not confidence '
              'intervals, calibrated adequacy thresholds, or evidence of biological effects.'))


if __name__ == '__main__':
    main()
