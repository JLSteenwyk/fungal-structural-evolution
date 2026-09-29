"""Plot descriptive ranges of the complete audited residual grid."""
import argparse
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--completion', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    evidence = json.loads(args.completion.read_text())
    assert evidence['status'] == 'complete_descriptive_summary_with_full_table_readback'
    receipt_path = Path(evidence['receipt_path'])
    assert sha(receipt_path) == evidence['receipt_sha256']
    receipt = json.loads(receipt_path.read_text())
    assert receipt == evidence['receipt']
    source = receipt_path.parent/'descriptive_ranges.tsv'
    assert sha(source) == receipt['artifacts'][source.name]
    table = pd.read_csv(source, sep='\t')
    assert len(table) == 35 and not table.duplicated(['tree', 'metric']).any()
    assert table.total_fits.eq(28808).all() and table.available_fits.eq(28808).all()
    assert table.unavailable_fits.eq(0).all()
    quantiles = ['minimum', 'q025', 'q25', 'median', 'q75', 'q975', 'maximum']
    assert np.isfinite(table[quantiles]).all().all()
    assert (np.diff(table[quantiles].to_numpy(), axis=1) >= 0).all()
    trees = ['mafft_guide', 'profile_guide', 'pmsf_mafft_profile',
             'pmsf_profile_mafft', 'pmsf_profile_profile']
    assert set(table.tree) == set(trees)
    labels = ['MAFFT guide', 'Profile guide', 'PMSF: MAFFT / profile',
              'PMSF: profile / MAFFT', 'PMSF: profile / profile']
    metrics = [('second_raw_moment', 'Second raw moment', 1),
               ('fourth_raw_moment', 'Fourth raw moment', 1),
               ('fraction_absolute_above_3', 'Residuals with |z| > 3 (%)', 100)]
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(source, args.output/source.name)
    fig, axes = plt.subplots(1, 3, figsize=(12, 5.5), sharey=True)
    plotted = []
    for ax, (metric, title, scale) in zip(axes, metrics):
        subset = table[table.metric.eq(metric)].set_index('tree').loc[trees]
        for y, (tree, row) in enumerate(subset.iterrows()):
            q = row[quantiles].to_numpy(float)*scale
            ax.plot([q[0], q[6]], [y, y], color='#c8c8c8', lw=1, zorder=1)
            ax.plot([q[1], q[5]], [y, y], color='#267b91', lw=2, zorder=2)
            ax.plot([q[2], q[4]], [y, y], color='#267b91', lw=7, solid_capstyle='butt', zorder=3)
            ax.plot(q[3], y, 'o', color='#172e38', ms=5, zorder=4)
            plotted.append(dict(tree=tree, metric=metric, display_scale=scale,
                                **dict(zip(quantiles, q.tolist()))))
        low = float(subset.minimum.min()*scale)
        high = float(subset.maximum.max()*scale)
        margin = (high-low)*.08
        ax.set_xlim(max(0, low-margin), high+margin)
        ax.set_title(title, fontsize=11)
        ax.set_yticks(range(5), labels)
        ax.set_ylim(4.6, -.6)
        ax.grid(axis='x', color='#eeeeee')
        ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
        assert ax.get_xlim()[0] <= low and ax.get_xlim()[1] >= high
    fig.suptitle('Residual diagnostics across five tree settings', fontsize=15)
    handles = [Line2D([0], [0], marker='o', color='#172e38', lw=0, label='Median'),
               Line2D([0], [0], color='#267b91', lw=7, label='25th–75th percentiles'),
               Line2D([0], [0], color='#267b91', lw=2, label='2.5th–97.5th percentiles'),
               Line2D([0], [0], color='#c8c8c8', lw=1, label='Minimum–maximum')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .13), ncol=4,
               frameon=False, fontsize=9)
    fig.text(.5, .045, '28,808 overlapping fits per tree; each fit weighted equally. No unavailable diagnostic fits.\n'
             'Ranges describe variation across fits, not confidence intervals or calibrated rejection thresholds.\n'
             'Marginal residuals use fitted covariance; correlations remain.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .23, 1, .94))
    for extension in ['png', 'pdf']:
        fig.savefig(args.output/f'residual_diagnostics.{extension}', dpi=180)
    plt.close(fig)
    pd.DataFrame(plotted).to_csv(args.output/'plotted_ranges.tsv', sep='\t', index=False)
    assert sha(source) == receipt['artifacts'][source.name]
    write_json(args.output/'receipt.json', dict(
        status='complete_residual_summary_figure_pending_visual_inspection',
        source_completion_sha256=sha(args.completion), source_receipt_sha256=sha(receipt_path),
        source_table_sha256=sha(source), script_sha256=sha(__file__), plotted_rows=15,
        full_source_rows=35, all_extrema_within_axes=True,
        artifacts={p.name: sha(p) for p in args.output.iterdir()},
        scope='Descriptive ranges across overlapping fits; no adequacy test or biological inference.'))


if __name__ == '__main__':
    main()
