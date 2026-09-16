#!/usr/bin/env python3
"""Plot every expanded ESMFold specification with its copy-review sensitivity."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['comparison', 'full', 'omission', 'output']:
        ap.add_argument('--' + name, type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    comparison = json.loads((a.comparison / 'receipt.json').read_text())
    if comparison['status'] != 'passed_copy_review_coupling_comparison':
        raise ValueError('Audited full/omission comparison required')
    if (comparison['full_markers'], comparison['omission_markers']) != (89, 88):
        raise ValueError('Unexpected figure cohorts')
    source = {}
    pins = {str(a.comparison / 'receipt.json'): sha(a.comparison / 'receipt.json')}
    for label, folder in [('full', a.full), ('omission', a.omission)]:
        receipt = json.loads((folder / 'receipt.json').read_text())
        if comparison['source_receipts'][str(folder)] != sha(folder / 'receipt.json'):
            raise ValueError('Comparison source mismatch')
        path = folder / 'contrast_summary.tsv'
        if sha(path) != receipt['artifacts'][path.name]:
            raise ValueError('Changed contrast summary')
        pins[str(path)] = sha(path)
        source[label] = pd.read_csv(path, sep='\t').set_index(['model_id', 'contrast'], verify_integrity=True)
    grid = list(itertools.product(['3di_af', '3di_af_empirical', '3di_llm'],
                                  ['gamma', 'freerate'], ['tien2013', 'miller1987'],
                                  ['coverage_confidence', 'composition_adjusted']))
    contrasts = [('aa_log1p_rate_at_RSA_0.25', 'Sequence-rate association\nat RSA = 0.25'),
                 ('rsa_at_aa_log1p_rate_0', 'RSA association\nat log(1 + AA rate) = 0'),
                 ('interaction', 'Sequence-rate × RSA\ninteraction')]
    colors = {'full': '#146C94', 'omission': '#C65C24'}
    labels = {'3di_af': '3Di-AF', '3di_af_empirical': '3Di-AF+F', '3di_llm': '3Di-LLM',
              'gamma': 'G4', 'freerate': 'R4', 'tien2013': 'Tien', 'miller1987': 'Miller',
              'coverage_confidence': 'Base', 'composition_adjusted': '+Comp'}
    fig, axes = plt.subplots(1, 3, figsize=(14, 10.5), sharey=True)
    plotted = []
    for ax, (contrast, title) in zip(axes, contrasts):
        ax.axvline(0, color='#777777', linestyle='--', linewidth=.8)
        for start in [0, 16]:
            ax.axhspan(start - .5, start + 7.5, color='#F1F4F6', zorder=0)
        for y, spec in enumerate(grid):
            model = '__'.join(spec)
            for label, offset in [('full', -.14), ('omission', .14)]:
                row = source[label].loc[(model, contrast)]
                point, lower, upper = (float(row[x]) for x in ['point_estimate', 'bootstrap_percentile_lower', 'bootstrap_percentile_upper'])
                if not np.isfinite([point, lower, upper]).all() or lower > upper:
                    raise ValueError('Invalid plotted interval')
                ax.plot([lower, upper], [y + offset, y + offset], color=colors[label], lw=1.1, alpha=.9)
                ax.plot(point, y + offset, 'o' if label == 'full' else 's', color=colors[label], markersize=3.3)
                plotted.append(dict(analysis=label, model_id=model, contrast=contrast,
                                    point_estimate=point, lower=lower, upper=upper, row=y))
        ax.set_title(title, fontsize=11, pad=12)
        ax.set_xlabel('Conditional coefficient', fontsize=10)
        ax.grid(axis='y', alpha=.18)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(axis='x', labelsize=9)
    axes[0].set_yticks(range(24), [' | '.join(labels[x] for x in spec) for spec in grid], fontsize=8)
    axes[0].set_ylim(23.7, -.7)
    fig.suptitle('Expanded ESMFold coupling: full cohort and copy-review omission', fontsize=15, y=.975)
    handles = [Line2D([0], [0], color=colors['full'], marker='o', markersize=5, label='Full: 89 markers, 22,205 sites'),
               Line2D([0], [0], color=colors['omission'], marker='s', markersize=5, label='Omission: 88 markers, 21,904 sites')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.60, .936), ncol=2, frameon=False, fontsize=10)
    fig.text(.025, .077, 'Lines: unadjusted 95% intervals from 2,000 whole-marker bootstrap draws per specification. All 24 specifications shown.', fontsize=9)
    fig.text(.025, .056, 'Base: coverage and confidence controls; +Comp: additional amino-acid composition controls. All fits include marker intercepts.', fontsize=9)
    fig.text(.025, .035, 'Omission excludes marker 4986044at2759. Conditional on estimated rates, fixed gene trees and sequence-derived predictions.', fontsize=9)
    fig.text(.025, .014, 'These are model-relative associations, not causal effects, physical displacement or confirmed single-ortholog inference.', fontsize=9)
    fig.subplots_adjust(left=.285, right=.985, bottom=.14, top=.85, wspace=.23)
    a.output.mkdir(parents=True)
    table = a.output / 'plotted_intervals.tsv'
    pd.DataFrame(plotted).to_csv(table, sep='\t', index=False)
    # Read exported coordinates back against independently indexed source rows.
    back = pd.read_csv(table, sep='\t')
    assert len(back) == 144 and not back.duplicated(['analysis', 'model_id', 'contrast']).any()
    for row in back.itertuples():
        raw = source[row.analysis].loc[(row.model_id, row.contrast)]
        assert np.allclose([row.point_estimate, row.lower, row.upper],
                           raw[['point_estimate', 'bootstrap_percentile_lower', 'bootstrap_percentile_upper']].to_numpy(float), rtol=0, atol=1e-14)
    for suffix in ['svg', 'png']:
        fig.savefig(a.output / ('expanded_copy_review_coupling.' + suffix), dpi=150)
    plt.close(fig)
    assert all(sha(Path(p)) == h for p, h in pins.items())
    receipt = dict(status='complete_expanded_copy_review_coupling_figure', plotted_rows=144,
                   plotted_numeric_values_checked=432, model_specifications=24,
                   source_pins=pins, script_sha256=sha(Path(__file__)),
                   artifacts={p.name: sha(p) for p in a.output.iterdir()},
                   interpretation='Every specification and both reviewed cohorts shown with unadjusted marker-bootstrap intervals at explicit fixed reference values. Visual inspection remains separate. Not a test of differences between cohorts.')
    (a.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
