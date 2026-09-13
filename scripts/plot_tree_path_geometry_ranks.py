#!/usr/bin/env python3
"""Plot all within-marker rank associations with paired model comparisons."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from audit_busco_gene_copies import sha
from assess_pae_sensitivity import checked_receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--source-label', required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use an immutable figure output')
    receipt = checked_receipt(a.summary)
    with (a.summary / 'marker_rank_associations.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    models = ['aa_tree_path_point', '3di_af_tree_path_point', '3di_af_empirical_tree_path_point', '3di_llm_tree_path_point']
    labels = ['AA\nLG+F+G4', '3Di\nAF+G4', '3Di\nAF+F+G4', '3Di\nLLM+G4']
    metrics = ['ca_superposition_rmsd_angstrom', 'local_distance_mean_absolute_change_angstrom', 'pae10_local_mean_absolute_change_angstrom']
    titles = ['Whole-protein Cα RMSD', 'Local distance change', 'PAE ≤10 Å local distance change']
    indexed = {(r['marker'], r['path_metric'], r['geometry_metric']): r for r in rows}
    markers = sorted({r['marker'] for r in rows})
    if len(rows) != len(indexed) or len(rows) != len(markers) * 12 or len(markers) != receipt['markers']:
        raise ValueError('Incomplete or duplicate marker/model/metric grid')
    # Same deterministic offset for each marker across models and panels.
    offsets = {m: (int(hashlib.sha256(m.encode()).hexdigest()[:8], 16) / (2**32 - 1) - .5) * .20 for m in markers}
    fig, axes = plt.subplots(1, 3, figsize=(12, 5.5), sharey=True)
    displayed = 0; omitted = 0; medians = {}
    for ax, metric, title in zip(axes, metrics, titles):
        vectors = []
        for marker in markers:
            selected = [indexed[marker, model, metric] for model in models]
            if len({(r['pairs'], r['omitted_geometry_pairs']) for r in selected}) != 1:
                raise ValueError('Model cohorts differ')
            values = np.array([float(r['rank_correlation']) if r['rank_correlation'] else np.nan for r in selected])
            if np.any(np.abs(values[np.isfinite(values)]) > 1):
                raise ValueError('Invalid correlation')
            vectors.append(values)
            x = np.arange(4) + offsets[marker]
            ax.plot(x, values, color='#586878', lw=.55, alpha=.16, zorder=1)
            finite = np.isfinite(values)
            ax.scatter(x[finite], values[finite], s=10, color='#42677d', alpha=.50, linewidths=0, zorder=2)
            displayed += int(finite.sum()); omitted += int((~finite).sum())
        values = np.array(vectors)
        if not np.isfinite(values).any(axis=0).all():
            raise ValueError('Entire model unestimable')
        median = np.nanmedian(values, axis=0); medians[metric] = dict(zip(models, median.tolist()))
        ax.scatter(np.arange(4), median, marker='_', s=350, linewidths=2.5, color='#b84523', zorder=3)
        ax.axhline(0, color='.6', lw=.7, ls='--')
        ax.set_xticks(np.arange(4), labels, fontsize=9)
        ax.set(title=title, ylim=(-1.04, 1.04), xlim=(-.4, 3.4))
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=.12)
    axes[0].set_ylabel('Within-marker Spearman correlation')
    fig.suptitle(f'{a.source_label}: tree distances and direct geometry across {len(markers)} markers', fontsize=13)
    fig.text(.5, .11, 'Blue points: individual markers; lines connect the same marker across models. Orange bars: equal-marker medians.', ha='center', fontsize=9)
    fig.text(.5, .06, 'Descriptive point estimates. Shared ancestry, branch uncertainty and prediction effects remain unadjusted.', ha='center', fontsize=9)
    fig.text(.5, .025, 'AF and LLM identify 3Di substitution models; all structures in this figure are from the stated prediction source.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .15, 1, .94))
    a.output.mkdir(parents=True)
    for ext in ['svg', 'pdf', 'png']:
        fig.savefig(a.output / ('tree_path_geometry_ranks.' + ext), dpi=180)
    plt.close(fig)
    out = {'status': 'complete_descriptive_figure', 'source_receipt_sha256': sha(a.summary / 'receipt.json'),
           'script_sha256': sha(Path(__file__)), 'source_label': a.source_label, 'markers': len(markers),
           'displayed_points': displayed, 'unestimable_omitted_points': omitted, 'median_correlations': medians,
           'interpretation': receipt['interpretation'], 'artifacts': {f.name: sha(f) for f in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
