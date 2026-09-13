#!/usr/bin/env python3
"""Plot descriptive tree-path versus structural-geometry benchmarks."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--benchmark', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt = checked_receipt(args.benchmark)
    if args.output.exists():
        raise FileExistsError('Use a new immutable figure output')
    rows = list(csv.DictReader((args.benchmark / 'path_geometry.tsv').open(), delimiter='\t'))
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 9.0), layout='constrained')
    panels = [
        ('aa_path_distance', '3di_af_path_distance', 'AA tree path (substitutions/site)', '3Di tree path (state substitutions/site)', 'Sequence and structural model distances'),
        ('uncorrected_3di_difference', '3di_af_path_distance', 'Observed 3Di mismatch fraction', '3Di tree path (state substitutions/site)', 'Observed mismatch and fitted path distance'),
        ('3di_af_path_distance', 'ca_superposition_rmsd_angstrom', '3Di tree path (state substitutions/site)', 'Whole-marker CA RMSD (Å)', 'Global geometry includes domain placement'),
        ('3di_af_path_distance', 'pae10_local_mean_absolute_change_angstrom', '3Di tree path (state substitutions/site)', 'PAE-filtered local distance change (Å)', 'Local geometry with directional PAE ≤10 Å')]
    counts = {}
    for ax, (xkey, ykey, xlabel, ylabel, title) in zip(axes.ravel(), panels):
        accepted = [r for r in rows if r[xkey] != '' and r[ykey] != '']
        x = np.array([float(r[xkey]) for r in accepted]); y = np.array([float(r[ykey]) for r in accepted])
        colors = [float(r['fraction_of_tree_columns']) for r in accepted]
        scatter = ax.scatter(x, y, c=colors, cmap='viridis', vmin=0, vmax=1, s=15, alpha=.65, edgecolors='none')
        if 'rmsd' in ykey:
            ax.set_yscale('symlog', linthresh=1)
        ax.set_xlim(left=0); ax.set_ylim(bottom=0)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(xlabel, fontsize=9); ax.set_ylabel(ylabel, fontsize=9)
        ax.text(.03, .95, f'{len(accepted)} taxon pairs', transform=ax.transAxes, va='top', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
        counts[ykey + '_vs_' + xkey] = len(accepted)
    fig.colorbar(scatter, ax=axes, shrink=.6, label='Fraction of tree alignment observed in both taxa')
    fig.suptitle('Fitted evolutionary paths and direct structural comparisons\nDependent observations; descriptive benchmark, not an acceleration test', fontsize=12)
    args.output.mkdir(parents=True)
    fig.savefig(args.output / 'tree_path_geometry.svg', metadata={'Date': None})
    fig.savefig(args.output / 'tree_path_geometry.png', dpi=180)
    result = {'status': 'complete_tree_path_geometry_figure',
        'benchmark_receipt_sha256': sha(args.benchmark / 'receipt.json'), 'script_sha256': sha(Path(__file__)),
        'points_by_panel': counts, 'source_marker_count': receipt['markers'],
        'interpretation': 'No regression, independent-observation correlation or significance test. Whole-marker RMSD includes uncertain relative domain placement. Geometry uses shared observed positions; tree fits use all retained marker positions and taxa. RMSD axis is linear below 1 Angstrom and logarithmic above it.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
