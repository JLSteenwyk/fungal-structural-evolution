#!/usr/bin/env python3
"""Plot descriptive paired measurements without treating pairs as independent tests."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparisons', type=Path, required=True)
    args = parser.parse_args()
    path = args.comparisons / 'pairwise_metrics.tsv'
    receipt = json.loads((args.comparisons / 'receipt.json').read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest() != receipt['artifacts']['pairwise_metrics.tsv']:
        raise ValueError('Changed comparison table')
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    selected = [r for r in rows if r['plddt_cutoff'] == '70']
    if not selected:
        raise ValueError('No qualifying pLDDT70 comparisons')
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, metric, label in zip(axes, ['ca_superposition_rmsd_angstrom', 'local_distance_mean_absolute_change_angstrom'],
                                ['Cα superposition RMSD (Å)', 'Mean local Cα distance change (Å)']):
        usable = [r for r in selected if r[metric] != '']
        ax.scatter([float(r['uncorrected_sequence_difference']) for r in usable],
                   [float(r[metric]) for r in usable], s=16, alpha=.35, color='#286c8e', linewidths=0)
        ax.set_xlabel('Uncorrected sequence difference at compared sites')
        ax.set_ylabel(label)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.15)
    fig.suptitle(f'Direct marker comparisons: {len(selected):,} taxon pairs across {len({r["marker"] for r in selected})} markers', fontsize=12)
    fig.text(.5, .02, 'Both models require pLDDT ≥70 at compared residues; at least 50 residues and half of shared positions.\nAcquisition snapshot; pairs share taxa and ancestry. No correlation test or branch-rate inference is shown.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .1, 1, .95))
    for suffix in ['png', 'svg', 'pdf']:
        fig.savefig(args.comparisons / f'direct_comparisons.{suffix}', dpi=180)
    plt.close(fig)
    print('Plotted', len(selected), 'descriptive comparison rows')


if __name__ == '__main__':
    main()
