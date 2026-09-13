#!/usr/bin/env python3
"""Plot descriptive PAE sensitivity and rank pairs for manual structural review."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from compare_marker_structures import ROOT, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparisons', type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads((args.comparisons / 'receipt.json').read_text())
    path = args.comparisons / 'pae_sensitivity.tsv'
    if sha(path) != receipt['artifacts'][path.name]:
        raise ValueError('Changed sensitivity results')
    with path.open() as handle:
        rows = [r for r in csv.DictReader(handle, delimiter='\t')
                if r['pae_cutoff_angstrom'] == '10' and r['residue_pair_scope'] == 'all_nonadjacent']
    qualified = [r for r in rows if r['confident_mean_absolute_distance_change_angstrom'] != '']
    x = np.array([float(r['unfiltered_mean_absolute_distance_change_angstrom']) for r in qualified])
    y = np.array([float(r['confident_mean_absolute_distance_change_angstrom']) for r in qualified])
    confidence = np.array([float(r['confident_fraction']) for r in qualified])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))
    points = axes[0].scatter(x, y, c=confidence, cmap='viridis', vmin=0, vmax=1, s=13, alpha=.8, linewidths=0)
    maximum = max(float(x.max()), float(y.max())) * 1.03
    axes[0].plot([0, maximum], [0, maximum], '--', color='.55', linewidth=.8)
    axes[0].set(xlabel='All eligible residue pairs: mean |Δ distance| (Å)',
                ylabel='PAE-confident pairs: mean |Δ distance| (Å)', title='Matched-site structural differences')
    axes[1].scatter([float(r['ca_superposition_rmsd_angstrom']) for r in qualified], confidence,
                    s=13, alpha=.55, color='#306998', linewidths=0)
    axes[1].set(xlabel='Global Cα superposition RMSD (Å)', ylabel='Fraction of residue pairs passing PAE filter',
                ylim=(-.025, 1.025), title='Confidence and global superposition')
    colorbar = fig.colorbar(points, ax=axes[0], fraction=.046, pad=.04)
    colorbar.set_label('Fraction passing PAE filter')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle(f'PAE ≤10 Å in both directions in both models · {len(qualified)} taxon pairs', fontsize=12)
    fig.text(.5, .02, 'pLDDT ≥70; nonadjacent residues. Descriptive, shared ancestry unadjusted; filters change pair composition.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .06, 1, .95))
    artifacts = {}
    for suffix in ['svg', 'png', 'pdf']:
        output = args.comparisons / f'pae_sensitivity.{suffix}'
        fig.savefig(output, dpi=200)
        artifacts[output.name] = sha(output)
    plt.close(fig)
    taxa_path = ROOT / 'metadata/analysis_manifest.tsv'
    with taxa_path.open() as handle:
        taxa = {r['taxon_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    ranked = []
    for rank, r in enumerate(sorted(rows, key=lambda r: (-float(r['ca_superposition_rmsd_angstrom']),
                              r['marker'], r['taxon_a'], r['taxon_b'])), 1):
        ranked.append(dict(review_rank=rank, species_a=taxa[r['taxon_a']]['species_name'],
                           species_b=taxa[r['taxon_b']]['species_name'], **r))
    output = args.comparisons / 'global_rmsd_review_order.tsv'
    with output.open('w') as handle:
        writer = csv.DictWriter(handle, list(ranked[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(ranked)
    artifacts[output.name] = sha(output)
    result = {'source_receipt_sha256': sha(args.comparisons / 'receipt.json'),
        'taxon_manifest_sha256': sha(taxa_path), 'script_sha256': sha(Path(__file__)),
        'taxon_pairs_plotted': len(qualified), 'pairs_without_confident_residue_pairs': len(rows) - len(qualified),
        'median_confident_residue_pair_fraction': float(np.median(confidence)),
        'median_unfiltered_mean_absolute_distance_change_angstrom': float(np.median(x)),
        'median_confident_mean_absolute_distance_change_angstrom': float(np.median(y)),
        'review_order_interpretation': 'All assessed pairs ordered by global RMSD for manual review, not statistical discoveries or independent evolutionary events.',
        'artifacts': artifacts}
    (args.comparisons / 'figure_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
