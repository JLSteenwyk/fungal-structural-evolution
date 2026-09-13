#!/usr/bin/env python3
"""Plot marker-level alignment coverage and residue-pair correspondence."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from build_species_matrix import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison', type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads((args.comparison / 'receipt.json').read_text())
    path = args.comparison / 'marker_correspondence.tsv'
    if sha(path) != receipt['artifacts'][path.name]:
        raise ValueError('Changed correspondence results')
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    profile = np.array([int(r['profile_columns']) for r in rows])
    mafft = np.array([int(r['mafft_columns']) for r in rows])
    informative = [r for r in rows if r['edge_jaccard'] != '']
    score = np.array([float(r['edge_jaccard']) for r in informative])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    axes[0].scatter(profile, mafft, color='#306998', s=22, alpha=.65, linewidths=0)
    maximum = max(profile.max(), mafft.max()) * 1.03
    axes[0].plot([0, maximum], [0, maximum], '--', color='.6', linewidth=.8)
    axes[0].set(xlabel='Profile columns retained at 50% occupancy', ylabel='MAFFT columns retained at 50% occupancy',
                title='Alignment coverage per marker')
    axes[1].hist(score, bins=np.linspace(0, 1, 21), color='#39877a', edgecolor='white')
    axes[1].set(xlabel='Residue-pair edge Jaccard similarity', ylabel='Number of markers', xlim=(0, 1),
                title='Correspondence on common retained residues')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle(f'{len(rows)} markers · HMM match states versus full-protein MAFFT alignments')
    fig.text(.5, .02, 'Agreement is not accuracy. Shared residue-pair assignments are conditional on coverage retained by both methods.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .06, 1, .95))
    artifacts = {}
    for extension in ['svg', 'png', 'pdf']:
        path = args.comparison / f'alignment_correspondence.{extension}'
        fig.savefig(path, dpi=200)
        artifacts[path.name] = sha(path)
    plt.close(fig)
    summary = {'source_receipt_sha256': sha(args.comparison / 'receipt.json'),
        'script_sha256': sha(Path(__file__)), 'markers': len(rows),
        'median_marker_edge_jaccard': float(np.median(score)),
        'min_marker_edge_jaccard': float(score.min()), 'max_marker_edge_jaccard': float(score.max()),
        'profile_columns': int(profile.sum()), 'mafft_columns': int(mafft.sum()), 'artifacts': artifacts}
    (args.comparison / 'figure_receipt.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
