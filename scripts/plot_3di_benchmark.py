#!/usr/bin/env python3
"""Plot the matched-site 3Di/geometry benchmark and report confidence-regime coverage."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from compare_marker_structures import sha
from assess_pae_sensitivity import checked_receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--benchmark', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt = checked_receipt(args.benchmark)
    rows = list(csv.DictReader((args.benchmark / 'benchmark.tsv').open(), delimiter='\t'))
    exclusions = list(csv.DictReader((args.benchmark / 'exclusions.tsv').open(), delimiter='\t'))
    if args.output.exists():
        raise FileExistsError('Use a new immutable figure directory')
    args.output.mkdir(parents=True)
    summaries = []
    for scope in ['whole_marker', 'domain']:
        for regime in receipt['regimes']:
            selected = [r for r in rows if r['scope'] == scope and r['confidence_regime'] == regime]
            entry = {'scope': scope, 'confidence_regime': regime, 'comparison_rows': len(selected),
                     'excluded_rows': sum(r['scope'] == scope and r['confidence_regime'] == regime for r in exclusions)}
            for field in ['retained_fraction', 'compared_residues', 'uncorrected_3di_difference', 'uncorrected_sequence_difference',
                          'ca_superposition_rmsd_angstrom', 'local_distance_mean_absolute_change_angstrom']:
                values = [float(r[field]) for r in selected if r[field] != '']
                entry['median_' + field] = float(np.median(values)) if values else ''
            summaries.append(entry)
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.7), layout='constrained')
    selected_regime = 'six_residue70_pae10'
    for i, scope in enumerate(['whole_marker', 'domain']):
        selected = [r for r in rows if r['scope'] == scope and r['confidence_regime'] == selected_regime]
        seq = np.array([float(r['uncorrected_sequence_difference']) for r in selected])
        three = np.array([float(r['uncorrected_3di_difference']) for r in selected])
        rmsd = np.array([float(r['ca_superposition_rmsd_angstrom']) for r in selected])
        retained = np.array([float(r['retained_fraction']) for r in selected])
        local = np.array([float(r['local_distance_mean_absolute_change_angstrom']) if r['local_distance_mean_absolute_change_angstrom'] else np.nan for r in selected])
        for ax, x, y in zip(axes[i], [seq, three, three], [three, rmsd, local]):
            points = ax.scatter(x, y, c=retained, cmap='viridis', vmin=.5, vmax=1., s=10, alpha=.65, linewidths=0)
            ax.spines[['top', 'right']].set_visible(False)
        axes[i, 0].set(xlabel='Amino-acid mismatch fraction', ylabel='3Di mismatch fraction', xlim=(-.02, 1.02), ylim=(-.02, 1.02))
        axes[i, 1].set(xlabel='3Di mismatch fraction', ylabel='Matched-site Cα RMSD (Å)', xlim=(-.02, 1.02))
        axes[i, 1].set_yscale('symlog', linthresh=1)
        axes[i, 1].set_ylim(bottom=0)
        axes[i, 2].set(xlabel='3Di mismatch fraction', ylabel='Local mean |Δ distance| (Å)', xlim=(-.02, 1.02))
        axes[i, 2].set_yscale('symlog', linthresh=1)
        axes[i, 2].set_ylim(bottom=0)
        for ax in axes[i]:
            ax.set_title(('Whole marker' if scope == 'whole_marker' else 'Conserved domain') + f' · {len(selected)} comparisons', fontsize=10)
    fig.colorbar(points, ax=axes.ravel().tolist(), shrink=.7, pad=.025).set_label('Fraction of original geometric sites retained')
    fig.suptitle('Coordinate-derived 3Di states and direct geometry at identical aligned residues\nSix feature residues pLDDT ≥70; all directional context PAE ≤10 Å\nDescriptive comparisons; ancestry and protein-family effects unadjusted', fontsize=11)
    for suffix in ['svg', 'png', 'pdf']:
        fig.savefig(args.output / ('3di_geometry.' + suffix), dpi=180)
    plt.close(fig)
    with (args.output / 'confidence_summary.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, list(summaries[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(summaries)
    result = {'benchmark_receipt_sha256': sha(args.benchmark / 'receipt.json'), 'script_sha256': sha(Path(__file__)),
        'plotted_regime': selected_regime, 'scope_regime_summary': summaries,
        'interpretation': 'Shared-site descriptive comparisons with dependent taxon pairs and family effects. No substitution correction, significance test, branch-rate validation or conversion from alphabet states to physical displacement.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
