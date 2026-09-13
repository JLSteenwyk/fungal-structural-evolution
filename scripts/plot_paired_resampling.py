#!/usr/bin/env python3
"""Plot block-length sensitivity of conditional AA and 3Di branch intervals."""
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
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checked_receipt(args.audit)
    if args.output.exists():
        raise FileExistsError('Use a new immutable figure directory')
    rows = list(csv.DictReader((args.audit / 'conditional_intervals.tsv').open(), delimiter='\t'))
    lookup = {(r['marker'], r['split_taxa'], r['fit'], int(r['block_length'])): r for r in rows}
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8))
    summary = {}
    for ax, label, title, color in zip(axes, ['aa', '3di'], ['Amino-acid branches', 'Structural-state branches (AF model)'], ['#27618c', '#ad5224']):
        widths = {block: [] for block in [1, 10, 30]}
        for key, row in lookup.items():
            if key[2:] != (label, 1):
                continue
            paired = [lookup[key[:3] + (block,)] for block in [1, 10, 30]]
            if any(r['interval_status'] != 'conditional_percentile_sensitivity' for r in paired):
                continue
            for block, r in zip([1, 10, 30], paired):
                width = float(r['percentile_97_5']) - float(r['percentile_2_5'])
                if width < 0:
                    raise ValueError('Invalid interval bounds')
                widths[block].append(width)
        x, y = np.array(widths[1]), np.array(widths[30])
        maximum = max(x.max(), y.max()) * 1.05
        ax.scatter(x, y, s=15, color=color, alpha=.5, edgecolors='none')
        ax.plot([0, maximum], [0, maximum], color='#555555', lw=.9, ls='--')
        ax.set_xscale('symlog', linthresh=.01); ax.set_yscale('symlog', linthresh=.01)
        ax.set_xlim(0, maximum); ax.set_ylim(0, maximum)
        ax.set_xlabel('Interval width: single-column resampling')
        ax.set_ylabel('Interval width: 30-column blocks')
        ax.set_title(title, fontsize=11)
        ax.text(.04, .96, f'{len(x)} matched branches\nDashed: equal widths', transform=ax.transAxes, va='top', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
        summary[label] = {'matched_branches': len(x), 'median_interval_width': {str(b): float(np.median(v)) for b, v in widths.items()},
                          'branches_wider_at_block30_than_block1': int(np.sum(y > x))}
    fig.suptitle('Branch estimates depend on resampling block length', fontsize=13)
    fig.text(.5, .015, 'Widths in expected state substitutions/site; axes linear near zero, logarithmic above 0.01.\nFixed topology and prediction source. Nonlocal feature dependence and other uncertainty remain unresolved.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .09, 1, .95))
    args.output.mkdir(parents=True)
    fig.savefig(args.output / 'paired_resampling.svg', metadata={'Date': None})
    fig.savefig(args.output / 'paired_resampling.png', dpi=180)
    receipt = {'status': 'complete_conditional_resampling_figure', 'audit_receipt_sha256': sha(args.audit / 'receipt.json'),
        'script_sha256': sha(Path(__file__)), 'summary': summary,
        'interpretation': 'Matched branch interval-width comparison across resampling blocks, not significance testing. Points share ancestry and are not independent. 95% percentile ranges are conditional sensitivities, not validated overall confidence bounds.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
