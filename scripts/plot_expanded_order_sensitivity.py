#!/usr/bin/env python3
"""Plot complete expanded-order diagnostics with fixed common-mask denominators."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    pins = {str(args.plan): sha(args.plan), **plan['pins']}

    def verify():
        for path, digest in pins.items():
            if sha(path) != digest:
                raise ValueError('Changed figure source: ' + path)

    verify()
    root = Path(plan['source'])
    receipt = json.loads((root / 'receipt.json').read_text())
    audit = json.loads(Path(plan['audit']).read_text())
    assert audit['status'] == 'passed_full_primary_order_sensitivity_readback'
    assert audit['source_receipt_sha256'] == sha(root / 'receipt.json')
    assert audit['pair_mask_rows'] == receipt['pair_mask_rows'] == 250403
    assert audit['mapping_counts'] == receipt['mapping_counts']
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
    data = pd.read_csv(root / 'pair_mask_sensitivity.tsv', sep='\t', float_precision='round_trip')
    assert len(data) == receipt['pair_mask_rows'] and not data.duplicated(['pair_key', 'mask']).any()
    common = data[data.in_both_mask_cohort]
    assert len(common) == 2 * receipt['common_mask_pairs']
    masks = ['full', 'plddt70']
    assert set(common[common['mask'].eq('full')].pair_key) == set(common[common['mask'].eq('plddt70')].pair_key)
    summary = []
    for mask, cohort in zip(masks, ['full_common', 'plddt70_common']):
        sub = common[common['mask'].eq(mask)]
        counts = sub.mapping_status.value_counts().to_dict()
        assert counts == receipt['mapping_counts'][cohort]
        summary.append(dict(mask=mask, pairs=len(sub), **counts))
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    frame = pd.DataFrame(summary)
    frame.to_csv(out / 'mapping_counts.tsv', sep='\t', index=False)
    pd.testing.assert_frame_equal(frame, pd.read_csv(out / 'mapping_counts.tsv', sep='\t'))
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.5))
    same = 100 * frame.different_same_count / frame.pairs
    different = 100 * frame.different_count / frame.pairs
    axes[0].bar([0, 1], same, color='#0072B2', label='Same count; changed mapping')
    axes[0].bar([0, 1], different, bottom=same, color='#E69F00', label='Changed count')
    for i, row in frame.iterrows():
        changed = row.different_same_count + row.different_count
        axes[0].text(i, 100 * changed / row.pairs + .009,
                     f'{changed:,}/{row.pairs:,}', ha='center', va='bottom', fontsize=8)
    axes[0].set_xticks([0, 1], ['Full protein', 'pLDDT ≥70'])
    axes[0].set_ylim(0, max(same + different) * 1.65)
    axes[0].set_ylabel('Pairs with changed residue mapping (%)')
    axes[0].set_title('A  Both-mask common cohort', loc='left')
    axes[0].legend(frameon=False, fontsize=7, loc='upper left')
    curve_counts = {}
    for ax, mask, letter in zip(axes[1:], masks, ['B', 'C']):
        sub = common[common['mask'].eq(mask) & common.mapping_status.ne('identical')]
        values = np.sort(sub.rmsd_recomputed_absolute_difference.to_numpy())
        assert len(values) and np.isfinite(values).all() and (values >= 0).all()
        curve_counts[mask] = len(values)
        ax.step(np.r_[values[0], values], np.r_[0, np.arange(1, len(values) + 1) / len(values)],
                where='post', color='#0072B2')
        ax.set_xscale('symlog', linthresh=1e-6)
        ax.set_xlim(left=0)
        ax.set_ylim(0, 1.03)
        ax.set_title(f'{letter}  Changed mappings: ' + ('full' if mask == 'full' else 'pLDDT ≥70'), loc='left')
        ax.set_xlabel('Absolute input-order RMSD difference (Å)')
        ax.set_ylabel('Cumulative fraction of changed mappings')
        ax.text(.04, .92, f'n = {len(values):,}', transform=ax.transAxes)
    fig.suptitle('Expanded whole-protein input-order sensitivity', fontsize=13)
    fig.text(.5, .015, 'Same 115,591 pairs across masks. B–C condition on changed mappings and retain zero differences. Descriptive; pairs are not independent events.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .065, 1, .95))
    for suffix in ['png', 'pdf', 'svg']:
        fig.savefig(out / ('expanded_order_sensitivity.' + suffix), dpi=180)
    plt.close(fig)
    verify()
    result = dict(status='complete_source_checked_expanded_order_figure_pending_visual_review',
                  plan_sha256=sha(args.plan), source_hashes=pins, common_pairs=receipt['common_mask_pairs'],
                  conditional_curve_counts=curve_counts,
                  artifacts={p.name: sha(p) for p in out.iterdir()},
                  scope='Every common-mask mapping count checked against the independent full audit; every changed-mapping RMSD difference plotted, including zeros. '
                        'No order selection, confidence interval, independence assumption or biological effect estimate.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
