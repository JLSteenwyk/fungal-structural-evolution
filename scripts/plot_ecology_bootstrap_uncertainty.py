#!/usr/bin/env python3
"""Plot audited ecological mapping uncertainty with explicit bootstrap denominators."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, required=True)
    args = parser.parse_args()
    root = Path('results/ecology/bootstrap-uncertainty-summary-20260927-v1')
    receipt = json.loads((root / 'receipt.json').read_text())
    proof_path = Path('metadata/ecology_bootstrap_uncertainty_summary_readback_20260927.json')
    proof = json.loads(proof_path.read_text())
    assert proof['status'] == 'passed_full_ecology_bootstrap_summary_readback'
    assert proof['source_receipt_sha256'] == sha(root / 'receipt.json')
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
    for suffix in ['.png', '.pdf', '.svg', '.tsv', '.receipt.json']:
        if args.prefix.with_suffix(suffix).exists():
            raise FileExistsError('Use a new figure prefix')
    split = pd.read_csv(root / 'split_uncertainty.tsv', sep='\t')
    distribution = pd.read_csv(root / 'tree_metric_distributions.tsv', sep='\t')
    conditions = [('profile_profile', 'original32_coding'), ('profile_mafft', 'original32_coding'),
                  ('profile_profile', 'ramaria_unknown'), ('profile_mafft', 'ramaria_unknown')]
    labels = ['Profile/profile\nOriginal', 'Profile/MAFFT\nOriginal',
              'Profile/profile\nRamaria unknown', 'Profile/MAFFT\nRamaria unknown']
    tips = ['F264124', 'F68786', 'F113071']
    tip_labels = ['Botryobasidium botryosum', 'Sphaerobolus stellatus', 'Ramaria rubella']
    assert set(split.loc[split.required_change > 0, 'canonical_side']) == set(tips)
    rows = []
    values = np.zeros((3, 4))
    for i, tip in enumerate(tips):
        for j, (source, coding) in enumerate(conditions):
            selected = split[(split.tree_source == source) & (split.coding == coding) & (split.canonical_side == tip)]
            assert len(selected) == 1
            r = selected.iloc[0]
            assert r.present == r.ensemble_trees == 1000
            values[i, j] = r.required_change / 1000
            rows.append(dict(panel='required', tree_source=source, coding=coding, taxon_id=tip,
                             edge_count='', trees=int(r.required_change), denominator=1000, fraction=values[i, j]))
    histograms = {}
    for source in ['profile_profile', 'profile_mafft']:
        subsets = []
        for coding in ['original32_coding', 'ramaria_unknown']:
            d = distribution[(distribution.tree_source == source) & (distribution.coding == coding) & (distribution.metric == 'optional_change')]
            counts = d.set_index('value').trees.reindex(range(17, 25), fill_value=0)
            assert counts.sum() == 1000
            subsets.append(counts)
            for edge_count, count in counts.items():
                rows.append(dict(panel='optional_distribution', tree_source=source, coding=coding,
                                 taxon_id='', edge_count=edge_count, trees=count, denominator=1000, fraction=count / 1000))
        assert subsets[0].equals(subsets[1]), 'Cannot combine differing coding distributions'
        histograms[source] = subsets[0].to_numpy()
    scores = distribution[distribution.metric == 'minimum_changes']
    assert len(scores) == 4 and (scores.trees == 1000).all()
    for r in scores.itertuples(index=False):
        assert r.value == (7 if r.coding == 'original32_coding' else 6)
        rows.append(dict(panel='minimum_score', tree_source=r.tree_source, coding=r.coding,
                         taxon_id='', edge_count=r.value, trees=r.trees, denominator=1000, fraction=1.0))
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={'width_ratios': [1.3, 1]})
    fig.subplots_adjust(left=.20, right=.98, top=.78, bottom=.30, wspace=.38)
    fig.suptitle('Ecological state changes: stable counts, uncertain locations', fontsize=16, x=.52, y=.96)
    fig.text(.20, .87, '2 guide ensembles × 1,000 bootstrap trees; root-free, equal-cost state mapping', fontsize=11)
    left.imshow(values, vmin=0, vmax=1, cmap='Blues', aspect='auto')
    left.set_yticks(range(3), tip_labels, fontstyle='italic')
    left.set_xticks(range(4), labels, fontsize=9)
    left.tick_params(length=0)
    left.set_title('A  Change required in every optimal mapping\nPercentage of bootstrap trees', loc='left', fontsize=11, pad=14)
    for i in range(3):
        for j in range(4):
            value = values[i, j]
            label = f'{100 * value:g}%' + ('*' if i == 2 and j >= 2 else '')
            left.text(j, i, label, ha='center', va='center', color='white' if value > .6 else '#222222', fontsize=11)
    x = np.arange(17, 25)
    for offset, source, label, color in [(-.19, 'profile_profile', 'Profile/profile', '#2166ac'),
                                         (.19, 'profile_mafft', 'Profile/MAFFT', '#d95f02')]:
        right.bar(x + offset, histograms[source] / 10, width=.38, label=label, color=color)
    right.set_xticks(x)
    right.set_xlabel('Edges where change is optional')
    right.set_ylabel('Bootstrap trees (%)')
    right.set_title('B  Ambiguity across complete trees\nSame distributions under both codings', loc='left', fontsize=11, pad=14)
    right.legend(frameon=False, fontsize=9)
    right.spines[['top', 'right']].set_visible(False)
    right.set_ylim(0, 65)
    fig.text(.03, .17, 'Every tree has minimum score 7 with original coding, or 6 with Ramaria unknown.', fontsize=11)
    fig.text(.03, .115, '* Ramaria is unassigned in this sensitivity coding; 0% does not establish ecological stasis.', fontsize=10)
    fig.text(.03, .065, 'Fractions describe mapping sensitivity, not posterior probabilities or independent origins. Optional edges are jointly constrained.', fontsize=9)
    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in ['png', 'pdf', 'svg']:
        fig.savefig(args.prefix.with_suffix('.' + extension), dpi=180, facecolor='white')
    plt.close(fig)
    pd.DataFrame(rows).to_csv(args.prefix.with_suffix('.tsv'), sep='\t', index=False)
    result = dict(status='complete_ecology_bootstrap_uncertainty_figure_pending_review',
                  source_receipt_sha256=sha(root / 'receipt.json'), source_readback_sha256=sha(proof_path),
                  script_sha256=sha(__file__), data_rows=len(rows),
                  artifacts={args.prefix.with_suffix('.' + ext).name: sha(args.prefix.with_suffix('.' + ext)) for ext in ['png', 'pdf', 'svg', 'tsv']},
                  scope='All ever-required terminal edges shown across four conditions; all optional-count probabilities including zero bins retained. Original and Ramaria-unknown histograms plotted once only after exact equality checks. Minimum scores annotate audited distributions. Not independent origins or ecological effects.')
    args.prefix.with_suffix('.receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
