#!/usr/bin/env python3
"""Plot descriptive assignment sensitivity and its model/hit denominators."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter, FuncFormatter
from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--completion', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    completion = json.loads(args.completion.read_text())
    receipt = Path(completion['receipt'])
    assert sha(receipt) == completion['receipt_sha256']
    source = receipt.parent / 'cluster_change_by_extension.tsv'
    assert sha(source) == completion['artifacts'][source.name]
    with source.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert len(rows) == len(completion['summaries']) == 6
    for row, expected in zip(rows, completion['summaries']):
        assert all(type(value)(row[key]) == value for key, value in expected.items())
    labels = [r['added_residue_bin'] for r in rows]
    n = [int(r['hits']) for r in rows]
    changed = [int(r['changed_cluster']) for r in rows]
    fractions = [k / total for k, total in zip(changed, n)]
    assert sum(n) == 868338 and sum(changed) == 130084
    assert all(f == float(r['changed_fraction']) for f, r in zip(fractions, rows))
    args.output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.7))
    fig.subplots_adjust(left=.09, right=.98, bottom=.23, top=.79, wspace=.35)
    fig.suptitle('Domain boundary choice and structural cluster assignment', y=.97, fontsize=15)
    fig.text(.5, .89, '868,338 model/hit pairs · expanded September 28 catalog', ha='center', color='#555555')
    bars = axes[0].bar(labels, fractions, color='#287C8E', width=.7)
    axes[0].set(ylabel='Different cluster assignments', ylim=(0, 1),
                xlabel='Added residues: envelope − alignment', title='A   Fraction of pairs')
    axes[0].yaxis.set_major_formatter(PercentFormatter(1))
    axes[0].bar_label(bars, labels=[f'{100*f:.1f}%' for f in fractions], padding=4, fontsize=9)
    bars = axes[1].bar(labels, n, color='#8697A7', width=.7)
    axes[1].set(ylabel='Model/hit pairs', ylim=(0, max(n)*1.2),
                xlabel='Added residues: envelope − alignment', title='B   Denominators')
    axes[1].yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value/1000:.0f}k'))
    axes[1].bar_label(bars, labels=[f'{value:,}' for value in n], padding=4, fontsize=8)
    for ax in axes:
        ax.set_axisbelow(True)
        ax.grid(axis='y', color='#E6E6E6', linewidth=.6)
    fig.text(.09, .07, 'Descriptive comparison within one partition; hits are not independent evolutionary replicates.\n'
             'Boundary bins are reporting categories. Assignment changes do not establish biological divergence.',
             fontsize=9, color='#555555', va='center')
    artifacts = {}
    for extension in ['png', 'pdf']:
        path = args.output / ('boundary_cluster_extension.'+extension)
        fig.savefig(path, dpi=200)
        artifacts[path.name] = sha(path)
    plt.close(fig)
    result = dict(status='rendered_boundary_cluster_extension_figure',
                  source_hashes={str(p): sha(p) for p in [args.completion, receipt, source]},
                  script_sha256=sha(__file__), artifacts=artifacts,
                  pairs=sum(n), changed_assignments=sum(changed), plotted_rows=rows,
                  scope='All six bins checked against completed table; ratios recomputed from integers. Rendering and arithmetic only, not an independent source join or inferential test.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
