#!/usr/bin/env python3
"""Plot audited, descriptive boundary-extension cluster disagreement."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    pins = {str(p): sha(p) for p in (args.table, args.receipt)}
    receipt = json.loads(args.receipt.read_text())
    if receipt['status'] != 'complete_boundary_extension_cluster_summary':
        raise ValueError('Summary incomplete')
    if sha(args.table) != receipt['artifacts']['extension_cluster_disagreement.tsv']:
        raise ValueError('Table hash mismatch')
    with args.table.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if [r['added_residues'] for r in rows] != ['0', '1-4', '5-9', '10-19', '20-49', '50+']:
        raise ValueError('Unexpected bins')
    pairs = [int(r['pairs']) for r in rows]
    different = [int(r['different_clusters']) for r in rows]
    if sum(pairs) != receipt['candidate_pairs']:
        raise ValueError('Pair count mismatch')
    for r, n, k in zip(rows, pairs, different):
        if n <= 0 or not 0 <= k <= n or abs(float(r['fraction_different']) - k/n) > 1e-14:
            raise ValueError('Invalid fraction')
    percentages = [100*k/n for k, n in zip(different, pairs)]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(10, 6.4))
    fig.subplots_adjust(left=.12, right=.97, bottom=.25, top=.80)
    fig.suptitle('Domain clustering is sensitive to boundary choice', y=.96, fontsize=17)
    fig.text(.12, .88, f'{sum(pairs):,} alignment/envelope pairs in one joint Foldseek partition', fontsize=11)
    bars = ax.bar(range(6), percentages, color='#346c91', width=.65)
    ax.set_xticks(range(6), ['0\n(identical)', '1–4', '5–9', '10–19', '20–49', '50+'])
    ax.set_xlabel('Residues added by the HMM envelope')
    ax.set_ylabel('Pairs assigned to different clusters (%)')
    ax.set_ylim(0, 100)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis='y', alpha=.25)
    for bar, pct, k, n in zip(bars, percentages, different, pairs):
        if abs(bar.get_height()-pct) > 1e-12:
            raise ValueError('Bar height mismatch')
        ax.text(bar.get_x()+bar.get_width()/2, pct+2, f'{pct:.1f}%\n{k:,}/{n:,}', ha='center', va='bottom', fontsize=9)
    fig.text(.12, .11, 'Labels: different-cluster pairs / all pairs in each bin. Identical intervals share one database entry.', fontsize=9)
    fig.text(.12, .065, 'Descriptive census: shared proteins, families and ancestry; no causal or evolutionary-event inference.', fontsize=9)
    args.output.mkdir(parents=True, exist_ok=False)
    for ext in ('svg', 'pdf', 'png'):
        fig.savefig(args.output / f'boundary_cluster_extension.{ext}', dpi=170, facecolor='white')
    plt.close(fig)
    svg = args.output / 'boundary_cluster_extension.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    for p, digest in pins.items():
        if sha(p) != digest:
            raise ValueError('Source changed')
    result = dict(status='complete_boundary_cluster_extension_figure_pending_visual_review',
                  candidate_pairs=sum(pairs), different_cluster_pairs=sum(different), bars_checked=6,
                  source_hashes=pins, script_sha256=sha(__file__),
                  artifacts={p.name: sha(p) for p in args.output.iterdir()},
                  scope='Descriptive proportions within a single joint partition; not independent observations, separate-run stability, causal effects or evolutionary events.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
