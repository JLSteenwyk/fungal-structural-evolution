#!/usr/bin/env python3
"""Plot the validated whole-protein candidate partition; no evolutionary inference."""
import argparse
from collections import Counter
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
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['composition', 'readback', 'summary', 'output']:
        ap.add_argument('--' + name, type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError('Use fresh output')
    table = a.composition / 'cluster_composition.tsv'
    receipt = a.composition / 'receipt.json'
    s = json.loads(a.summary.read_text())
    r = json.loads(a.readback.read_text())
    pr = json.loads(receipt.read_text())
    if (r['status'] != 'passed_full_cluster_composition_source_readback'
            or sha(a.readback) != s['readback_receipt_sha256']
            or sha(receipt) != s['producer_receipt_sha256']
            or r['producer_receipt_sha256'] != sha(receipt)
            or sha(table) != s['table_sha256']
            or pr['artifacts'][table.name] != sha(table)):
        raise ValueError('Source validation chain differs')
    pins = {str(p): sha(p) for p in [table, receipt, a.readback, a.summary]}
    bins = [(1, 1, '1'), (2, 4, '2–4'), (5, 9, '5–9'),
            (10, 49, '10–49'), (50, 99, '50–99'), (100, None, '≥100')]
    counts = {key: Counter() for key in ['models', 'taxa']}
    totals = Counter()
    seen = set()
    with table.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            rep = row['representative']
            if rep in seen:
                raise ValueError('Duplicate representative')
            seen.add(rep)
            values = {k: int(row[k]) for k in ['models', 'proteins', 'taxa']}
            if min(values.values()) < 1 or values['proteins'] < max(values['models'], values['taxa']):
                raise ValueError('Invalid cluster counts')
            totals.update(values)
            for key in counts:
                for lo, hi, label in bins:
                    if lo <= values[key] and (hi is None or values[key] <= hi):
                        counts[key][label] += 1
                        break
    n = len(seen)
    if (n != r['clusters'] or n != s['counts']['clusters']
            or totals['models'] != r['models'] or totals['proteins'] != r['protein_links']
            or n - counts['models']['1'] != s['counts']['multiple_models']
            or n - counts['taxa']['1'] != s['counts']['multiple_taxa']
            or any(sum(c.values()) != n for c in counts.values())):
        raise ValueError('Source totals differ')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=.09, right=.97, bottom=.25, top=.76, wspace=.3)
    fig.suptitle('Whole-protein structural candidate clusters', x=.06, ha='left', y=.97, fontsize=17)
    fig.text(.06, .89, f'{n:,} clusters · {totals["models"]:,} models · {totals["proteins"]:,} protein links', fontsize=11)
    data = []
    for ax, key, title, color in zip(axes, counts, ['A  Models per cluster', 'B  Taxa per cluster'], ['#346c91', '#bd7035']):
        vals = [counts[key][label] for _, _, label in bins]
        bars = ax.bar(range(len(bins)), vals, color=color, width=.65)
        ax.set_yscale('log'); ax.set_ylim(.8, max(vals)*5)
        ax.set_xticks(range(len(bins)), [x[2] for x in bins])
        ax.set_ylabel('Number of clusters (log scale)')
        ax.set_title(title, loc='left', pad=12)
        ax.spines[['top', 'right']].set_visible(False)
        for bar, (_, _, label), value in zip(bars, bins, vals):
            if bar.get_height() != value:
                raise ValueError('Bar differs from source count')
            ax.text(bar.get_x()+bar.get_width()/2, max(value, 1)*1.15,
                    f'{value:,}', ha='center', va='bottom', fontsize=9)
            data.append(dict(dimension=key, bin=label, clusters=value, denominator=n, percent=100*value/n))
    fig.text(.06, .13, 'Bins are disjoint; both panels include all candidate clusters. A model may map to multiple proteins or taxa.', fontsize=10)
    fig.text(.06, .075, 'Descriptive atlas census: clusters are not orthogroups or evolutionary events. Confidence filtering remains pending.', fontsize=10, color='#555555')
    a.output.mkdir(parents=True)
    with (a.output/'plot_data.tsv').open('w') as f:
        w = csv.DictWriter(f, list(data[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(data)
    for ext in ['svg', 'pdf', 'png']:
        fig.savefig(a.output/('whole_proteome_cluster_census.'+ext), dpi=170, facecolor='white')
    plt.close(fig)
    p = a.output/'whole_proteome_cluster_census.svg'
    p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines())+'\n')
    if any(sha(p) != h for p, h in pins.items()):
        raise ValueError('Source changed during plotting')
    out = dict(status='complete_descriptive_cluster_census_pending_visual_review',
               clusters=n, models=totals['models'], protein_links=totals['proteins'],
               source_hashes=pins, script_sha256=sha(__file__), bars_checked=len(data),
               artifacts={p.name:sha(p) for p in a.output.iterdir()},
               scope='Full composition-table census and all bar heights checked; no confidence, homology, orthology or evolutionary inference.')
    (a.output/'receipt.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
