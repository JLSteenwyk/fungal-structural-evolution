#!/usr/bin/env python3
"""Plot domain-conditioned geometry without treating taxon pairs as independent events."""
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
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source-label', default='')
    args = parser.parse_args()
    receipt = json.loads((args.comparisons / 'receipt.json').read_text())
    path = args.comparisons / 'domain_comparisons.tsv'
    if sha(path) != receipt['artifacts'][path.name]:
        raise ValueError('Changed comparisons')
    rows = list(csv.DictReader(path.open(), delimiter='\t'))
    if not rows:
        raise ValueError('No comparisons to plot')
    if args.output.exists():
        raise FileExistsError('Use an immutable new figure output')
    args.output.mkdir(parents=True)
    x = np.array([float(r['domain_sites_under_whole_marker_fit_rmsd_angstrom']) for r in rows])
    y = np.array([float(r['domain_own_fit_rmsd_angstrom']) for r in rows])
    seq = np.array([float(r['domain_uncorrected_sequence_difference']) for r in rows])
    conf = np.array([float(r['pae10_local_pair_fraction']) if r['pae10_local_pair_fraction'] else np.nan for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.9))
    finite = np.isfinite(conf)
    dense = len(rows) > 5000
    size, opacity = (3, .2) if dense else (16, .7)
    for ax, values in zip(axes, [x, seq]):
        ax.scatter(values[~finite], y[~finite], s=size, color='.65', alpha=opacity, rasterized=dense)
        points = ax.scatter(values[finite], y[finite], c=conf[finite], cmap='viridis_r',
                            vmin=0, vmax=1, s=size, alpha=opacity, linewidths=0, rasterized=dense)
        ax.set_yscale('symlog', linthresh=1)
        ax.set_ylabel('Domain fitted independently: Cα RMSD (Å)')
        ax.spines[['top', 'right']].set_visible(False)
    maximum = max(x.max(), y.max()) * 1.05
    axes[0].plot([0, maximum], [0, maximum], '--', color='.55', lw=.8)
    axes[0].set_xscale('symlog', linthresh=1)
    axes[0].set(xlabel='Same domain sites under whole-marker fit: RMSD (Å)', title='Global placement and internal geometry')
    axes[1].set(xlabel='Sequence difference at the same domain sites', title='Sequence and domain structure', xlim=(-.01, 1.01))
    fig.colorbar(points, ax=axes[1], fraction=.047, pad=.04).set_label('Local residue-pair fraction passing PAE ≤10 Å')
    fig.suptitle(f'{args.source_label + chr(32) if args.source_label else chr(32)}{len(rows):,} domain–taxon-pair comparisons · {receipt["distinct_pfam_domains"]} Pfam domains', fontsize=12)
    fig.text(.5, .035, 'Separate fitting necessarily improves fit. Dependent descriptive observations; ancestry and prediction-source effects unadjusted.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .075, 1, .96))
    for suffix in ['svg', 'png', 'pdf']:
        fig.savefig(args.output / ('domain_comparisons.' + suffix), dpi=180)
    plt.close(fig)
    taxa_path = ROOT / 'metadata/analysis_manifest.tsv'
    taxa = {r['taxon_id']: r['species_name'] for r in csv.DictReader(taxa_path.open(), delimiter='\t')}
    ranked = []
    for rank, row in enumerate(sorted(rows, key=lambda r: (-float(r['domain_fit_improvement_angstrom']),
                                 r['marker'], r['taxon_a'], r['taxon_b'], r['pfam_accession'])), 1):
        ranked.append({'manual_review_rank': rank, 'species_a': taxa[row['taxon_a']],
                       'species_b': taxa[row['taxon_b']]} | row)
    with (args.output / 'domain_fit_review_order.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, list(ranked[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(ranked)
    result = {'source_receipt_sha256': sha(args.comparisons / 'receipt.json'),
        'script_sha256': sha(Path(__file__)), 'taxon_manifest_sha256': sha(taxa_path),
        'source_label': args.source_label, 'rasterized_points': dense, 'rows_plotted': len(rows), 'rows_without_local_pairs': int((~finite).sum()),
        'median_domain_sites_global_fit_rmsd_angstrom': float(np.median(x)),
        'median_domain_own_fit_rmsd_angstrom': float(np.median(y)),
        'median_pae10_local_fraction': float(np.nanmedian(conf)),
        'review_order': 'All comparisons ranked by the improvement from refitting identical domain sites; a manual inspection order, not independent discoveries or evidence of motion.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
