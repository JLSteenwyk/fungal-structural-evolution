#!/usr/bin/env python3
"""Summarize observed broad-eukaryotic marker recovery across the full sample."""
import csv
import hashlib
import json
import statistics
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    manifest = ROOT / 'metadata/analysis_manifest.tsv'
    qc_path = ROOT / 'metadata/busco_eukaryota_qc.tsv'
    taxa = read_table(manifest)
    qc = {r['taxon_id']: r for r in read_table(qc_path)}
    if set(qc) != {r['taxon_id'] for r in taxa}:
        raise ValueError('QC must cover precisely the complete analysis manifest')
    groups = defaultdict(list)
    for taxon in taxa:
        row = qc[taxon['taxon_id']]
        group = taxon['lineage'].split(';')[0] if taxon['study_role'] == 'ingroup' else 'Outgroups'
        groups[group].append((taxon, row))
    names = sorted(g for g in groups if g != 'Outgroups') + ['Outgroups']
    summary = []
    for group in names:
        values = [float(row['complete_percent']) for _, row in groups[group]]
        summary.append({'group': group, 'taxa': len(values), 'minimum': min(values),
                        'median': statistics.median(values), 'maximum': max(values)})
    out = ROOT / 'results/qc'
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'marker_recovery_by_group.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, list(summary[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(summary)
    fig, axes = plt.subplots(1, 2, figsize=(12, 9), sharey=True,
                             gridspec_kw={'width_ratios': [1, 2.5]})
    positions = list(range(len(names)))
    colors = ['#286c8e' if g != 'Outgroups' else '#ba6a2f' for g in names]
    axes[0].barh(positions, [len(groups[g]) for g in names], color=colors)
    axes[0].set_yticks(positions, names)
    axes[0].invert_yaxis()
    axes[0].set_xlabel('Sampled taxa')
    axes[0].set_xlim(0, max(len(v) for v in groups.values()) * 1.18)
    for i, group in enumerate(names):
        axes[0].text(len(groups[group]) + 2, i, str(len(groups[group])), va='center', fontsize=9)
        records = sorted(groups[group], key=lambda pair: pair[0]['taxon_id'])
        values = [float(row['complete_percent']) for _, row in records]
        # Deterministic jitter separates observations without implying phylogenetic distances.
        jitter = [(int(hashlib.sha256(t['taxon_id'].encode()).hexdigest()[:8], 16) / 0xffffffff - .5) * .48
                  for t, _ in records]
        axes[1].scatter(values, [i + j for j in jitter], s=15, alpha=.55, color=colors[i], linewidths=0)
        axes[1].plot(statistics.median(values), i, marker='|', color='black', markersize=14, markeredgewidth=2)
    axes[1].set_xlim(-2, 102)
    axes[1].set_xlabel('Complete eukaryota_odb12.2 markers (%)')
    axes[1].grid(axis='x', alpha=.2)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Full-dataset sampling and broad marker recovery\n501 fungi + 25 outgroups; raw proteomes; 125 markers', fontsize=14)
    fig.text(.5, .025, 'Each point is one taxon; black ticks show group medians.\nLow recovery can reflect marker loss, divergence or annotation limitations; no taxa are excluded here.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .065, 1, .94))
    for suffix in ['png', 'pdf', 'svg']:
        fig.savefig(out / f'marker_recovery.{suffix}', dpi=180)
    plt.close(fig)
    hashes = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    (ROOT / 'metadata/full_dataset_qc_summary.json').write_text(json.dumps({
        'taxa': len(taxa), 'groups': summary, 'lineage_dataset': 'eukaryota_odb12.2',
        'manifest_sha256': hashes(manifest), 'qc_sha256': hashes(qc_path),
        'matplotlib_version': matplotlib.__version__,
        'artifacts': [{'path': str(p.relative_to(ROOT)), 'sha256': hashes(p)} for p in sorted(out.glob('marker_recovery*'))],
        'interpretation': 'Descriptive raw-proteome marker recovery, not an automatic exclusion rule or independent genome completeness estimate.'}, indent=2) + '\n')
    print('Full-dataset QC figure and summary written for', len(taxa), 'taxa')


if __name__ == '__main__':
    main()
