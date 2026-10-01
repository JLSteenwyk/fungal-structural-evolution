#!/usr/bin/env python3
"""Publish fully bound crossed-run RF tables and a standalone checked figure."""
import argparse
import csv
import json
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--completion', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text()); c = json.loads(args.completion.read_text()); root = Path(plan['output'])
    assert c['status'] == 'complete_verified_full_four_run_pmsf_ML_consensus_sensitivity' and len(c['services']) == 2
    for path, digest in c['source_hashes'].items(): assert sha(path) == digest
    assert c['source_hashes'][str(args.plan)] == sha(args.plan)
    with (root / 'comparisons.tsv').open() as f: comparisons = list(csv.DictReader(f, delimiter='\t'))
    with (root / 'split_presence.tsv').open() as f: presence = list(csv.DictReader(f, delimiter='\t'))
    splits = {}
    for row in presence:
        if row['present'] == 'True': splits.setdefault(row['view'], set()).add(tuple(json.loads(row['split_taxa_json'])))
    assert len(splits) == 8 and all(len(s) == 523 for s in splits.values())
    runs = list(plan['runs']); tags = ['P/P', 'P/M', 'M/P', 'M/M']; matrices = {}; within = np.zeros(4, dtype=int)
    indexed = {(r['view_a'], r['view_b']): r for r in comparisons}; assert len(indexed) == 16
    for kind in ['ml', 'consensus']:
        matrix = np.zeros((4, 4), dtype=int)
        for i in range(4):
            for j in range(i + 1, 4):
                a, b = runs[i] + ':' + kind, runs[j] + ':' + kind
                distance = len(splits[a].symmetric_difference(splits[b]))
                assert int(indexed[a, b]['rf_distance']) == distance
                matrix[i, j] = matrix[j, i] = distance
        matrices[kind] = matrix
    for i, run in enumerate(runs):
        a, b = run + ':ml', run + ':consensus'; within[i] = len(splits[a] ^ splits[b])
        assert int(indexed[a, b]['rf_distance']) == within[i]
    plt.rcParams['svg.fonttype'] = 'none'; plt.rcParams['font.size'] = 10
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), gridspec_kw={'width_ratios': [1, 1, .8]}, layout='constrained')
    vmax = max(int(m.max()) for m in matrices.values()); labels = {}
    for ax, kind, title in zip(axes[:2], ['ml', 'consensus'], ['Maximum-likelihood trees', 'Bootstrap-consensus trees']):
        mat = matrices[kind]; im = ax.imshow(mat, vmin=0, vmax=vmax, cmap='Blues')
        ax.set_xticks(range(4), tags); ax.set_yticks(range(4), tags); ax.set_title(title); ax.set_xlabel('Alignment / guide')
        for i in range(4):
            for j in range(4):
                gid = f'rf-{kind}-{i}-{j}'; label = ax.text(j, i, str(mat[i, j]), ha='center', va='center', color='white' if mat[i, j] > vmax * .55 else 'black'); label.set_gid(gid); labels[gid] = str(mat[i, j])
        assert np.array_equal(im.get_array(), mat)
    colorbar = fig.colorbar(im, ax=list(axes[:2]), shrink=.75); colorbar.set_label('RF distance (internal splits)')
    bars = axes[2].bar(tags, within, color='#4d779e'); axes[2].set_ylim(0, max(12, int(within.max()) + 3)); axes[2].set_title('ML versus consensus'); axes[2].set_ylabel('RF distance'); axes[2].set_xlabel('Alignment / guide')
    for i, bar in enumerate(bars):
        assert bar.get_height() == int(within[i]); label = axes[2].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + .25, str(within[i]), ha='center'); gid = f'rf-within-{i}'; label.set_gid(gid); labels[gid] = str(within[i])
    fig.suptitle('Full crossed PMSF sensitivity: 526 taxa in every tree\nP = profile; M = MAFFT. Distances describe topology differences, not evolutionary rates.', fontsize=12)
    folder = Path('docs/figures'); folder.mkdir(exist_ok=True); png = folder / 'pmsf_four_run_RF_sensitivity_20261001.png'; svg = folder / 'pmsf_four_run_RF_sensitivity_20261001.svg'
    assert not png.exists() and not svg.exists(); fig.savefig(png, dpi=200); fig.savefig(svg); plt.close(fig)
    exported = ET.parse(svg).getroot()
    for gid, value in labels.items():
        nodes = exported.findall(f'.//*[@id="{gid}"]'); assert len(nodes) == 1
        assert ''.join(nodes[0].itertext()).strip() == value
    table = Path('docs/tables/pmsf_four_run_split_comparisons_20261001.tsv'); boundary = Path('docs/tables/pmsf_four_run_role_boundary_20261001.tsv')
    assert not table.exists() and not boundary.exists(); shutil.copyfile(root / 'comparisons.tsv', table); shutil.copyfile(root / 'rooting_boundary.tsv', boundary)
    assert sha(table) == c['source_hashes'][str(root / 'comparisons.tsv')] and sha(boundary) == c['source_hashes'][str(root / 'rooting_boundary.tsv')]
    result = dict(status='passed_bound_full_four_run_pmsf_RF_figure_and_tables', plan_sha256=sha(args.plan), completion_sha256=sha(args.completion),
                  matrices={k: v.tolist() for k, v in matrices.items()}, within_run_RF=within.tolist(), checked_SVG_value_labels=len(labels), comparison_table_rows=16, boundary_table_rows=8,
                  artifacts={str(p): sha(p) for p in [png, svg, table, boundary]}, script_sha256=sha(__file__), scientific_eligibility=False,
                  scope='Full source closure rechecked. All12cross-run and4within-run RF distances independently recounted from the entire verified split-presence grid; both plotted matrices,4bar heights and36exported SVG value labels checked. Standalone PNG/SVG and full16/8row tables bound. No pixelwise raster validation, inferential significance, chosenroot or accepted species framework.')
    with args.output.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
