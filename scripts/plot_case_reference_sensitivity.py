#!/usr/bin/env python3
"""Plot every verified case/reference sensitivity setting, retaining zero cases."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path('results/experimental_structures/whole-domain-case-reference-robustness-20260927-v1')
    proof = Path('metadata/case_experimental_reference_robustness_completed_20260927.json')
    closure = json.loads(proof.read_text())
    receipt = root / 'receipt.json'
    assert closure['producer_receipt_sha256'] == sha(receipt)
    assert closure['readback_receipt_sha256'] == sha(closure['readback_receipt_path'])
    assert closure['status'] == 'complete_full_case_reference_robustness_readback'
    source = root / 'case_summary.tsv'
    assert sha(source) == json.loads(receipt.read_text())['artifacts'][source.name]
    with source.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 936
    case_fields = ['family', 'gene_a', 'gene_b', 'pfam_accession', 'species_name']
    cases = list(dict.fromkeys(tuple(row[k] for k in case_fields) for row in rows))
    assert len(cases) == 13
    metrics = ['whole', 'domain', 'outside_independent', 'outside_domain_anchored']
    titles = ['Whole protein', 'Domain', 'Outside: independent fit', 'Outside: domain-anchored fit']
    states = ['same_direction', 'opposite_direction', 'variable_or_within_margin']
    colors = ['#247a94', '#cf6c32', '#b9bdc2']
    lookup = {}
    for row in rows:
        key = tuple(row[k] for k in case_fields) + (row['metric'], row['screen'], row['margin_angstrom'])
        assert key not in lookup
        lookup[key] = row
    output = Path('results/figures/case-reference-sensitivity-20260927-v1')
    output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    exported = []
    checked = 0
    with PdfPages(output / 'all_reference_sensitivity_settings.pdf') as pdf:
        for n in [30, 50]:
            for coverage in [50, 70, 90]:
                for margin in ['0.0', '0.01', '0.1']:
                    screen = f'n{n}_c{coverage}'
                    fig, axes = plt.subplots(1, 4, figsize=(19, 8), sharey=True)
                    fig.subplots_adjust(left=.20, right=.995, top=.79, bottom=.23, wspace=.08)
                    fig.suptitle('How experimental reference choice affects duplicate structural contrasts', x=.02, ha='left', y=.98, fontsize=17)
                    fig.text(.02, .925, f'Minimum {n} common residues and {coverage}% regional coverage per fungal protein; direction margin {margin} Å', fontsize=11)
                    fig.text(.02, .884, 'All two-boundary × full/pLDDT ≥70 variants must qualify. Labels: qualified units / entities / dependence components.', fontsize=10)
                    labels = [f'{case[-1]}  ·  {case[0]}' for case in cases]
                    for ax, metric, title in zip(axes, metrics, titles):
                        for y, case in enumerate(cases):
                            row = lookup[case + (metric, screen, margin)]
                            total = int(row['qualified_units'])
                            counts = [int(row[state + '_units']) for state in states]
                            assert sum(counts) == total
                            left = 0.
                            for state, color, count in zip(states, colors, counts):
                                fraction = count / total if total else 0.
                                bar = ax.barh(y, fraction, left=left, height=.62, color=color)[0]
                                assert abs(bar.get_width() - fraction) < 1e-12
                                checked += 1
                                exported.append(dict(zip(case_fields, case), metric=metric, screen=screen, margin_angstrom=margin, category=state, count=count, qualified_units=total, fraction=fraction if total else '', qualified_entities=row['qualified_entities'], qualified_dependency_components=row['qualified_dependency_components']))
                                left += fraction
                            if total:
                                label = f"{total}/{row['qualified_entities']}/{row['qualified_dependency_components']}"
                                ax.text(1.025, y, label, va='center', fontsize=8)
                            else:
                                ax.text(.02, y, 'No qualifying units', va='center', fontsize=8, color='#666666')
                        ax.set_title(title, fontsize=10, loc='left', pad=12)
                        ax.set_xlim(0, 1.45)
                        ax.set_xticks([0, .5, 1], ['0%', '50%', '100%'])
                        ax.set_xlabel('Share of qualified reference units')
                        ax.set_yticks(range(13), labels)
                        ax.set_ylim(12.7, -.7)
                        ax.spines[['top', 'right', 'left']].set_visible(False)
                        ax.tick_params(axis='y', length=0)
                    fig.legend([Patch(color=c) for c in colors], ['Same direction', 'Opposite direction', 'Variable or within margin'], loc='lower left', bbox_to_anchor=(.20, .125), ncol=3, frameon=False)
                    fig.text(.02, .087, 'Directions compare A-minus-B RMSD using fungal versus experimental references on the same sequence-anchored residue maps.', fontsize=10)
                    fig.text(.02, .055, 'Units share chains, sequences and publications; proportions are descriptive, not independent replication, significance or experimental validation.', fontsize=10)
                    fig.text(.02, .023, 'Selected 13 cases; zero coverage is not evidence of no structural difference. Domain orientation and prediction uncertainty remain relevant.', fontsize=10)
                    pdf.savefig(fig)
                    if screen == 'n50_c70' and margin == '0.1':
                        for ext in ['png', 'svg', 'pdf']:
                            fig.savefig(output / ('reference_sensitivity_n50_c70_margin01.' + ext), dpi=140)
                    plt.close(fig)
    table = output / 'plotted_counts.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, list(exported[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(exported)
    with table.open() as handle:
        readback = list(csv.DictReader(handle, delimiter='\t'))
    assert readback == [{k: str(v) for k, v in row.items()} for row in exported]
    assert checked == len(exported) == 2808
    result = dict(status='complete_figure_pending_visual_review',pages=18,case_metric_settings=936,bars_checked=checked,source_hashes={str(p):sha(p) for p in [proof, receipt, source]},script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in output.iterdir()},scope='All 13 cases, four metrics, six screens and three margins; every plotted fraction and zero-denominator label checked. Descriptive dependent reference-unit proportions only.')
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
