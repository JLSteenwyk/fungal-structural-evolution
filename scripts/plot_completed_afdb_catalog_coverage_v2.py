#!/usr/bin/env python3
"""Plot verified full-universe catalog availability without confidence claims."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    args.output.mkdir(exist_ok=False)
    cp = Path('metadata/completed_afdb_catalog_comparison_completed_20261005_v1.json')
    closure = json.loads(cp.read_text())
    assert closure['status'] == 'complete_verified_full_afdb_catalog_comparison'
    pins = {}; bind(pins, cp); bind(pins, Path(__file__))
    source = Path(closure['public_change_table'])
    bind(pins, source, closure['source_hashes'][str(source)])
    with source.open() as handle: rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == len({r['taxon_id'] for r in rows}) == 526
    summary = defaultdict(Counter)
    for r in rows:
        n = int(r['representative_proteins'])
        summary[r['study_role']].update(taxa=1, representative_proteins=n,
            old=int(r['proteins_with_model_old']), new=int(r['proteins_with_model_new']))
    assert summary['ingroup']['taxa'] == 501 and summary['outgroup']['taxa'] == 25
    for role, values in summary.items():
        expected = closure['study_role_change'][role]
        assert values['old'] == expected['proteins_with_model_old']
        assert values['new'] == expected['proteins_with_model_new']
        assert values['representative_proteins'] == expected['representative_proteins']
    all_values = sum(summary.values(), Counter())
    assert all_values['old'] == closure['old_links'] and all_values['new'] == closure['new_links']
    assert all_values['representative_proteins'] == closure['representative_proteins']
    groups = [('All taxa', all_values), ('Fungi', summary['ingroup']), ('Outgroups', summary['outgroup'])]
    table = args.output / 'protein_weighted_coverage.tsv'
    with table.open('x') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['group', 'taxa', 'representative_proteins', 'old_links', 'new_links', 'old_fraction', 'new_fraction'])
        for name, values in groups:
            n = values['representative_proteins']
            writer.writerow([name, values['taxa'], n, values['old'], values['new'], values['old']/n, values['new']/n])
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.3))
    for role, label, color, marker in [('ingroup', '501 fungi', '#2271B2', 'o'), ('outgroup', '25 outgroups', '#D55E00', '^')]:
        subset = [r for r in rows if r['study_role'] == role]
        old = [100*int(r['proteins_with_model_old'])/int(r['representative_proteins']) for r in subset]
        new = [100*int(r['proteins_with_model_new'])/int(r['representative_proteins']) for r in subset]
        axes[0].scatter(old, new, s=22 if role == 'ingroup' else 45, alpha=.6 if role == 'ingroup' else .85,
                        color=color, marker=marker, label=label, edgecolors='none')
    axes[0].plot([0, 100], [0, 100], color='.4', linewidth=1, linestyle='--')
    axes[0].set(xlim=(-3, 103), ylim=(-3, 103), xlabel='September 28 linked proteins (%)',
                ylabel='October 5 linked proteins (%)', title='A. Every taxon retained')
    axes[0].legend(loc='lower right', frameon=False)
    for side, offset, color, label in [('old', -.19, '#999999', 'September 28'), ('new', .19, '#009E73', 'October 5')]:
        values = [100*c[side]/c['representative_proteins'] for _, c in groups]
        bars = axes[1].bar([i+offset for i in range(3)], values, width=.36, color=color, label=label)
        axes[1].bar_label(bars, labels=[f'{v:.1f}%' for v in values], padding=3, fontsize=9)
    axes[1].set(ylim=(0, 100), ylabel='Linked representative proteins (%)',
                title='B. Coverage weighted by protein count')
    axes[1].set_xticks(range(3), [name for name, _ in groups])
    axes[1].legend(frameon=False, loc='upper left')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=.12)
        ax.set_axisbelow(True)
    fig.suptitle('Existing AFDB model availability across 526 study taxa', fontsize=14)
    fig.text(.06, .025, '5,815,847 representative proteins; exact full-length sequence links. All taxa, including zero coverage, retained.\n'
             'Before confidence/PAE filtering and ESMFold integration. Gains are newly linked existing models; no new predictions.', fontsize=9)
    fig.tight_layout(rect=(0, .105, 1, .94))
    public = []
    for extension in ['png', 'pdf']:
        path = args.output / ('catalog_coverage.' + extension)
        fig.savefig(path, dpi=180)
        copy = Path('docs/figures/completed_afdb_catalog_coverage_20261005_v2.' + extension)
        with copy.open('xb') as handle: handle.write(path.read_bytes())
        assert sha(copy) == sha(path)
        bind(pins, path); bind(pins, copy); public.append(str(copy))
    plt.close(fig)
    bind(pins, table); verify(pins)
    result = dict(status='completed_verified_full_afdb_catalog_coverage_figure',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, representative_proteins=5815847,
        taxon_points=526, protein_weighted_groups=3, public_figures=public,
        source_hashes=pins, scientific_eligibility=False, gpu=False, new_predictions=0,
        scope='All 526 source taxon rows plotted, including zeros; pooled coverage divides total links by '
              'total representative proteins per group. Source counts agree with full independent comparison '
              'closure. No confidence, homology, prediction accuracy, evolutionary effect or uncertainty inference.')
    with args.receipt.open('x') as handle: json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__': main()
