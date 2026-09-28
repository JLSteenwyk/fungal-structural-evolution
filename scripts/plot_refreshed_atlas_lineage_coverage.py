"""Describe frozen AFDB coverage by manifest lineage; no inferential tests."""
import csv
import json
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from ancestral_chain_attempt import write_json
from readback_whole_proteome_catalog import sha


def read(path):
    with open(path) as h:
        return list(csv.DictReader(h, delimiter='\t'))


def write(path, rows):
    with path.open('w') as h:
        writer = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    assert read(path) == [{k: str(v) for k, v in row.items()} for row in rows]


def main():
    base = Path('results/structures/whole-proteome-catalog-change-20260928-v1')
    receipt = json.loads((base / 'receipt.json').read_text())
    audit_path = Path('metadata/whole_proteome_catalog_comparison_readback_20260928.json')
    audit = json.loads(audit_path.read_text())
    assert audit['status'] == 'passed_full_catalog_comparison_source_replay'
    assert audit['source_receipt_sha256'] == sha(base / 'receipt.json')
    source = base / 'taxon_coverage_change.tsv'
    assert sha(source) == receipt['artifacts'][source.name]
    manifest_path = Path('metadata/analysis_manifest.tsv')
    catalog_plan = json.loads(Path('metadata/whole_proteome_structure_catalog_plan_20260928.json').read_text())
    assert sha(manifest_path) == catalog_plan['pins'][str(manifest_path)]
    manifest = {r['taxon_id']: r for r in read(manifest_path)}
    taxa, groups = [], defaultdict(list)
    for row in read(source):
        m = manifest[row['taxon_id']]
        assert row['species_name'] == m['species_name'] and row['study_role'] == m['study_role']
        lineage = m['lineage'].split(';')[0]
        assert lineage
        row = dict(row, manifest_lineage=lineage)
        taxa.append(row); groups[row['study_role'], lineage].append(row)
    assert len(taxa) == len({r['taxon_id'] for r in taxa}) == 526
    summaries = []
    for (role, lineage), rows in sorted(groups.items()):
        total = sum(int(r['representative_proteins']) for r in rows)
        result = dict(study_role=role, manifest_lineage=lineage, taxa=len(rows), proteins=total)
        for side in ('old', 'new'):
            linked = sum(int(r['proteins_with_model_' + side]) for r in rows)
            values = [int(r['proteins_with_model_' + side]) / int(r['representative_proteins']) for r in rows]
            result.update({side + '_linked_proteins': linked,
                           side + '_protein_weighted_fraction': linked / total,
                           side + '_taxon_median_fraction': float(np.median(values)),
                           side + '_taxon_min_fraction': min(values),
                           side + '_taxon_max_fraction': max(values),
                           side + '_zero_coverage_taxa': sum(v == 0 for v in values)})
        result['newly_linked_proteins'] = sum(int(r['new_catalog_link']) for r in rows)
        summaries.append(result)
    assert sum(r['proteins'] for r in summaries) == 5815847
    assert sum(r['new_linked_proteins'] for r in summaries) == 1955694
    assert sum(r['newly_linked_proteins'] for r in summaries) == 636181
    out = Path('results/figures/refreshed-atlas-lineage-coverage-20260928-v2')
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'taxon_coverage.tsv', taxa)
    write(out / 'lineage_coverage.tsv', summaries)
    fig, axes = plt.subplots(1, 2, figsize=(13, 11), sharey=True)
    colors = ['#9a9a9a', '#197b8e']
    labels = []
    for i, row in enumerate(summaries):
        labels.append(f"{row['manifest_lineage']} (n={row['taxa']})" + (' *' if row['study_role'] == 'outgroup' else ''))
        old, new = [100 * row[s + '_protein_weighted_fraction'] for s in ('old', 'new')]
        axes[0].plot([old, new], [i, i], color='#bbbbbb', lw=1.2)
        for value, color, label, marker in zip([old, new], colors, ['September 22', 'September 28'], ['o', 's']):
            axes[0].scatter(value, i, color=color, marker=marker, s=25, label=label if i == 0 else None, zorder=3)
        members = sorted(groups[row['study_role'], row['manifest_lineage']], key=lambda x: x['taxon_id'])
        offsets = np.linspace(-.23, .23, len(members)) if len(members) > 1 else [0]
        values = [100 * float(r['coverage_fraction_new']) for r in members]
        axes[1].scatter(values, i + np.array(offsets), s=9, color=colors[1], alpha=.45, linewidths=0)
        axes[1].scatter(100 * row['new_taxon_median_fraction'], i, s=40, marker='|', color='#222222', zorder=4)
    axes[0].set_yticks(range(len(labels)), labels, fontsize=9)
    axes[0].set_ylim(len(labels) - .3, -.7)
    for ax in axes:
        ax.set_xlim(-2, 102); ax.set_xticks([0, 25, 50, 75, 100])
        ax.grid(axis='x', color='#eeeeee'); ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_title('Protein-weighted coverage by lineage')
    axes[1].set_title('September 28 coverage of each taxon')
    axes[0].set_xlabel('Proteins with an exact-sequence model (%)')
    axes[1].set_xlabel('Proteins with an exact-sequence model (%)')
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc='lower center', bbox_to_anchor=(.5, .08), ncol=2, frameon=False, fontsize=9)
    fig.suptitle('Structural coverage remains uneven across the sampled lineages', fontsize=14)
    fig.text(.5, .025, 'Frozen AFDB catalogs; local ESMFold models excluded. Points on right: all 526 taxa; black ticks: taxon medians.\nLineage labels follow the project manifest; * denotes outgroups. Taxon counts are entries, not verified unique species.\nCatalog availability precedes confidence filtering and is not evidence of biological structural conservation.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .12, 1, .97))
    fig.savefig(out / 'lineage_coverage.png', dpi=180)
    fig.savefig(out / 'lineage_coverage.pdf'); plt.close(fig)
    result = dict(status='completed_lineage_coverage_figure_pending_visual_inspection',
                  taxa=526, lineage_groups=len(summaries),
                  zero_coverage_taxa=sum(r['new_zero_coverage_taxa'] for r in summaries),
                  sources={str(p): sha(p) for p in [source, base / 'receipt.json', audit_path, manifest_path]},
                  script_sha256=sha(__file__),
                  artifacts={p.name: sha(p) for p in out.iterdir()},
                  scope='Descriptive full-sampling coverage, protein-weighted fractions and taxon distributions. Manifest categories are not independent evolutionary contrasts; no tests or confidence intervals.')
    write_json(out / 'receipt.json', result)
    write_json(Path('metadata/refreshed_atlas_lineage_coverage_figure_20260928.json'), result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
