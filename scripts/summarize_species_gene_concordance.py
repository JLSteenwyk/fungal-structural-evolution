#!/usr/bin/env python3
"""Join all closed gCF branches to original support and plot every available value.

Descriptive concordance differs from clade support. These summaries do not
accept a species root or assign causes to discordance. Zero-denominator values
remain unavailable; no branch is dropped from the exported table.
"""
import argparse
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--completion', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    completion = json.loads(args.completion.read_text())
    assert completion['status'] == 'complete_verified_full_native_species_gene_concordance'
    assert completion['branch_marker_cells'] == 1046000 and completion['branch_summary_rows'] == 8368
    assert completion['exact_process_journals_checked'] == 2
    archive_path = Path(completion['full_hash_archive'])
    assert sha(archive_path) == completion['full_hash_archive_sha256']
    archive = json.loads(archive_path.read_text())
    assert len(archive['services']) == 2 and archive['summary']['branch_marker_cells'] == 1046000
    bindings = dict(archive['source_hashes'])
    bind(bindings, args.completion)
    bind(bindings, archive_path)
    bind(bindings, args.plan)
    bind(bindings, __file__)
    verify(bindings)
    plan = json.loads(args.plan.read_text())
    assert args.completion == Path(plan['completion']) and args.output == Path(plan['output'])
    assert plan['resources']['gpu'] is False and plan['resources']['paid_resources'] is False
    for path, digest in plan['pins'].items(): bind(bindings, path, digest)
    producer = json.loads(Path(plan['source_plan']).read_text())
    species = json.loads(Path(producer['species_plan']).read_text())
    manifest = list(csv.DictReader(Path(species['manifest']).open(), delimiter='\t'))
    taxa = sorted(row['taxon_id'] for row in manifest)
    index = {taxon: i for i, taxon in enumerate(taxa)}
    universe = (1 << len(taxa)) - 1
    supports = {}
    for label, spec in species['runs'].items():
        path = Path(spec['audit']) / 'branch_support.tsv'
        bind(bindings, path)
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                mask = sum(1 << index[taxon] for taxon in json.loads(row['split_taxa_json']))
                key = label, row['tree'], min(mask, universe ^ mask)
                assert key not in supports
                supports[key] = row
    assert len(supports) == 4184
    reader = Path(completion['independent_readback']).parent
    summary_path = reader / 'independent_branch_concordance.jsonl.gz'
    rows, seen = [], set()
    groups = defaultdict(list)
    with gzip.open(summary_path, 'rt') as handle:
        for line in handle:
            row = json.loads(line)
            key = row['species_run'], row['tree_type'], int(row['canonical_split_mask'], 16)
            support = supports[key]
            identity = key + (row['marker_alignment'],)
            assert identity not in seen
            seen.add(identity)
            n = row['decisive']
            assert n + row['not_decisive'] == 125
            assert row['concordant'] + row['alternative_1'] + row['alternative_2'] + row['residual_discordance'] == n
            value = 100 * row['concordant'] / n if n else None
            sh = float(support['sh_alrt_percent']) if support['sh_alrt_percent'] else None
            boot = float(support['empirical_ufboot_percent'])
            exported = dict(species_run=row['species_run'], tree_type=row['tree_type'], marker_alignment=row['marker_alignment'],
                            canonical_split_id=row['canonical_split_id'], split_taxa_json=json.dumps(row['split_taxa'], separators=(',', ':')),
                            original_sh_alrt=sh, original_empirical_ufboot=boot,
                            decisive_genes=n, not_decisive_genes=row['not_decisive'], concordant_genes=row['concordant'],
                            canonical_nni1_genes=row['alternative_1'], canonical_nni2_genes=row['alternative_2'],
                            residual_discordant_genes=row['residual_discordance'], gcf_percent=value)
            rows.append(exported)
            groups[row['run']].append(exported)
    assert len(rows) == len(seen) == 8368 and len(groups) == 16
    assert set(seen) == {key + (alignment,) for key in supports for alignment in ['profile', 'mafft']}
    args.output.mkdir(parents=True, exist_ok=False)
    table_path = args.output / 'all_branch_support_and_concordance.tsv'
    with table_path.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    summaries = []
    for name, group in groups.items():
        values = [row['gcf_percent'] for row in group if row['gcf_percent'] is not None]
        summaries.append(dict(run=name, branches=len(group), available_gcf=len(values), zero_decisive=len(group) - len(values),
                              median_gcf_percent=statistics.median(values),
                              min_gcf_percent=min(values), max_gcf_percent=max(values),
                              decisive_gene_range=[min(row['decisive_genes'] for row in group), max(row['decisive_genes'] for row in group)],
                              ufboot95_with_available_gcf=sum(row['original_empirical_ufboot'] >= 95 and row['gcf_percent'] is not None for row in group),
                              ufboot95_and_gcf_below50=sum(row['original_empirical_ufboot'] >= 95 and row['gcf_percent'] is not None and row['gcf_percent'] < 50 for row in group)))
    fig, axes = plt.subplots(4, 2, figsize=(9.6, 12), sharex=True, sharey=True, constrained_layout=True)
    colors = dict(profile='#2878B5', mafft='#C86D22')
    for i, label in enumerate(species['runs']):
        for j, kind in enumerate(['ml', 'consensus']):
            ax = axes[i, j]
            for alignment in ['profile', 'mafft']:
                group = groups[label + '-' + kind + '-' + alignment]
                available = [row for row in group if row['gcf_percent'] is not None]
                ax.scatter([row['original_empirical_ufboot'] for row in available],
                           [row['gcf_percent'] for row in available], s=10, alpha=.4,
                           color=colors[alignment], linewidths=0, label=alignment.upper() + ' markers')
            ax.axvline(95, color='#aaaaaa', linewidth=.6, linestyle=':')
            ax.set_title(label.replace('_', ' / ') + ' — ' + ('ML' if kind == 'ml' else 'consensus'), fontsize=10)
            ax.set_xlim(-2, 102)
            ax.set_ylim(-2, 102)
            ax.set_xticks([0, 25, 50, 75, 100])
            ax.set_yticks([0, 25, 50, 75, 100])
            if i == 3: ax.set_xlabel('Original empirical ultrafast bootstrap (%)')
            if j == 0: ax.set_ylabel('Gene concordance factor (%)')
    axes[0, 0].legend(loc='upper left', fontsize=8)
    fig.suptitle('Marker concordance across candidate fungal species trees', fontsize=14)
    paths = [args.output / 'species_gcf_support.pdf', args.output / 'species_gcf_support.png']
    for path in paths: fig.savefig(path, dpi=180)
    plt.close(fig)
    summary_file = args.output / 'run_summaries.json'
    summary_file.write_text(json.dumps(summaries, indent=2) + '\n')
    verify(bindings)
    result = dict(status='complete_descriptive_full_closed_species_gcf_support_summary', branches=8368,
                  native_runs=16, run_summaries=summaries, source_hashes=bindings,
                  artifacts={path.name: sha(path) for path in [table_path, summary_file] + paths},
                  scientific_eligibility=False, scope=plan['scope'])
    with (args.output / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], branches=len(rows), summaries=summaries)), flush=True)


if __name__ == '__main__':
    main()
