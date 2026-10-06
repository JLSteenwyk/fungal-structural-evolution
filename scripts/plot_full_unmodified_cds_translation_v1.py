#!/usr/bin/env python3
"""Descriptive fixed-code diagnostics from closed full translation evidence."""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from reference_measurement_union_sources import bind, verify
from ancestral_chain_attempt import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['completion', 'producer', 'output', 'receipt']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    completion, producer = [json.loads(p.read_text()) for p in [args.completion, args.producer]]
    assert completion['status'] == 'complete_verified_all526_unmodified_cds_fixed_code_translation'
    assert completion['source_hashes'][str(args.producer)] == sha(args.producer)
    assert completion['taxa'] == 526 and completion['target_records'] == 5923039
    assert all(t['actual_terminal_exit_code'] == 0 for t in completion['original_transports'])
    verify(completion['source_hashes'])
    roles = ['inherited', 'snapshot_nuclear', 'snapshot_mitochondrial']
    categories = ['translated_exact_match_to_normalized_protein', 'translated_mismatch_to_normalized_protein',
                  'not_translated_non_triplet_length', 'translated_no_linked_normalized_protein',
                  'unspecified_code_not_tested']
    labels = ['Exact protein match', 'Protein mismatch', 'Nontriplet; not translated',
              'No linked protein', 'Code unspecified; not tested']
    colors = ['#277DA1', '#F3722C', '#777777', '#C4A5CF', '#DDDDDD']
    counts = completion['target_role_status_counts']
    values = np.array([[counts.get(role + ':' + category, 0) for category in categories] for role in roles], dtype=np.int64)
    assert int(values.sum()) == sum(counts.values())
    assert all(int(row.sum()) == 5923039 for row in values)
    taxa = {r['taxon_id']: r for r in producer['taxa_reports']}
    selected = ['F241526', 'F54195', 'F5486']
    species = ['Candida africana', 'Ascoidea rubescens', 'Candida viswanathii']
    context = json.loads(Path('metadata/full_genetic_code_context_completed_20261006_v1.json').read_text())
    context_table = Path('metadata/full_genetic_code_context_taxon_dispositions_20261006_v1.tsv')
    assert context['source_hashes'][str(context_table)] == sha(context_table)
    with context_table.open() as handle:
        rows = {r['taxon_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    source_flags, mismatches = [], []
    for taxon in selected:
        source_flags.append(taxa[taxon]['source_products'])
        assert rows[taxon]['source_products'] == str(taxa[taxon]['source_products'])
        assert taxa[taxon]['source_products'] == taxa[taxon]['target_records']
        assert taxa[taxon]['product_role_status_counts'].get('inherited:no_original_target', 0) == 0
        expected = int(rows[taxon]['products:differs_from_snapshot_nuclear_and_mitochondrial_codes_requires_review'])
        assert expected == source_flags[-1]
        mismatches.append(taxa[taxon]['target_role_status_counts'].get('snapshot_nuclear:translated_mismatch_to_normalized_protein', 0))
    args.output.mkdir(parents=True)
    data_path = args.output / 'diagnostic_counts.tsv'
    with data_path.open('x', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['panel', 'role_or_taxon', 'classification', 'count', 'denominator', 'units'])
        for role, row in zip(roles, values):
            for category, value in zip(categories, row):
                writer.writerow(['fixed_code_targets', role, category, int(value), 5923039, 'original_CDS_targets'])
        for taxon, total, mismatch in zip(selected, source_flags, mismatches):
            writer.writerow(['flagged_entries', taxon, 'inherited_code_differs_from_both_snapshot_assignments', total, total, 'original_products'])
            writer.writerow(['flagged_entries', taxon, 'snapshot_nuclear_translation_mismatch', mismatch, total, 'original_CDS_targets'])
    figure, axes = plt.subplots(2, 1, figsize=(10.8, 7.5), gridspec_kw={'height_ratios': [1.1, 1]})
    start = np.zeros(3)
    for index, (label, color) in enumerate(zip(labels, colors)):
        axes[0].barh(np.arange(3), values[:, index] / 1e6, left=start, color=color, label=label)
        start += values[:, index] / 1e6
    axes[0].set_yticks(np.arange(3), ['Saved inherited code', 'Snapshot nuclear code', 'Snapshot mitochondrial code'])
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 6.1)
    axes[0].set_xlabel('Original CDS targets (millions); same full denominator for each role')
    axes[0].set_title('A  Fixed-code comparisons across all 526 entries', loc='left')
    axes[0].legend(loc='upper center', bbox_to_anchor=(0.55, -0.32), ncol=3, frameon=False, fontsize=8)
    positions = np.arange(3)
    axes[1].barh(positions - 0.16, source_flags, height=0.3, color='#777777', label='Products with an inherited-code assignment disagreement')
    axes[1].barh(positions + 0.16, mismatches, height=0.3, color='#F3722C', label='Original CDSs mismatching the protein under snapshot nuclear code')
    for position, total, mismatch in zip(positions, source_flags, mismatches):
        axes[1].text(total + 120, position - 0.16, f'{total:,}', va='center', fontsize=9)
        axes[1].text(mismatch + 120, position + 0.16, f'{mismatch:,}', va='center', fontsize=9)
    axes[1].set_yticks(positions, species)
    axes[1].invert_yaxis()
    axes[1].set_xlim(0, max(source_flags) * 1.22)
    axes[1].xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
    axes[1].set_xlabel('Records; each selected entry has one original CDS per product')
    axes[1].set_title('B  Code-assignment disagreement does not imply every translation differs', loc='left')
    axes[1].legend(loc='upper center', bbox_to_anchor=(0.6, -0.31), frameon=False, fontsize=8)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    figure.subplots_adjust(left=0.25, right=0.96, top=0.94, bottom=0.17, hspace=1.15)
    figure.text(0.03, 0.015, 'Descriptive diagnostics only: codes are fixed comparisons; neither compartment nor biological code is adopted.\n'
                'Mismatches are not proved annotation errors, evolutionary changes or selection results.', fontsize=9)
    output_paths = [args.output / ('fixed_code_diagnostics.' + extension) for extension in ['png', 'pdf']]
    for path in output_paths:
        figure.savefig(path, dpi=180)
    plt.close(figure)
    pins = {}
    for path, value in completion['source_hashes'].items():
        bind(pins, path, value)
    for path in [args.completion, args.producer, context_table,
                 Path('metadata/full_genetic_code_context_completed_20261006_v1.json'), Path(__file__),
                 data_path, *output_paths]:
        bind(pins, path)
    result = dict(status='complete_descriptive_full_fixed_code_translation_figure',
        checked_utc=datetime.now(timezone.utc).isoformat(), targets=5923039, taxa=526,
        plotted_target_role_status_counts=values.tolist(), plotted_roles=roles, plotted_categories=categories,
        flagged_taxa=selected, inherited_assignment_disagreement_products=source_flags,
        snapshot_nuclear_mismatch_targets=mismatches, source_hashes=pins,
        figures=[str(p) for p in output_paths], data=str(data_path), scientific_eligibility=False,
        genetic_code_admission=False, biological_codon_eligibility=False, gpu=False, new_predictions=0,
        scope='Source-bound descriptive plot of full independently closed fixed-code target totals and three '
              'flagged entry record counts. Original products and CDSs have distinct units, with one-to-one '
              'target counts in these three entries. No error, compartment, functional or evolutionary conclusion.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
