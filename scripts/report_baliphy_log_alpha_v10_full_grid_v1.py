#!/usr/bin/env python3
"""Publish all closed V10 diagnostics without accepting short-run posteriors."""
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ancestral_chain_attempt import sha
from readback_baliphy_log_alpha_v10 import HEADER, load
from reference_measurement_union_sources import bind, verify


def table(path, rows):
    with path.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    with path.open() as handle:
        observed = list(csv.DictReader(handle, delimiter='\t'))
    assert observed == [{k: str(v) for k, v in row.items()} for row in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'producer', 'reader', 'reader-transport', 'output', 'receipt'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    plan, producer, reader, transport = [json.loads(p.read_text()) for p in
        (args.plan, args.producer, args.reader, args.reader_transport)]
    assert producer['status'] == 'complete_all24_V10_full_comparison_dispositions_pending_readback'
    assert reader['status'] == 'complete_all24_V10_independent_paired_files_latent_states_and_native_readback'
    assert transport['validation_sha256'] == sha(args.reader)
    assert transport['original_tool_terminal_exit_code'] == 0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    assert reader['full_input_logger_noninterference_passed']
    assert reader['full_latent_trace_integrity_passed']
    assert reader['roles'] == reader['complete_diagnostic_roles'] == 24
    assert reader['diagnostic_rows'] == reader['checked_rows'] == 504
    assert reader['paired_file_comparisons'] == 144
    assert reader['scientific_eligibility'] is reader['posterior_qualified'] is False
    verify(transport['source_hashes'])
    pins = dict(transport['source_hashes'])
    root = Path(plan['output'])
    dispositions_path = root / 'dispositions.json'
    pairs_path = root / 'paired_scientific_outputs.json'
    dispositions = {r['chain_id']: r for r in json.loads(dispositions_path.read_text())}
    pairs = {r['chain_id']: r for r in json.loads(pairs_path.read_text())}
    audits = {r['chain_id']: r for r in reader['diagnostic_audits']}
    assert set(dispositions) == set(pairs) == set(audits) and len(pairs) == 24
    records, summaries, plotted = [], [], {}
    for cid, disposition in sorted(dispositions.items()):
        pair, audit = pairs[cid], audits[cid]
        assert pair['all_original_files_byte_identical'] and disposition['exit_code'] == 0
        assert audit['complete_native_trace_checked'] and not audit['errors']
        diagnostic = Path(pair['diagnostic_path'])
        bind(pins, diagnostic, pair['diagnostic_sha256'])
        assert sha(diagnostic) == pair['diagnostic_sha256']
        values = [load(line) for line in diagnostic.read_text().splitlines()]
        assert values[0] == HEADER
        values = values[1:]
        assert [r['iter'] for r in values] == list(range(21))
        overflow_iterations = []
        latent = []
        for row in values:
            state = row['parameters//']['S1/']
            x = Decimal(state['latentLogAlpha'])
            assert x.is_finite()
            latent.append(x)
            special = isinstance(state['derivedAlpha'], str)
            if special:
                assert state['derivedAlpha'] == '__project_scalar_v6__:positive_infinity'
                overflow_iterations.append(row['iter'])
            records.append(dict(chain_id=cid, effective_input_group=disposition['effective_input_group'],
                prior_label=disposition['prior_label'], chain_role=disposition['chain_role'],
                iteration=row['iter'], latent_log_alpha=str(x),
                derived_alpha_state='positive_infinity' if special else 'finite',
                derived_alpha=str(state['derivedAlpha']),
                latent_log_density=str(state['latentLogDensity']),
                category_rates_json=json.dumps([str(r) for r in state['categoryRates']], separators=(',', ':')),
                original_integrity_disposition=disposition['status'], scientific_eligibility=False,
                posterior_qualified=False))
        assert len(overflow_iterations) == audit['overflow_rows']
        assert overflow_iterations == [r['iteration'] for r in audit['overflow_observations']]
        summaries.append(dict(chain_id=cid, effective_input_group=disposition['effective_input_group'],
            prior_label=disposition['prior_label'], chain_role=disposition['chain_role'],
            native_exit_code=disposition['exit_code'], scientific_files_byte_identical=pair['all_original_files_byte_identical'],
            diagnostic_rows=len(values), independently_checked_rows=audit['checked_rows'],
            overflow_rows=len(overflow_iterations), overflow_iterations_json=json.dumps(overflow_iterations, separators=(',', ':')),
            minimum_latent_log_alpha=str(min(latent)), maximum_latent_log_alpha=str(max(latent)),
            original_integrity_disposition=disposition['status'], scientific_eligibility=False, posterior_qualified=False))
        plotted[(disposition['effective_input_group'], disposition['prior_label'], disposition['chain_role'])] = (latent, overflow_iterations)
    assert len(records) == 504 and sum(r['overflow_rows'] for r in summaries) == reader['overflow_rows']
    args.output.mkdir(parents=True)
    role_table = args.output / 'role_dispositions.tsv'
    row_table = args.output / 'all_latent_rows.tsv'
    table(role_table, summaries)
    table(row_table, records)
    inputs = sorted({r['effective_input_group'] for r in summaries})
    priors = ['broad', 'centered', 'package']
    assert len(inputs) == 2 and set(plotted) == {(i, p, c) for i in inputs for p in priors for c in range(1, 5)}
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True, sharey=True, constrained_layout=True)
    colors = ['#277DA1', '#43AA8B', '#F8961E', '#9B5DE5']
    for i, group in enumerate(inputs):
        for j, prior in enumerate(priors):
            ax = axes[i, j]
            for chain, color in zip(range(1, 5), colors):
                latent, special = plotted[group, prior, chain]
                ys = [float(x) for x in latent]
                ax.plot(range(21), ys, color=color, lw=1.3, label='Chain ' + str(chain))
                if special:
                    ax.scatter(special, [ys[k] for k in special], marker='x', s=48, color='#C1121F', zorder=5)
            ax.set_yscale('symlog', linthresh=10)
            ax.grid(alpha=.2)
            ax.set_title(prior.capitalize() + ' prior; input ' + group[:8])
            ax.set_xticks([0, 5, 10, 15, 20])
            if i == 1:
                ax.set_xlabel('Saved diagnostic iteration')
            if j == 0:
                ax.set_ylabel('Latent log(alpha), symmetric log scale')
    axes[0, 2].legend(frameon=False, fontsize=8)
    fig.suptitle('Full 24-run latent-alpha diagnostics: all saved states retained', fontsize=13)
    fig.supxlabel('Red crosses: derived alpha exceeds finite-double representation; finite latent state remains visible.\n'
        'Short-run diagnostics, not posterior estimates. Inputs are protein-tip sets; chains are not biological replicates.', fontsize=9)
    figures = []
    for suffix in ('png', 'pdf'):
        path = args.output / ('latent_alpha_diagnostics.' + suffix)
        fig.savefig(path, dpi=180)
        figures.append(str(path))
        bind(pins, path)
    plt.close(fig)
    for path in (args.plan, args.producer, args.reader, args.reader_transport, dispositions_path, pairs_path,
                 role_table, row_table, Path(__file__)):
        bind(pins, path)
    result = dict(status='complete_closed_all24_V10_diagnostic_tables_and_figure',
        checked_utc=datetime.now(timezone.utc).isoformat(), roles=24, effective_inputs=2,
        paired_scientific_file_comparisons=144, diagnostic_rows=504,
        overflow_rows=reader['overflow_rows'], roles_with_overflow=sum(r['overflow_rows'] > 0 for r in summaries),
        native_zero_roles=producer['native_zero_exit_roles'], integrity_checked_roles=producer['integrity_checked_roles'],
        special_value_review_roles=24-producer['integrity_checked_roles'], figures=figures,
        role_table=str(role_table), all_row_table=str(row_table), source_hashes=pins,
        scientific_eligibility=False, posterior_qualified=False, old_review_arrays_admitted=False,
        historical_latent_values_recovered=False, gpu=False, new_predictions=0, new_native_runs=0,
        all_eight_aims_incomplete=True,
        scope='All 24 closed paired native outcomes, all 144 unchanged scientific files, and all 504 independently '
              'verified latent diagnostics retained in readable tables and a descriptive plot. Finite latent states '
              'and overflow tags explain these newly captured states only. Original special-value reviews remain '
              'excluded; no burn-in chosen, posterior admitted, historical latent recovered or biological effect inferred.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
