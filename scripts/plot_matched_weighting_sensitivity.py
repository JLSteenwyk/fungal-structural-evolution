"""Plot the complete verified unadjusted matched-domain sensitivity grid."""
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    roots = [Path('results/structural_comparisons/matched-domain-record-summaries-20260927-v1'),
             Path('results/structural_comparisons/matched-record-sensitivity-20260927-v1')]
    proofs = [Path('metadata/matched_domain_record_summary_readback_20260927.json'),
              Path('metadata/matched_record_sensitivity_readback_20260927.json')]
    names = ['record_summary.tsv', 'order_weighting_sensitivity.tsv']
    statuses = ['passed_full_matched_domain_record_summary_readback',
                'passed_full_matched_record_sensitivity_readback']
    sources = {}
    tables = []
    for root, proof, name, status in zip(roots, proofs, names, statuses):
        receipt = json.loads((root / 'receipt.json').read_text())
        audit = json.loads(proof.read_text())
        assert audit['status'] == status and audit['source_receipt_sha256'] == sha(root / 'receipt.json')
        assert sha(root / name) == receipt['artifacts'][name]
        for path in (root / 'receipt.json', proof, root / name):
            sources[str(path)] = sha(path)
        tables.append(pd.read_csv(root / name, sep='\t'))
    raw, sensitivity = tables
    assert len(raw) == 82944 and len(sensitivity) == 62208
    weights = ['record', 'family_equal', 'taxon_equal']
    labels = ['Equal records', 'Equal families', 'Equal taxa']
    quantile_rows = []
    for weight in weights:
        values = raw[f'rmsd_difference_{weight}_mean']
        assert values.notna().all() and np.isfinite(values).all()
        for q in (0., .25, .5, .75, 1.):
            quantile_rows.append(dict(weighting=weight, quantile=q, value=values.quantile(q), settings=len(values)))
    quantiles = pd.DataFrame(quantile_rows)
    scenarios = [f'S{i:02d}' for i in range(1, 55)]
    counts = []
    for weight in weights:
        for scenario in scenarios:
            part = sensitivity[sensitivity.weighting.eq(weight) & sensitivity.scenario_id.eq(scenario)]
            assert len(part) == 384
            assert part.rmsd_order_sign_status.isin(['positive_all_orders', 'negative_all_orders']).all()
            negative = int(part.rmsd_order_sign_status.eq('negative_all_orders').sum())
            counts.append(dict(weighting=weight, scenario_id=scenario,
                               negative_groups=negative, total_groups=len(part), negative_fraction=negative / len(part)))
    counts = pd.DataFrame(counts)
    out = Path('docs/figures/matched_weighting_sensitivity_20260927')
    paths = [out.with_suffix('.' + ext) for ext in ('png', 'pdf', 'svg')]
    qp = Path(str(out) + '_quantiles.tsv')
    cp = Path(str(out) + '_negative_groups.tsv')
    rp = Path(str(out) + '.receipt.json')
    assert not any(p.exists() for p in paths + [qp, cp, rp])
    quantiles.to_csv(qp, sep='\t', index=False)
    counts.to_csv(cp, sep='\t', index=False)
    # Independent scalar order-statistic interpolation from the original TSV.
    values = {w: [] for w in weights}
    for row in csv.DictReader((roots[0] / names[0]).open(), delimiter='\t'):
        for weight in weights:
            values[weight].append(float(row[f'rmsd_difference_{weight}_mean']))
    for weight in weights:
        values[weight].sort()
    for row in csv.DictReader(qp.open(), delimiter='\t'):
        ordered = values[row['weighting']]
        position = float(row['quantile']) * (len(ordered) - 1)
        lo, hi = math.floor(position), math.ceil(position)
        expected = ordered[lo] + (position - lo) * (ordered[hi] - ordered[lo])
        assert math.isclose(float(row['value']), expected, rel_tol=1e-12, abs_tol=1e-14)
        assert int(row['settings']) == len(ordered)
    independent_total, independent_negative = Counter(), Counter()
    for row in csv.DictReader((roots[1] / names[1]).open(), delimiter='\t'):
        key = row['weighting'], row['scenario_id']
        independent_total[key] += 1
        independent_negative[key] += row['rmsd_order_sign_status'] == 'negative_all_orders'
    for row in csv.DictReader(cp.open(), delimiter='\t'):
        key = row['weighting'], row['scenario_id']
        assert int(row['total_groups']) == independent_total[key]
        assert int(row['negative_groups']) == independent_negative[key]
        assert float(row['negative_fraction']) == independent_negative[key] / independent_total[key]
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(2, 1, figsize=(12, 7.2), gridspec_kw={'height_ratios': [1.2, 1]})
    colors = ['#31688e', '#b65b30', '#238b65']
    for i, (weight, color) in enumerate(zip(weights, colors)):
        q = quantiles[quantiles.weighting.eq(weight)].set_index('quantile').value
        axes[0].plot([q[0.], q[1.]], [i, i], lw=2, color=color)
        axes[0].plot([q[.25], q[.75]], [i, i], lw=10, color=color, solid_capstyle='butt')
        axes[0].scatter([q[.5]], [i], color='white', edgecolor=color, s=65, zorder=3)
    axes[0].axvline(0, color='#333333', lw=1, linestyle='--')
    axes[0].set_yticks(range(3), labels)
    axes[0].set_ylim(2.6, -.6)
    axes[0].set_xlabel('Unadjusted RMSD contrast: duplicate pair − matched background (Å)')
    axes[0].set_title('A  Full range, interquartile range and median across all 82,944 settings', loc='left', fontsize=11)
    axes[0].spines[['top', 'right']].set_visible(False)
    axes[0].grid(axis='x', alpha=.15)
    heat = counts.pivot(index='weighting', columns='scenario_id', values='negative_fraction').loc[weights, scenarios]
    im = axes[1].imshow(heat, aspect='auto', cmap='Oranges', vmin=0, vmax=max(.2, float(heat.to_numpy().max())))
    axes[1].set_yticks(range(3), labels)
    axes[1].set_xticks(range(54), [str(i) for i in range(1, 55)], fontsize=7, rotation=90)
    axes[1].set_xlabel('Matching scenario (S01–S54; 384 groups per scenario and weighting)')
    axes[1].set_title('B  Fraction of groups with a negative contrast in every input order', loc='left', fontsize=11)
    colorbar = fig.colorbar(im, ax=axes[1], pad=.015, fraction=.025)
    colorbar.ax.yaxis.set_major_formatter(PercentFormatter(1))
    fig.suptitle('Matched-domain contrasts depend on weighting', fontsize=16, y=.98)
    fig.text(.07, .025, 'All guides, matching policies, masks, domain boundaries and coverage screens are retained. Settings are dependent.\nRanges describe sensitivity, not confidence intervals. Contrasts are unadjusted; phylogenetic models and calibration remain pending.', fontsize=9)
    fig.tight_layout(rect=[0, .105, 1, .955], h_pad=2.2)
    for path in paths:
        fig.savefig(path, dpi=180)
    plt.close(fig)
    for path, expected in sources.items():
        assert sha(path) == expected
    receipt = dict(status='complete_verified_matched_weighting_sensitivity_figure',
        raw_setting_rows=len(raw), sensitivity_rows=len(sensitivity),
        quantile_values_checked=len(quantiles), scenario_weighting_cells_checked=len(counts),
        negative_groups_by_weighting={w: int(counts[counts.weighting.eq(w)].negative_groups.sum()) for w in weights},
        source_sha256=sources, script_sha256=sha(__file__),
        artifacts={str(p): sha(p) for p in paths + [qp, cp]},
        scope='All settings retained; exported quantiles independently reconstructed by scalar interpolation and all negative-group counts by CSV counters. Descriptive weighting sensitivity, not independent replications, confidence intervals, adjusted duplication effects or biological significance.')
    rp.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k not in ('source_sha256', 'artifacts')}, indent=2))


if __name__ == '__main__':
    main()
