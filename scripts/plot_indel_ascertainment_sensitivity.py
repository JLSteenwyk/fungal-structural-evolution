"""Receipt-bound descriptive figure; repeated models are not independent events."""
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    source = Path('results/ancestral/indel-ascertainment-comparison-20260927-v1')
    out = Path('results/figures/indel-ascertainment-sensitivity-20260928-v1')
    receipt = json.loads((source / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_156_ascertainment_probability_comparisons'
    for name, digest in receipt['artifacts'].items():
        assert sha(source / name) == digest, name
    with gzip.open(source / 'character_comparisons.tsv.gz', 'rt') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    with (source / 'encoding_summary.tsv').open() as f:
        summaries = list(csv.DictReader(f, delimiter='\t'))
    assert len(summaries) == receipt['encodings'] == 156
    assert len({r['encoding'] for r in summaries}) == 156
    assert len(rows) == receipt['node_characters']
    keys = {(r['encoding'], r['level'], r['character']) for r in rows}
    assert len(keys) == len(rows)
    groups = defaultdict(list)
    for r in rows:
        p, q = float(r['all_taxa_gap_probability']), float(r['observed_mask_gap_probability'])
        # Retain raw roundoff; do not clip probabilities or conceal excursions.
        tol = 64 * np.finfo(float).eps
        assert np.isfinite([p, q]).all() and -tol <= p <= 1+tol and -tol <= q <= 1+tol
        assert abs(abs(p-q)-float(r['absolute_difference'])) < 1e-14
        assert int(r['map_switch']) == int((p >= .5) != (q >= .5))
        assert int(r['opposed_states_both_at_least_090']) == int((p >= .9 and q <= .1) or (q >= .9 and p <= .1))
        groups[r['encoding']].append(r)
    for s in summaries:
        rs = groups.pop(s['encoding'], [])
        if not rs:
            assert int(s['node_characters']) == 0
            continue
        assert len(rs) == int(s['node_characters'])
        assert sum(int(r['map_switch']) for r in rs) == int(s['map_switches'])
        assert sum(int(r['opposed_states_both_at_least_090']) for r in rs) == int(s['strong_conflicts'])
        assert abs(max(float(r['absolute_difference']) for r in rs)-float(s['maximum_difference'])) < 1e-14
    assert not groups
    assert sum(int(r['map_switch']) for r in rows) == receipt['map_switches']
    assert sum(int(r['opposed_states_both_at_least_090']) for r in rows) == receipt['strong_conflicts']
    assert abs(max(float(r['absolute_difference']) for r in rows)-receipt['maximum_probability_difference']) < 1e-14
    families = defaultdict(lambda: [0, 0, 0.])
    for r in rows:
        a = families[r['encoding'].split('-')[0]]
        a[0] += 1
        a[1] += int(r['map_switch'])
        a[2] = max(a[2], float(r['absolute_difference']))
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'family_summary.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['family', 'node_character_comparisons', 'map_switches', 'maximum_difference'])
        for family, values in sorted(families.items()):
            w.writerow([family, *values])
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), layout='constrained')
    for policy, color in [('terminal_gap', '#137c8b'), ('terminal_unknown', '#bc5a21')]:
        rs = [r for r in rows if r['encoding'].endswith(policy)]
        axes[0].scatter([float(r['all_taxa_gap_probability']) for r in rs],
                        [float(r['observed_mask_gap_probability']) for r in rs],
                        s=4, alpha=.25, color=color, rasterized=True, label=policy.replace('_', ' '))
    axes[0].plot([0, 1], [0, 1], color='black', lw=.7)
    axes[0].set(xlabel='All-taxa conditioning: P(exact gap)', ylabel='Observed-mask conditioning: P(exact gap)',
                xlim=(-.02, 1.02), ylim=(-.02, 1.02), title='A  Conditional ancestral probabilities', aspect='equal')
    axes[0].legend(markerscale=3, frameon=False)
    names = sorted(families)
    rates = [100 * families[n][1] / families[n][0] for n in names]
    axes[1].barh(names, rates, color='#137c8b')
    for i, n in enumerate(names):
        a = families[n]
        axes[1].text(rates[i]+.12, i, f'{a[1]}/{a[0]:,}', va='center', fontsize=8)
    axes[1].set(xlabel='Comparisons switching at P = 0.5 (%)', xlim=(0, max(rates)+5),
                title='B  Sensitivity by protein family')
    axes[1].invert_yaxis()
    fig.suptitle('Indel ascertainment sensitivity — fixed fitted parameters\nRepeated nodes and model settings; not independent evolutionary events', fontsize=11)
    for ext in ['png', 'pdf']:
        fig.savefig(out / f'indel_ascertainment.{ext}', dpi=180)
    plt.close(fig)
    record = {
        'status': 'all_character_rows_and_encoding_aggregates_verified_figure_produced',
        'source_receipt_sha256': sha(source / 'receipt.json'), 'script_sha256': sha(__file__),
        'node_character_comparisons': len(rows), 'map_switches': receipt['map_switches'],
        'strong_conflicts': receipt['strong_conflicts'], 'families': len(families),
        'raw_probability_bounds_excursions': sum(float(r[k]) < 0 or float(r[k]) > 1 for r in rows for k in ['all_taxa_gap_probability', 'observed_mask_gap_probability']),
        'bounds_tolerance': float(64 * np.finfo(float).eps),
        'artifacts': {str(p): sha(p) for p in sorted(out.iterdir()) if p.suffix in ['.png', '.pdf', '.tsv']},
        'scope': 'Descriptive sensitivity only. Three empty encodings retained in source. No model ranking, independent-event counts, parameter uncertainty or compatible ancestral sequence claim.'}
    (out / 'receipt.json').write_text(json.dumps(record, indent=2)+'\n')
    Path('metadata/indel_ascertainment_figure_20260928.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k != 'artifacts'}, indent=2))


if __name__ == '__main__':
    main()
