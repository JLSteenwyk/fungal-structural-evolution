#!/usr/bin/env python3
"""Publish complete closed order/context counts and a descriptive coverage figure."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from full_triad_context_geometry_sources import DEFINITIONS, MASKS, POLICIES
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def closed(path, status):
    c = json.loads(Path(path).read_text()); assert c['status'] == status and c['exact_process_journals_checked'] == 2
    bindings = {str(path): sha(path)}
    for field in ['producer_receipt', 'independent_readback', 'full_hash_archive']: bind(bindings, c[field], c[field + '_sha256'])
    r, a, archive = [json.loads(Path(c[field]).read_text()) for field in ['producer_receipt', 'independent_readback', 'full_hash_archive']]
    assert a['producer_receipt_sha256'] == sha(c['producer_receipt']) and len(archive['services']) == 2
    assert len(archive['source_hashes']) == c['bound_source_hashes']
    for name, digest in r['artifacts'].items(): bind(bindings, Path(c['producer_receipt']).parent / name, digest)
    return c, r, bindings


def main():
    rp = 'metadata/full_triad_order_robustness_completed_20261001.json'; cp = 'metadata/full_triad_context_geometry_completed_20261001.json'
    robust, _, bindings = closed(rp, 'complete_verified_full_triad_all_order_robustness')
    context, cr, other = closed(cp, 'complete_verified_full_triad_context_geometry'); bindings.update(other); bind(bindings, __file__)
    assert robust['correspondence_work_triads'] == context['measured_triads'] == 27056
    assert robust['robustness_groups'] == context['measured_robustness_groups'] == 108224
    table = Path(context['producer_receipt']).parent / 'context_geometry_counts.tsv'
    with table.open() as handle: rows = list(csv.DictReader(handle, delimiter='\t')); columns = list(rows[0])
    assert len(rows) == context['summary_rows'] == 720
    index = {}
    for row in rows:
        key = tuple(row[k] for k in ['guide', 'design', 'mask', 'mapping_definition', 'screen', 'policy']); assert key not in index
        index[key] = row
        assert 0 <= int(row['passed_contexts']) <= int(row['parent_eligible_contexts']) <= int(row['source_contexts'])
        assert int(row['source_contexts']) == context['guide_contexts'][row['guide']]
    screens = [s['id'] for s in cr['screens']]
    for g in ['profile', 'mafft']:
        for d in ['availability', 'sequence_first']:
            for definition in DEFINITIONS:
                for sid in screens:
                    full, p70, both = [index[g, d, m, definition, sid, POLICIES[2]] for m in MASKS]
                    assert int(both['passed_contexts']) <= min(int(full['passed_contexts']), int(p70['passed_contexts']))
                    for mask in MASKS:
                        q = [int(index[g, d, mask, definition, sid, p]['passed_contexts']) for p in POLICIES]
                        assert q[2] <= q[1] <= q[0] and q[4] <= q[2] <= q[3]
    dest_table = Path('docs/tables/full_triad_context_geometry_counts_20261001.tsv'); assert not dest_table.exists(); dest_table.write_bytes(table.read_bytes())
    order_table = Path('docs/tables/full_triad_all_order_screen_counts_20261001.tsv'); assert not order_table.exists()
    order_rows = []
    for mask in MASKS[:2]:
        for definition in DEFINITIONS:
            for sid in screens:
                key = ':'.join([mask, definition, sid]); all_count = robust['all_order_screen_pass_counts'][key]; any_count = robust['any_order_screen_pass_counts'][key]
                assert 0 <= all_count <= any_count <= robust['correspondence_work_triads']
                order_rows.append(dict(mask=mask, mapping_definition=definition, screen=sid, ordered_physical_triads=robust['correspondence_work_triads'], all_eight_orders_pass=all_count, any_order_pass=any_count))
    with order_table.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(order_rows[0]), delimiter='\t'); writer.writeheader(); writer.writerows(order_rows)
    cohorts = [('profile', 'availability'), ('mafft', 'availability'), ('profile', 'sequence_first'), ('mafft', 'sequence_first')]
    labels = ['Profile\navailability', 'MAFFT\navailability', 'Profile\nsequence first', 'MAFFT\nsequence first']
    colors = ['#377eb8', '#4daf4a', '#984ea3']; x = np.arange(len(cohorts)); width = .24
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True); series = []; checked_bars = 0
    for ax, definition in zip(axes, DEFINITIONS):
        for i, mask in enumerate(MASKS):
            values = [int(index[g, d, mask, definition, 'n50_c70', POLICIES[2]]['passed_contexts']) for g, d in cohorts]
            bars = ax.bar(x + (i - 1) * width, values, width, label={'full': 'Full structure', 'plddt70': 'pLDDT ≥70', 'both_masks': 'Both masks'}[mask], color=colors[i])
            for bar, value in zip(bars, values):
                assert bar.get_height() == value; checked_bars += 1
                ax.annotate(f'{value:,}', (bar.get_x() + bar.get_width() / 2, value), xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=7)
            series.append(dict(mapping_definition=definition, mask=mask, counts=values))
        ax.set_xticks(x, labels); ax.set_title(definition.replace('_', ' ').capitalize()); ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=.2); ax.set_axisbelow(True)
    axes[0].set_ylabel('Source contexts passing all eight alignment orders')
    axes[1].legend(loc='upper right', fontsize=8)
    fig.suptitle('Duplication/reference geometry coverage: fixed lexical reference, both native guides')
    fig.text(.5, .015, '≥50 shared residues and ≥70% of each original protein; original parent gates retained.\nGuide/design/core counts overlap and are not independent evolutionary events.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .065, 1, .95)); prefix = Path('docs/figures/full_triad_context_geometry_coverage_20261001')
    assert checked_bars == 24
    artifacts = {}
    for suffix in ['png', 'svg', 'pdf']:
        path = prefix.with_suffix('.' + suffix); assert not path.exists(); fig.savefig(path, dpi=180); artifacts[str(path)] = sha(path)
    plt.close(fig)
    with dest_table.open() as handle: assert list(csv.DictReader(handle, delimiter='\t')) == rows
    with order_table.open() as handle: assert [{k: str(v) for k, v in r.items()} for r in order_rows] == list(csv.DictReader(handle, delimiter='\t'))
    for p in [dest_table, order_table]: artifacts[str(p)] = sha(p)
    verify(bindings)
    result = dict(status='published_complete_closed_full_triad_order_and_context_counts', full_context_summary_rows=720, full_physical_screen_summary_rows=24,
                  figure_bar_heights_checked=checked_bars, figure_series=series, source_hashes=bindings, artifacts=artifacts, scientific_eligibility=False,
                  scope='Complete closed720 context/screen/policy rows and24 distinct physical all/any-order screen rows. Figure uses only fixed lexical/native-both policy at50residues/70% original coverage; every24barheight/label derives from closed full SQL-checked counts. Full guide/design/core/mask/source-parent denominators retained. Repeated contexts and physical triples are not independent events; quality coverage does not establish duplication effects or selection.')
    with Path('metadata/full_triad_geometry_context_counts_published_20261001.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__': main()
