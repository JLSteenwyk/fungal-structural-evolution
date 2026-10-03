#!/usr/bin/env python3
"""Publish the closed whole-attempt recovery overlay for all original quartets/cutoffs."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
import numpy as np
from background_measurement_union_sources import closed_source
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

CUTS = ['250', '500']
PRIORS = ['package', 'centered', 'broad']
SCALAR = ['passes_scalar_screen_only', 'scalar_mixing_requires_review', 'constant_chain_requires_review']
CATEGORY = ['observed_state_indicator_screens_pass_only', 'categorical_mixing_requires_review', 'no_observed_state_variation_requires_review']
FIELDS = ['model_input_identity', 'effective_input_group', 'prior', 'family', 'proteins',
    'original_configuration_ids', 'chain_ids', 'seeds', 'checked_chains', 'failed_chains',
    'allocation_notice_chains', 'cutoff', 'disposition', 'scalar_draws_per_chain',
    'saved_state_draws_per_chain', 'scalar_variables', 'scalar_pass_only',
    'scalar_mixing_review', 'scalar_constant_review', 'every_scalar_screen_pass',
    'length_variables', 'length_pass_only', 'length_mixing_review', 'length_constant_review',
    'every_length_screen_pass', 'categorical_coordinates', 'categorical_patterns',
    'coordinate_indicator_pass_only', 'coordinate_mixing_review', 'coordinate_no_variation_review',
    'pattern_indicator_pass_only', 'pattern_mixing_review', 'pattern_no_variation_review',
    'scientific_eligibility']


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); a = p.parse_args()
    plan = json.loads(a.plan.read_text()); bindings = dict(plan['pins']); bind(bindings, a.plan)
    assert not any(Path(plan['figure_stem'] + '.' + ext).exists() for ext in ['png', 'svg', 'pdf'])
    assert not Path(plan['table']).exists() and not Path(plan['completion']).exists()
    pc = closed_source(plan['diagnostic_completion'], 'complete_verified_full_baliphy_recovery_diagnostics',
        'complete_verified_full_baliphy_recovery_diagnostics_archive', 2, bindings)
    hc = pc
    original_bindings = {str(Path(q).resolve()): d for q, d in bindings.items()}
    def doc(path, expected_hash=None):
        path = Path(path); digest = original_bindings[str(path.resolve())]
        if expected_hash is not None: assert digest == expected_hash
        bind(bindings, path, digest)
        return json.loads(path.read_text())
    producer = doc(plan['producer_plan']); jobs = doc(producer['jobs'])
    assert len(jobs) == 1620 and len({j['chain']['chain_id'] for j in jobs}) == 1620
    groups = defaultdict(list)
    for j in jobs: groups[j['config']['model_input_identity']].append(j)
    assert len(groups) == 405 and all(len(v) == 4 for v in groups.values())
    recovery = doc(plan['recovery_completion'])
    overlay = doc(recovery['overlay'], recovery['overlay_sha256'])
    overlay_rows = {r['chain']['chain_id']: r for r in overlay['rows']}
    assert len(overlay_rows) == len(overlay['rows']) == 1620
    assert all(overlay_rows[j['chain']['chain_id']]['chain'] == j['chain'] for j in jobs)
    assert overlay['original_samples_concatenated'] is False
    assert overlay['scientific_eligibility'] is False
    ns = {cid: row['selected_disposition']['status'] for cid, row in overlay_rows.items()}
    original_horizon = doc(plan['horizon_completion'])
    haccount = doc(Path(original_horizon['full_hash_archive']).parent / 'accounting.json')
    original_notices = {x['chain_id'] for x in haccount['allocation_notices']}
    assert len(original_notices) == 5
    recovery_plan = doc(recovery['source_plan'], recovery['source_plan_sha256'])
    audit = doc(Path(recovery_plan['output']) / 'receipt.json')
    retried = {x['chain_id'] for x in audit['attempts']}
    assert len(retried) == len(audit['attempts']) == 3
    notices = (original_notices - retried) | {x['chain_id'] for x in audit['attempts'] if x['allocation_warnings'] > 0}
    diagnostics = doc(pc['producer_receipt'], pc['producer_receipt_sha256'])
    assert set(diagnostics['groups']) == set(groups)
    controllers = {name: {'groups': {group: dict(info.get(name, {}), status=info['status'])
        for group, info in diagnostics['groups'].items()}} for name in ['scalar', 'length', 'categorical']}
    rows = []; aggregation = Counter(); complete_by_prior = Counter(); missing_by_prior = Counter()
    for group, js in sorted(groups.items()):
        ch = js[0]['chain']; ids = sorted(j['chain']['chain_id'] for j in js)
        assert len({j['chain']['prior_label'] for j in js}) == 1 and ch['prior_label'] in PRIORS
        assert len({j['chain']['seed'] for j in js}) == 4
        for field in ['effective_input_group', 'family', 'proteins', 'original_configuration_ids']:
            assert all(j['chain'][field] == ch[field] for j in js)
        checked = sum(ns[c] == 'all_saved_alignments_and_candidate_nodes_checked' for c in ids)
        assert all(ns[c] in ['all_saved_alignments_and_candidate_nodes_checked', 'failed'] for c in ids)
        base = dict(model_input_identity=group, effective_input_group=ch['effective_input_group'], prior=ch['prior_label'],
            family=ch['family'], proteins=ch['proteins'], original_configuration_ids=json.dumps(ch['original_configuration_ids'], separators=(',', ':')),
            chain_ids=json.dumps(ids, separators=(',', ':')), seeds=json.dumps(sorted(j['chain']['seed'] for j in js)),
            checked_chains=checked, failed_chains=4-checked, allocation_notice_chains=len(notices & set(ids)), scientific_eligibility=0)
        states = {name: c['groups'][group] for name, c in controllers.items()}
        if checked != 4:
            assert states['scalar']['status'] == 'unresolved_failed_native_chain_retained'
            assert states['length']['status'] == 'unresolved_failed_native_chain_retained'
            assert states['categorical']['status'] == 'unresolved_failed_native_chain_retained'
            missing_by_prior[ch['prior_label']] += 1
            for cutoff in CUTS: rows.append({**base, 'cutoff': cutoff, 'disposition': 'unresolved_native_memory_failure'})
            continue
        complete_by_prior[ch['prior_label']] += 1
        assert states['scalar']['status'] == 'complete_scalar_length_category_screens_not_posterior_qualification'
        assert states['length']['status'] == 'complete_scalar_length_category_screens_not_posterior_qualification'
        assert states['categorical']['status'] == 'complete_scalar_length_category_screens_not_posterior_qualification'
        receipts = {name: doc(x['receipt'], x['receipt_sha256']) for name, x in states.items()}
        categorical = receipts['categorical']; report = doc(categorical['report'], categorical['report_sha256'])
        for cutoff in CUTS:
            s = receipts['scalar']['outputs'][cutoff]; sd = doc(s['path'], s['sha256'])
            l = receipts['length']['outputs'][cutoff]; ld = doc(l['path'], l['sha256'])
            scalar = Counter(v['status'] for v in sd['variables'].values())
            length = Counter(v['status'] for v in ld['variables'].values())
            assert set(scalar) <= set(SCALAR) and set(length) <= set(SCALAR)
            assert len(sd['variables']) == 37 and len(ld['variables']) == 4
            assert all(v['draws_per_chain'] == 1000-int(cutoff) for v in sd['variables'].values())
            assert all(v['draws_per_chain'] == (75 if cutoff == '250' else 50) for v in ld['variables'].values())
            cat = doc(Path(categorical['report']).parent / ('discard-' + cutoff) / 'summary.json')
            coordinates = Counter(cat['coordinate_status_counts']); patterns = Counter(cat['pattern_status_counts'])
            assert set(coordinates) <= set(CATEGORY) and set(patterns) <= set(CATEGORY)
            assert sum(coordinates.values()) == cat['coordinates'] == report['outputs'][cutoff]['coordinates']
            assert sum(patterns.values()) == cat['patterns'] == report['outputs'][cutoff]['patterns']
            assert cat['retained_samples_per_chain'] == (75 if cutoff == '250' else 50)
            row = dict(base, cutoff=cutoff, disposition='diagnostics_complete_requires_mixing_review',
                scalar_draws_per_chain=1000-int(cutoff), saved_state_draws_per_chain=cat['retained_samples_per_chain'],
                scalar_variables=sum(scalar.values()), scalar_pass_only=scalar[SCALAR[0]], scalar_mixing_review=scalar[SCALAR[1]],
                scalar_constant_review=scalar[SCALAR[2]], every_scalar_screen_pass=int(scalar[SCALAR[0]] == 37),
                length_variables=sum(length.values()), length_pass_only=length[SCALAR[0]], length_mixing_review=length[SCALAR[1]],
                length_constant_review=length[SCALAR[2]], every_length_screen_pass=int(length[SCALAR[0]] == 4),
                categorical_coordinates=cat['coordinates'], categorical_patterns=cat['patterns'],
                coordinate_indicator_pass_only=coordinates[CATEGORY[0]], coordinate_mixing_review=coordinates[CATEGORY[1]],
                coordinate_no_variation_review=coordinates[CATEGORY[2]], pattern_indicator_pass_only=patterns[CATEGORY[0]],
                pattern_mixing_review=patterns[CATEGORY[1]], pattern_no_variation_review=patterns[CATEGORY[2]])
            rows.append(row)
            for diagnostic, counts in [('scalar', scalar), ('length', length), ('coordinate', coordinates), ('pattern', patterns)]:
                aggregation.update({(cutoff, ch['prior_label'], diagnostic, status): n for status, n in counts.items()})
    assert len(rows) == 810 and sum(complete_by_prior.values()) == 403 and sum(missing_by_prior.values()) == 2
    assert all(complete_by_prior[x] + missing_by_prior[x] == 135 for x in PRIORS)
    for cutoff in CUTS:
        for state in SCALAR:
            assert sum(aggregation[cutoff, p, 'scalar', state] for p in PRIORS) == hc['scalar_status_counts'].get(cutoff+':'+state, 0)
            assert sum(aggregation[cutoff, p, 'length', state] for p in PRIORS) == pc['length_status_counts'].get(cutoff+':'+state, 0)
        for diagnostic, field in [('coordinate', 'coordinate_status_counts'), ('pattern', 'pattern_status_counts')]:
            for state in CATEGORY: assert sum(aggregation[cutoff, p, diagnostic, state] for p in PRIORS) == pc[field].get(cutoff+':'+state, 0)
        selected = [r for r in rows if r['cutoff'] == cutoff and r['disposition'].startswith('diagnostics_complete')]
        assert sum(r['every_scalar_screen_pass'] for r in selected) == hc['quartets_passing_every_scalar'][cutoff] == 0
        assert sum(r['every_length_screen_pass'] for r in selected) == pc['quartets_passing_every_length_scalar'][cutoff] == 0
        assert sum(r['categorical_coordinates'] for r in selected) == pc['categorical_coordinate_counts'][cutoff]
        assert sum(r['categorical_patterns'] for r in selected) == pc['categorical_pattern_counts'][cutoff]
    print('All 810 original quartet/cutoff rows reconciled; verifying full closed sources', flush=True)
    verify(bindings)
    table = Path(plan['table']); assert not table.exists(); table.parent.mkdir(exist_ok=True, parents=True)
    with table.open('x') as f:
        writer = csv.DictWriter(f, FIELDS, delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    with table.open() as f:
        saved = list(csv.DictReader(f, delimiter='\t'))
    assert saved == [{field: str(r.get(field, '')) for field in FIELDS} for r in rows]
    panels = []; expected_text = []
    fig, axes = plt.subplots(4, 2, figsize=(13, 13))
    labels = ['Passes marginal screen only', 'Mixing requires review', 'Constant/no observed variation; review']
    for row, (diagnostic, title) in enumerate([('scalar', 'Scalar variables'), ('length', 'Candidate-node lengths'),
                                              ('pattern', 'Distinct temporal patterns'), ('coordinate', 'Node/residue anchors')]):
        states = SCALAR if diagnostic in ['scalar', 'length'] else CATEGORY
        for column, cutoff in enumerate(CUTS):
            ax = axes[row, column]
            counts = np.array([[aggregation[cutoff, prior, diagnostic, state] for prior in PRIORS] for state in states])
            totals = counts.sum(axis=0); assert np.all(totals > 0)
            fractions = counts / totals * 100; left = np.zeros(3)
            for values, color, label in zip(fractions, ['#469c9c', '#d77a45', '#8c87b5'], labels):
                ax.barh(np.arange(3), values, left=left, color=color, label=label)
                for i, (value, offset) in enumerate(zip(values, left)):
                    if value >= 7:
                        text = f'{value:.1f}%'; ax.text(offset+value/2, i, text, ha='center', va='center', color='white', fontsize=9); expected_text.append(text)
                left += values
            for i, total in enumerate(totals):
                text = f'n={int(total):,}'; ax.text(102, i, text, va='center', fontsize=9); expected_text.append(text)
            ax.set_title(title + ' / discard through ' + cutoff, fontsize=11)
            ax.set_yticks(range(3), ['Package', 'Centered', 'Broad']); ax.invert_yaxis()
            ax.set_xlim(0, 125); ax.set_xticks([0,25,50,75,100]); ax.set_xlabel('Percent of reported diagnostic units')
            ax.spines[['top','right']].set_visible(False)
            panels.append(dict(diagnostic=diagnostic, cutoff=cutoff, priors=PRIORS, states=states, counts=counts.tolist(), totals=totals.tolist()))
    handles, legends = axes[0,0].get_legend_handles_labels(); fig.legend(handles, legends, loc='lower center', ncol=1, frameon=False, bbox_to_anchor=(.5,.025))
    fig.suptitle('Recovered ancestral first-horizon diagnostics: 403 completed four-chain groups\nNo group passes every scalar or candidate-length screen', fontsize=16)
    fig.text(.5,.012,'Two native memory-failed groups retained in the full table. Units are dependent; constant states do not establish convergence.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.105,1,.94))
    figures = []
    for extension in ['png','svg','pdf']:
        path = Path(plan['figure_stem'] + '.' + extension); assert not path.exists(); path.parent.mkdir(exist_ok=True, parents=True)
        fig.savefig(path, dpi=160); figures.append(path)
    plt.close(fig)
    svg = figures[1]; tree = ET.parse(svg); ns = '{http://www.w3.org/2000/svg}'
    text = ' '.join(''.join(e.itertext()) for e in tree.iter(ns+'text'))
    assert all(t in text for t in expected_text)
    metadata = ET.Element(ns+'metadata', id='full-ancestral-diagnostic-counts'); metadata.text = json.dumps(panels, separators=(',', ':'))
    tree.getroot().append(metadata); tree.write(svg, encoding='utf-8', xml_declaration=True)
    assert json.loads(ET.parse(svg).find(ns+'metadata[@id="full-ancestral-diagnostic-counts"]').text) == panels
    result = dict(status='published_full_recovered_horizon_ancestral_diagnostic_table_and_figure',
        selected_checked_chains=1618, unresolved_failed_chains=2, recovered_whole_chains=1,
        original_samples_concatenated=False, selected_allocation_notice_chains=len(notices),
        original_quartets=405, complete_quartets=403, unresolved_quartets=2, cutoff_rows=810, priors=PRIORS,
        complete_by_prior=dict(complete_by_prior), unresolved_by_prior=dict(missing_by_prior), plotted_panels=8,
        every_source_count_reconciled=True, full_table_serialized_values_checked=True, svg_labels_and_metadata_checked=True,
        source_completions=[plan['diagnostic_completion'],plan['recovery_completion']],
        closed_source_bindings_reverified=len(bindings), original_journals_in_immediate_diagnostic_closure=2,
        source_hashes={str(q):sha(q) for q in [a.plan,Path(__file__),plan['diagnostic_completion'],plan['recovery_completion']]},
        artifacts={str(q):sha(q) for q in [table,*figures]}, scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['completion']).open('x') as f: f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__ == '__main__': main()
