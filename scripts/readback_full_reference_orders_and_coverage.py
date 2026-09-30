#!/usr/bin/env python3
"""Independently reconstruct every full-reference order, metric and exact coverage decision."""
import argparse
import csv
import json
import math
from collections import Counter
from decimal import Decimal, ROUND_CEILING
from pathlib import Path

from reference_order_coverage_sources import load_sources, groups
from reference_measurement_union_sources import bind, verify, pair_ends
from run_ortholog_pair_guide_comparison import sha


def reconstruct(pair, mask, records, ends, models, inputs, screens):
    assert len(records) == 2 and sorted(r['order'] for r, _ in records) == [0, 1]
    row = dict(pair_key=pair, mask=mask, model_a=ends[0][0], version_a=ends[0][1],
               model_b=ends[1][0], version_b=ends[1][1], length_a=models[ends[0]]['length'], length_b=models[ends[1]]['length'])
    values = {}; usable = {}
    metric_names = ['aligned_length', 'rmsd_recomputed', 'sequence_identity_exact', 'tm_a', 'tm_b',
                    'retained_input_coverage_a', 'retained_input_coverage_b', 'joint_plddt70_fraction',
                    'original_coverage_a', 'original_coverage_b']
    for r, digest in records:
        order = r['order']; prefix = f'order{order}_'
        assert r['pair_key'] == pair and r['mask'] == mask and pair_ends(r) == ends
        directed = [(x['model_id'], x['version']) for x in r['directed_endpoints']]
        assert directed == (ends if order == 0 else list(reversed(ends))) and r['source_order'] in [0, 1]
        status = r['source_native_status']; n = r['source_numeric']; g = r['source_geometry']
        assert status in ['aligned', 'input_unavailable', 'native_error', 'parse_error', 'timeout']
        expected_reasons = []
        if status == 'aligned':
            assert n is not None and g is not None and r['native_metrics'] is not None
            for original in [n, g]:
                assert (original['pair_key'], original['mask'], int(original['order'])) == (pair, mask, r['source_order'])
            assert n['aligned_length'] == g['aligned_length'] and n['rmsd_status'] == g['rmsd_status']
            count = Decimal(n['aligned_length']); assert count == count.to_integral_value()
            assert int(count) == r['native_metrics']['aligned_length']
            assert 0 < count <= min(inputs[(*end, mask)]['retained_residues'] for end in ends)
            assert all(inputs[(*end, mask)]['status'] == 'ready' for end in ends)
            if n['rmsd_status'] != 'within_printed_rounding': expected_reasons.append('rmsd_discrepancy')
            if count <= 2: expected_reasons.append('fewer_than_three_pairs')
            if g['geometry_status'] != 'unique_at_numeric_tolerance': expected_reasons.append('nonunique_rotation')
        else:
            assert n is None and g is None and r['native_metrics'] is None
            expected_reasons = [status]
        assert expected_reasons == r['numerical_exclusion_reasons'] and r['numerical_usable'] == (len(expected_reasons) == 0)
        usable[order] = len(expected_reasons) == 0
        row[prefix + 'status'] = 'aligned' if usable[order] else 'excluded_numerically' if status == 'aligned' else status
        row.update({prefix + field: value for field, value in {
            'source': r['selected_source'], 'source_order': r['source_order'], 'source_checkpoint': r['source_checkpoint'],
            'source_checkpoint_sha256': r['source_checkpoint_sha256'], 'union_row_sha256': digest,
            'native_status': status, 'numerical_exclusion_reasons': ';'.join(expected_reasons),
            'geometry_status': g['geometry_status'] if g else '', 'rmsd_status': n['rmsd_status'] if n else '',
            'native_aligned_length': n['aligned_length'] if n else ''}.items()})
        for side, end in zip(['a', 'b'], ends):
            row[prefix + 'retained_residues_' + side] = inputs[(*end, mask)]['retained_residues']
            row[prefix + 'input_status_' + side] = inputs[(*end, mask)]['status']
        if usable[order]:
            values[order] = {field: float(n[field]) for field in ['aligned_length', 'rmsd_recomputed', 'sequence_identity_exact', 'joint_plddt70_fraction']}
            for side, end in zip(['a', 'b'], ends):
                native_side = ['left', 'right'][directed.index(end)]
                values[order]['tm_' + side] = float(n['tm_' + native_side + '_native'])
                values[order]['retained_input_coverage_' + side] = float(n['coverage_' + native_side])
                values[order]['original_coverage_' + side] = float(count / Decimal(models[end]['length']))
            assert all(math.isfinite(v) and v >= 0 for v in values[order].values())
            assert all(values[order][k] <= 1 for k in metric_names if k not in ['aligned_length', 'rmsd_recomputed'])
        for name in metric_names: row[prefix + name] = values[order][name] if usable[order] else ''
    for name in metric_names:
        row[name + '_order_absolute_difference'] = abs(values[0][name] - values[1][name]) if all(usable.values()) else ''
    number = sum(usable.values())
    row['order_summary_status'] = {2: 'both_orders_numerically_usable', 1: 'one_order_numerically_usable', 0: 'neither_order_numerically_usable'}[number]
    for spec in screens:
        why = [f'order{o}_not_numerically_usable' for o in [0, 1] if not usable[o]]
        aligned = [int(Decimal(records[ix][0]['source_numeric']['aligned_length'])) for ix in range(2) if records[ix][0]['numerical_usable']]
        if any(n < spec['minimum_aligned_residues'] for n in aligned): why.append('short_alignment')
        required = max(int((Decimal(models[end]['length']) * Decimal(str(spec['minimum_original_coverage']))).to_integral_value(rounding=ROUND_CEILING)) for end in ends)
        if any(n < required for n in aligned): why.append('low_original_protein_coverage')
        row[spec['id'] + '_pass'] = int(not why); row[spec['id'] + '_exclusions'] = ';'.join(why)
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); assert not args.output.exists(); plan = json.loads(args.plan.read_text())
    source, pairs, models, inputs, upstream, bindings = load_sources(plan); bind(bindings, args.plan)
    root = Path(plan['output']); receipt_path = root / 'receipt.json'; r = json.loads(receipt_path.read_text())
    assert r['status'] == 'complete_full_reference_order_and_original_coverage_pending_independent_readback'
    assert r['plan_sha256'] == sha(args.plan) and r['screens'] == plan['screens']
    assert r['source_hashes'] == bindings
    artifacts_expected = {'pair_mask_order_coverage.tsv'}; assert set(r['artifacts']) == artifacts_expected
    for name, digest in r['artifacts'].items(): bind(bindings, root / name, digest)
    bind(bindings, receipt_path); verify(bindings)
    seen, decisions, counts, passes, excluded = set(), {}, Counter(), Counter(), Counter()
    sources, statuses, numeric_counts = Counter(), Counter(), Counter(); maxima = None
    with (root / 'pair_mask_order_coverage.tsv').open() as handle:
        actual = csv.DictReader(handle, delimiter='\t')
        for key, records in groups(source):
            assert key not in seen and key[0] in pairs and key[1] in ['full', 'plddt70']; seen.add(key)
            expected = reconstruct(*key, records, pair_ends(pairs[key[0]]), models, inputs, plan['screens'])
            row = next(actual); assert set(actual.fieldnames) == set(expected) and len(actual.fieldnames) == len(expected)
            for field, value in expected.items():
                if isinstance(value, float):
                    assert math.isfinite(float(row[field])) and math.isclose(float(row[field]), value, rel_tol=1e-15, abs_tol=1e-15), (key, field)
                else: assert row[field] == str(value), (key, field)
            if maxima is None: maxima = {k: None for k in expected if k.endswith('_order_absolute_difference')}
            for field in maxima:
                if expected[field] != '': maxima[field] = expected[field] if maxima[field] is None else max(maxima[field], expected[field])
            counts[key[1] + ':' + expected['order_summary_status']] += 1; decisions[key] = {}
            for spec in plan['screens']:
                name = spec['id']; decisions[key][name] = bool(expected[name + '_pass']); passes[key[1] + ':' + name] += expected[name + '_pass']
                for reason in expected[name + '_exclusions'].split(';'):
                    if reason: excluded[key[1] + ':' + name + ':' + reason] += 1
            for state, _ in records:
                sources[state['selected_source']] += 1; statuses[key[1] + ':' + state['source_native_status']] += 1
                numeric_counts['usable' if state['numerical_usable'] else ';'.join(state['numerical_exclusion_reasons'])] += 1
            if len(seen) % 10000 == 0: print('Independent full reference pair/mask coverage', len(seen), '/', 2 * len(pairs), flush=True)
        assert next(actual, None) is None
    assert seen == {(key, mask) for key in pairs for mask in ['full', 'plddt70']}
    joint = {s['id']: sum(decisions[pair, 'full'][s['id']] and decisions[pair, 'plddt70'][s['id']] for pair in pairs) for s in plan['screens']}
    summary = dict(full_pairs=len(pairs), directed_dispositions=2 * len(seen), pair_mask_rows=len(seen), pair_screen_decisions=len(seen) * len(plan['screens']),
                   order_summary_counts=dict(counts), maximum_order_differences=maxima, pair_pass_counts=dict(passes), pair_exclusion_counts=dict(excluded),
                   pair_both_masks_pass_counts=joint, source_dispositions=dict(sources), native_status_counts=dict(statuses), numerical_counts=dict(numeric_counts))
    assert all(r[k] == v for k, v in summary.items())
    for field in ['source_dispositions', 'native_status_counts', 'numerical_counts']: assert summary[field] == upstream[field]
    verify(bindings)
    result = dict(status='passed_full_reference_order_and_original_coverage_readback', plan_sha256=sha(args.plan),
                  producer_receipt_sha256=sha(receipt_path), checker_sha256=sha(__file__), **summary, source_hashes=bindings,
                  scientific_eligibility=False, scope='Every pair/mask/order source identity, native status/exclusion, numerical metric, endpoint normalization, retained and full-protein coverage, blank unavailable field and order difference independently reconstructed from the closed full union. All six pass/exclusion decisions replayed with decimal ceilings against original lengths, all aggregates and full/both-mask intersections checked. Shared helper is source/proof I/O only; producer metric normalization and filtering functions are not used. No biological reference/asymmetry or calibrated effect acceptance.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
