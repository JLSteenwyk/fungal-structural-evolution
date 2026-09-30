#!/usr/bin/env python3
"""Normalize both reference measurement orders and apply full-protein coverage screens."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from reference_order_coverage_sources import load_sources, groups
from reference_measurement_union_sources import bind, verify, pair_ends
from summarize_duplication_alignment_orders import summarize_pair, METRICS
from screen_background_whole_protein_coverage import reasons
from run_ortholog_pair_guide_comparison import sha


def summarize(pair, mask, records, ends, models, inputs):
    states = {r['order']: (r, digest) for r, digest in records}
    assert len(records) == len(states) == 2 and set(states) == {0, 1}
    statuses, numeric, extra = {}, {}, {}
    for order in [0, 1]:
        r, digest = states[order]
        assert (r['pair_key'], r['mask']) == (pair, mask) and pair_ends(r) == ends
        directed = ends if order == 0 else ends[::-1]
        assert r['directed_endpoints'] == [dict(model_id=m, version=v) for m, v in directed]
        assert r['source_order'] in [0, 1]
        n, g = r['source_numeric'], r['source_geometry']
        status = r['source_native_status']; why = r['numerical_exclusion_reasons']
        assert status in ['aligned', 'input_unavailable', 'native_error', 'parse_error', 'timeout']
        assert r['numerical_usable'] == (not why)
        if status == 'aligned':
            assert n is not None and g is not None and r['native_metrics'] is not None
            assert (n['pair_key'], n['mask'], int(n['order'])) == (pair, mask, r['source_order'])
            assert (g['pair_key'], g['mask'], int(g['order'])) == (pair, mask, r['source_order'])
            assert n['aligned_length'] == g['aligned_length'] and n['rmsd_status'] == g['rmsd_status']
            expected = []
            if n['rmsd_status'] != 'within_printed_rounding': expected.append('rmsd_discrepancy')
            if int(n['aligned_length']) < 3: expected.append('fewer_than_three_pairs')
            if g['geometry_status'] != 'unique_at_numeric_tolerance': expected.append('nonunique_rotation')
            assert why == expected
            assert int(n['aligned_length']) == r['native_metrics']['aligned_length']
            assert 0 < int(n['aligned_length']) <= min(inputs[(*end, mask)]['retained_residues'] for end in ends)
            assert all(inputs[(*end, mask)]['status'] == 'ready' for end in ends)
        else:
            assert n is None and g is None and r['native_metrics'] is None and why == [status]
        statuses[order] = ('aligned' if not why else 'excluded_numerically') if status == 'aligned' else status
        if not why: numeric[order] = n
        prefix = f'order{order}_'
        for name, value in [('source', r['selected_source']), ('source_order', r['source_order']),
                            ('source_checkpoint', r['source_checkpoint']), ('source_checkpoint_sha256', r['source_checkpoint_sha256']),
                            ('union_row_sha256', digest), ('native_status', status),
                            ('numerical_exclusion_reasons', ';'.join(why)),
                            ('geometry_status', g['geometry_status'] if g else ''),
                            ('rmsd_status', n['rmsd_status'] if n else ''),
                            ('native_aligned_length', n['aligned_length'] if n else '')]: extra[prefix + name] = value
        for side, end in zip(['a', 'b'], ends):
            extra[prefix + 'retained_residues_' + side] = inputs[(*end, mask)]['retained_residues']
            extra[prefix + 'input_status_' + side] = inputs[(*end, mask)]['status']
            extra[prefix + 'original_coverage_' + side] = int(n['aligned_length']) / models[end]['length'] if not why else ''
    row = summarize_pair(pair, mask, statuses, numeric)
    row = {k.replace('coverage_', 'retained_input_coverage_'): v for k, v in row.items()}
    row['order_summary_status'] = row['order_summary_status'].replace('aligned', 'numerically_usable')
    row.update(model_a=ends[0][0], version_a=ends[0][1], model_b=ends[1][0], version_b=ends[1][1],
               length_a=models[ends[0]]['length'], length_b=models[ends[1]]['length'], **extra)
    for side in ['a', 'b']:
        row['original_coverage_' + side + '_order_absolute_difference'] = (
            abs(row['order0_original_coverage_' + side] - row['order1_original_coverage_' + side])
            if len(numeric) == 2 else '')
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text())
    source, pairs, models, inputs, upstream, bindings = load_sources(plan); bind(bindings, args.plan)
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    seen, decisions, order_counts, passes, excluded = set(), {}, Counter(), Counter(), Counter()
    source_counts, native_counts, numeric_counts = Counter(), Counter(), Counter()
    maxima = {k.replace('coverage_', 'retained_input_coverage_') + '_order_absolute_difference': None for k in METRICS}
    maxima.update({f'original_coverage_{s}_order_absolute_difference': None for s in ['a', 'b']})
    with (out / 'pair_mask_order_coverage.tsv').open('x') as handle:
        writer = None
        for key, records in groups(source):
            pair, mask = key
            assert key not in seen and pair in pairs and mask in ['full', 'plddt70']
            seen.add(key)
            row = summarize(pair, mask, records, pair_ends(pairs[pair]), models, inputs)
            work = dict(row, measurement_disposition='new_model_pair', **{'fit_' + k: v for k, v in row.items()})
            decisions[key] = {}
            for spec in plan['screens']:
                name = spec['id']; why = reasons(work, spec)
                row[name + '_pass'] = int(not why); row[name + '_exclusions'] = ';'.join(why)
                decisions[key][name] = not why; passes[mask + ':' + name] += not why
                for reason in why: excluded[mask + ':' + name + ':' + reason] += 1
            for r, _ in records:
                source_counts[r['selected_source']] += 1; native_counts[mask + ':' + r['source_native_status']] += 1
                numeric_counts['usable' if r['numerical_usable'] else ';'.join(r['numerical_exclusion_reasons'])] += 1
            order_counts[mask + ':' + row['order_summary_status']] += 1
            for field in maxima:
                if row[field] != '': maxima[field] = row[field] if maxima[field] is None else max(maxima[field], row[field])
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t', lineterminator='\n'); writer.writeheader()
            writer.writerow(row)
            if len(seen) % 10000 == 0: print('Full reference pair/mask coverage', len(seen), '/', 2 * len(pairs), flush=True)
    assert seen == {(pair, mask) for pair in pairs for mask in ['full', 'plddt70']}
    for field, actual in [('source_dispositions', source_counts), ('native_status_counts', native_counts), ('numerical_counts', numeric_counts)]:
        assert dict(actual) == upstream[field]
    joint = {s['id']: sum(decisions[key, 'full'][s['id']] and decisions[key, 'plddt70'][s['id']] for key in pairs) for s in plan['screens']}
    verify(bindings)
    result = dict(status='complete_full_reference_order_and_original_coverage_pending_independent_readback',
                  plan_sha256=sha(args.plan), full_pairs=len(pairs), directed_dispositions=2 * len(seen), pair_mask_rows=len(seen),
                  pair_screen_decisions=len(seen) * len(plan['screens']), screens=plan['screens'],
                  order_summary_counts=dict(order_counts), maximum_order_differences=maxima,
                  pair_pass_counts=dict(passes), pair_exclusion_counts=dict(excluded), pair_both_masks_pass_counts=joint,
                  source_dispositions=dict(source_counts), native_status_counts=dict(native_counts), numerical_counts=dict(numeric_counts),
                  artifacts={'pair_mask_order_coverage.tsv': sha(out / 'pair_mask_order_coverage.tsv')},
                  source_hashes=bindings, scientific_eligibility=False,
                  scope='All full-reference pair/mask/order states normalized to current ledger endpoints, while original source orders, checkpoint and row hashes, native flags and excluded raw aligned lengths remain explicit. Both orders required for each unchanged 30/50-residue and 50/70/90-percent screen. Original full protein lengths are denominators under both masks; retained-input coverage separately labeled. No favorable order, source reselection, metric averaging, missing-as-zero, native rerun, biological reference acceptance or asymmetry inference. Full-context/event linkage, common triads, sequence-locked/domain/PAE, prediction uncertainty, phylogeny and calibration remain downstream.')
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}, indent=2), flush=True)


if __name__ == '__main__': main()
