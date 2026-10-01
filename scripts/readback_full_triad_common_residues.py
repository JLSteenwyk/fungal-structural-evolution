#!/usr/bin/env python3
"""Independently reconstruct every full triad map from native strings and original positions."""
import argparse
import gzip
import itertools
import json
from collections import Counter
from functools import lru_cache
from itertools import zip_longest
from pathlib import Path
from full_triad_common_sources import load_design_inputs, load_measured_edges
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

SUMMARY_FIELDS = ['target_contexts', 'reference_tie_records', 'duplicate_reference_links', 'unique_ordered_model_triads',
                  'correspondence_work_triads', 'mask_order_states', 'needed_models', 'all_input_dispositions', 'input_source_counts',
                  'common_reference_residue_occurrences', 'cycle_consistent_residue_occurrences', 'counts', 'three_pair_pass_state_counts']


def reconstruct(source, cov, inputs, pair, mask, order, bindings):
    path = Path(source['source_checkpoint']); bind(bindings, path, source['source_checkpoint_sha256']); assert sha(path) == source['source_checkpoint_sha256']
    raw = json.loads(path.read_text()); expected_status = source['source_native_status']
    assert (raw['pair_key'], raw['mask'], raw['order'], raw['status']) == (pair, mask, source['source_order'], expected_status)
    assert raw['plan_sha256'] == source['source_plan_sha256'] and raw['input_manifest_sha256'] == source['source_input_manifest_sha256']
    actual = [(i['model_id'], int(i['version'])) for i in raw['inputs']]
    first = cov['model_a'], int(cov['version_a']); second = cov['model_b'], int(cov['version_b'])
    assert actual == ([first, second] if order == 0 else [second, first])
    rows = [inputs[(*model, mask)] for model in actual]; reasons = list(source['numerical_exclusion_reasons'])
    assert source['numerical_usable'] == (expected_status == 'aligned' and reasons == [])
    assert cov[f'order{order}_native_status'] == expected_status
    assert cov[f'order{order}_numerical_exclusion_reasons'] == ';'.join(reasons)
    for old, current in zip(raw['inputs'], rows):
        assert old['status'] == current['status']
        if current['status'] == 'ready': assert old['sha256'] == current['sha256']
    pairs = []
    if expected_status == 'aligned':
        assert all(r['status'] == 'ready' for r in rows)
        left, right = raw['metrics']['alignment_left'], raw['metrics']['alignment_right']
        assert len(left) == len(right) and left.replace('-', '') == rows[0]['sequence'] and right.replace('-', '') == rows[1]['sequence']
        # Independent iterator recovery rather than indexed residue counters.
        positions = [iter(r['original_positions']) for r in rows]
        for aa, bb in zip(left, right):
            a = next(positions[0]) if aa != '-' else None
            b = next(positions[1]) if bb != '-' else None
            if a is not None and b is not None: pairs.append((a, b))
        assert all(next(p, None) is None for p in positions)
        assert len(pairs) == raw['metrics']['aligned_length'] and len(set(pairs)) == len(pairs)
        assert len({p[0] for p in pairs}) == len({p[1] for p in pairs}) == len(pairs)
    else:
        assert expected_status in ['input_unavailable', 'native_error', 'parse_error', 'timeout']; reasons.append(expected_status)
        if expected_status == 'input_unavailable': assert any(r['status'] != 'ready' for r in rows)
    provenance = dict(source_checkpoint=source['source_checkpoint'], source_checkpoint_sha256=source['source_checkpoint_sha256'], source_order=source['source_order'],
                      source_plan_sha256=source['source_plan_sha256'], source_input_manifest_sha256=source['source_input_manifest_sha256'], source_native_status=expected_status,
                      numerical_usable=source['numerical_usable'], numerical_exclusion_reasons=source['numerical_exclusion_reasons'], selected_source=source['selected_source'],
                      current_order=order, directed_endpoints=[list(e) for e in actual], order_status=cov[f'order{order}_status'], mapping_exclusions=reasons)
    return actual, pairs, provenance


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); assert not args.output.exists(); plan = json.loads(args.plan.read_text())
    triads, inputs, design, source_counts, bindings = load_design_inputs(plan); bind(bindings, args.plan)
    checkpoints, coverage = load_measured_edges(plan, triads, bindings)
    out = Path(plan['output']); rp = out / 'receipt.json'; r = json.loads(rp.read_text())
    assert r['status'] == 'complete_full_triad_original_residue_mapping_pending_independent_readback' and r['plan_sha256'] == sha(args.plan)
    assert set(r['artifacts']) == {'common_residue_maps.jsonl.gz'} and r['full_source_design'] == plan['triad_completion'] and r['full_source_design_sha256'] == sha(plan['triad_completion'])
    original_bindings = dict(bindings); bind(bindings, rp); bind(bindings, out / 'common_residue_maps.jsonl.gz', r['artifacts']['common_residue_maps.jsonl.gz']); verify(bindings)
    @lru_cache(maxsize=2048)
    def native(kind, pair, mask, order):
        return reconstruct(checkpoints[(kind, pair, mask, order)], coverage[(kind, pair, mask)], inputs, pair, mask, order, original_bindings)
    expected_grid = ((triad, mask, orders) for triad in triads for mask in ['full', 'plddt70'] for orders in itertools.product([0, 1], repeat=3))
    counts, screen_counts = Counter(), Counter(); states = common_count = cycle_count = 0
    with gzip.open(out / 'common_residue_maps.jsonl.gz', 'rt') as handle:
        for expected, line in zip_longest(expected_grid, handle):
            assert expected is not None and line is not None
            triad, mask, orders = expected; actual = json.loads(line); models = [tuple(m) for m in triad['models']]
            directions = [(models[0], models[1]), (models[2], models[0]), (models[2], models[1])]
            maps, traces, passes, exclusions = [], [], [], []
            for label, kind, desired, order in zip(['ab', 'ar', 'br'], ['primary', 'reference', 'reference'], directions, orders):
                pair = triad['edges'][label]['pair_key']; endpoints, raw_pairs, trace = native(kind, pair, mask, order)
                assert set(endpoints) == set(desired)
                directed = raw_pairs if endpoints == list(desired) else [(b, a) for a, b in raw_pairs]
                maps.append(dict(directed)); traces.append(trace); cov = coverage[(kind, pair, mask)]
                flags, why = {}, {}
                for screen in plan['screens']:
                    sid = screen['id']; flags[sid] = cov[sid + '_pass'] == '1'; assert cov[sid + '_pass'] in ['0', '1']
                    why[sid] = cov[sid + '_exclusions'].split(';') if cov[sid + '_exclusions'] else []
                    assert flags[sid] == (why[sid] == [])
                passes.append(flags); exclusions.append(why)
            # Rebuild reference intersection and AB consistency with an
            # explicit relation, not the producer's dictionary helper.
            ab, ar, br = maps
            common = sorted([[a, br[ref], ref] for ref, a in ar.items() if ref in br], key=lambda t: t[2])
            ab_relation = set(ab.items()); cycle = [t for t in common if (t[0], t[1]) in ab_relation]
            reasons = sorted(set(itertools.chain.from_iterable(t['mapping_exclusions'] for t in traces)))
            if len(common) == 0: reasons += ['no_common_reference_residues']
            elif len(cycle) == 0: reasons += ['no_cycle_consistent_residues']
            all_pass = {screen['id']: sum(int(f[screen['id']]) for f in passes) == 3 for screen in plan['screens']}
            blocked = sum(len(t['mapping_exclusions']) for t in traces) > 0
            core_status = {}
            for name, triples in [('reference_common', common), ('cycle_consistent', cycle)]:
                if blocked: core_status[name] = 'source_excluded'
                elif len(triples) <= 2: core_status[name] = 'fewer_than_three_common_residues'
                else: core_status[name] = 'pending_common_coordinate_geometry'
            expected_row = dict(triad_id=triad['triad_id'], role_order=['a', 'b', 'reference'], models=triad['models'], mask=mask, orders=list(orders),
                                edge_order=['ab', 'reference_to_a', 'reference_to_b'], edge_provenance=traces,
                                original_lengths=[inputs[(*model, 'full')]['original_length'] for model in models],
                                retained_input_lengths=[inputs[(*model, mask)]['retained_residues'] for model in models],
                                edge_pair_screen_pass=passes, edge_pair_screen_exclusions=exclusions, all_three_pair_screen_pass=all_pass,
                                reference_common_triples=common, cycle_consistent_triples=cycle, common_reference_count=len(common), cycle_consistent_count=len(cycle),
                                common_core_fit_input_status=core_status, mapping_exclusions=reasons)
            assert actual == expected_row
            common_count += len(common); cycle_count += len(cycle); states += 1
            counts[mask + ':' + (';'.join(reasons) if reasons else 'recorded_without_edge_exclusions')] += 1
            for sid, passed in all_pass.items(): screen_counts[mask + ':' + sid] += int(passed)
            if states % 8000 == 0: print('Independent full native/common-residue states', states, '/', plan['expected']['potential_correspondence_states'], flush=True)
    assert states == plan['expected']['potential_correspondence_states'] == 16 * len(triads)
    assert r['source_hashes'] == original_bindings
    for path, h in original_bindings.items(): bind(bindings, path, h)
    summary = dict(target_contexts=design['target_contexts'], reference_tie_records=design['reference_tie_records'], duplicate_reference_links=design['duplicate_reference_links'],
                   unique_ordered_model_triads=design['unique_ordered_model_triads'], correspondence_work_triads=len(triads), mask_order_states=states,
                   needed_models=len(inputs) // 2, all_input_dispositions=plan['expected']['all_input_dispositions'], input_source_counts=source_counts,
                   common_reference_residue_occurrences=common_count, cycle_consistent_residue_occurrences=cycle_count, counts=dict(counts), three_pair_pass_state_counts=dict(screen_counts))
    assert set(summary) == set(SUMMARY_FIELDS) and all(r[k] == v for k, v in summary.items()); verify(bindings)
    result = dict(status='passed_full_triad_original_residue_mapping_readback', plan_sha256=sha(args.plan), producer_receipt_sha256=sha(rp),
                  checker_sha256=sha(__file__), **summary, source_hashes=bindings, scientific_eligibility=False,
                  scope='Every full scheduled triad/mask/eight-order state independently rebuilt from raw checkpoint strings and original residue-position iterators, actual current/original source directions, checksums, sequence and written-input identity. Reference intersections and AB relation cycles independently checked along with all status/numerical/missing/coverage flags, full grids and aggregates. Entire original design remains immutably linked, including unscheduled/missing/excluded contexts. Shares source/proof I/O only; no producer mapping or intersection code. No common-coordinate fit, accepted phylogeny/biological orthology/prediction accuracy or calibrated duplication effect.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
