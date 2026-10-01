#!/usr/bin/env python3
"""Map the full source-ready triad grid using original residue positions and both orders."""
import argparse
import gzip
import hashlib
import itertools
import json
import shutil
from collections import Counter
from functools import lru_cache
from pathlib import Path
from full_triad_common_sources import load_design_inputs, load_measured_edges
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def recover_native(source, coverage, inputs, pair, mask, order, bindings):
    path = source['source_checkpoint']; bind(bindings, path, source['source_checkpoint_sha256']); assert sha(path) == source['source_checkpoint_sha256']
    raw = json.loads(Path(path).read_text()); status = source['source_native_status']
    assert raw['pair_key'] == pair and raw['mask'] == mask and raw['order'] == source['source_order'] and raw['status'] == status
    assert raw['plan_sha256'] == source['source_plan_sha256'] and raw['input_manifest_sha256'] == source['source_input_manifest_sha256']
    actual = [(i['model_id'], int(i['version'])) for i in raw['inputs']]
    ends = [(coverage['model_' + s], int(coverage['version_' + s])) for s in ['a', 'b']]
    assert actual == (ends if order == 0 else ends[::-1]) and len(set(actual)) == 2
    rows = [inputs[(*model, mask)] for model in actual]
    for old, current in zip(raw['inputs'], rows):
        assert old['status'] == current['status']
        if current['status'] == 'ready': assert old['sha256'] == current['sha256']
    reasons = list(source['numerical_exclusion_reasons'])
    assert source['numerical_usable'] == (status == 'aligned' and not reasons)
    assert coverage[f'order{order}_native_status'] == status
    assert coverage[f'order{order}_numerical_exclusion_reasons'].split(';') == reasons if reasons else coverage[f'order{order}_numerical_exclusion_reasons'] == ''
    if status != 'aligned':
        assert status in ['input_unavailable', 'native_error', 'parse_error', 'timeout']
        if status not in reasons: reasons.append(status)
        if status == 'input_unavailable': assert any(r['status'] != 'ready' for r in rows)
        mapping = []
    else:
        assert all(r['status'] == 'ready' for r in rows)
        strings = [raw['metrics']['alignment_left'], raw['metrics']['alignment_right']]
        assert len(strings[0]) == len(strings[1]) and all(s.replace('-', '') == r['sequence'] for s, r in zip(strings, rows))
        left = right = 0; mapping = []
        for aa, bb in zip(*strings):
            if aa != '-' and bb != '-': mapping.append((rows[0]['original_positions'][left], rows[1]['original_positions'][right]))
            left += aa != '-'; right += bb != '-'
        assert left == len(rows[0]['sequence']) and right == len(rows[1]['sequence']) and len(mapping) == raw['metrics']['aligned_length']
        assert len(set(a for a, _ in mapping)) == len(set(b for _, b in mapping)) == len(mapping)
    provenance = {k: source[k] for k in ['source_checkpoint', 'source_checkpoint_sha256', 'source_order', 'source_plan_sha256', 'source_input_manifest_sha256', 'source_native_status', 'numerical_usable', 'numerical_exclusion_reasons', 'selected_source']}
    provenance['current_order'] = order; provenance['directed_endpoints'] = [list(e) for e in actual]
    provenance['order_status'] = coverage[f'order{order}_status']; provenance['mapping_exclusions'] = reasons
    return actual, mapping, provenance


def common_residues(ab, ar, br):
    common = [[ar[p], br[p], p] for p in sorted(set(ar) & set(br))]
    cycle = [t for t in common if ab.get(t[0]) == t[1]]
    return common, cycle


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text())
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2 ** 30
    triads, inputs, design, source_counts, bindings = load_design_inputs(plan); bind(bindings, args.plan)
    checkpoints, coverage = load_measured_edges(plan, triads, bindings)
    @lru_cache(maxsize=2048)
    def native(kind, pair, mask, order):
        return recover_native(checkpoints[(kind, pair, mask, order)], coverage[(kind, pair, mask)], inputs, pair, mask, order, bindings)
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False); counts = Counter(); common_count = cycle_count = states = 0; screen_counts = Counter()
    with gzip.open(out / 'common_residue_maps.jsonl.gz', 'xt', compresslevel=1) as handle:
        for number, triad in enumerate(triads, 1):
            ids = [tuple(m) for m in triad['models']]
            assert triad['triad_id'] == hashlib.sha256(json.dumps(ids, separators=(',', ':')).encode()).hexdigest()
            desired = [ids[:2], [ids[2], ids[0]], [ids[2], ids[1]]]
            for mask in ['full', 'plddt70']:
                for orders in itertools.product([0, 1], repeat=3):
                    maps, provenance, edge_flags, edge_exclusions = [], [], [], []
                    for label, kind, ends, order in zip(['ab', 'ar', 'br'], ['primary', 'reference', 'reference'], desired, orders):
                        key = triad['edges'][label]['pair_key']; actual, mapping, trace = native(kind, key, mask, order)
                        assert set(actual) == set(ends)
                        directed = mapping if actual == ends else [(b, a) for a, b in mapping]
                        maps.append(dict(directed)); provenance.append(trace)
                        cov = coverage[(kind, key, mask)]; flags, exclusions = {}, {}
                        for s in plan['screens']:
                            sid = s['id']; flags[sid] = bool(int(cov[sid + '_pass'])); exclusions[sid] = cov[sid + '_exclusions'].split(';') if cov[sid + '_exclusions'] else []
                            assert flags[sid] == (not exclusions[sid])
                        edge_flags.append(flags); edge_exclusions.append(exclusions)
                    common, cycle = common_residues(*maps)
                    reasons = sorted({why for trace in provenance for why in trace['mapping_exclusions']})
                    if not common: reasons.append('no_common_reference_residues')
                    elif not cycle: reasons.append('no_cycle_consistent_residues')
                    passes = {s['id']: all(f[s['id']] for f in edge_flags) for s in plan['screens']}
                    source_excluded = any(t['mapping_exclusions'] for t in provenance)
                    core_status = ['source_excluded' if source_excluded else ('fewer_than_three_common_residues' if len(t) < 3 else 'pending_common_coordinate_geometry') for t in [common, cycle]]
                    value = dict(triad_id=triad['triad_id'], role_order=['a', 'b', 'reference'], models=triad['models'], mask=mask,
                                 orders=list(orders), edge_order=['ab', 'reference_to_a', 'reference_to_b'], edge_provenance=provenance,
                                 original_lengths=[inputs[(*model, 'full')]['original_length'] for model in ids],
                                 retained_input_lengths=[inputs[(*model, mask)]['retained_residues'] for model in ids],
                                 edge_pair_screen_pass=edge_flags, edge_pair_screen_exclusions=edge_exclusions, all_three_pair_screen_pass=passes,
                                 reference_common_triples=common, cycle_consistent_triples=cycle, common_reference_count=len(common), cycle_consistent_count=len(cycle),
                                 common_core_fit_input_status=dict(reference_common=core_status[0], cycle_consistent=core_status[1]),
                                 mapping_exclusions=reasons)
                    handle.write(json.dumps(value, separators=(',', ':')) + '\n'); states += 1; common_count += len(common); cycle_count += len(cycle)
                    counts[mask + ':' + (';'.join(reasons) if reasons else 'recorded_without_edge_exclusions')] += 1
                    for sid, passed in passes.items(): screen_counts[mask + ':' + sid] += int(passed)
            if number % 500 == 0: print('Full triads original-residue maps', number, '/', len(triads), flush=True)
    assert states == plan['expected']['potential_correspondence_states'] == 16 * len(triads)
    # Keep the complete source design (including unscheduled/missing/excluded
    # contexts) by immutable artifact locators rather than repeating it 16 times.
    verify(bindings)
    result = dict(status='complete_full_triad_original_residue_mapping_pending_independent_readback', plan_sha256=sha(args.plan),
                  target_contexts=design['target_contexts'], reference_tie_records=design['reference_tie_records'], duplicate_reference_links=design['duplicate_reference_links'],
                  unique_ordered_model_triads=design['unique_ordered_model_triads'], correspondence_work_triads=len(triads), mask_order_states=states,
                  needed_models=len(inputs) // 2, all_input_dispositions=plan['expected']['all_input_dispositions'], input_source_counts=source_counts,
                  common_reference_residue_occurrences=common_count, cycle_consistent_residue_occurrences=cycle_count, counts=dict(counts), three_pair_pass_state_counts=dict(screen_counts),
                  full_source_design=plan['triad_completion'], full_source_design_sha256=sha(plan['triad_completion']), source_hashes=bindings,
                  artifacts={'common_residue_maps.jsonl.gz': sha(out / 'common_residue_maps.jsonl.gz')}, scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}, indent=2), flush=True)


if __name__ == '__main__': main()
