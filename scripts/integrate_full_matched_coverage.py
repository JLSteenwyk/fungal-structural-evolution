#!/usr/bin/env python3
"""Keep all frozen selected/unmatched controls while applying complete pair coverage."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import shutil
from collections import Counter
from pathlib import Path
from full_matched_coverage_sources import load, MASKS, ATTRITION_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def quality(path, screens, expected_pairs):
    records = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['pair_key'], row['mask']; assert key not in records and row['mask'] in MASKS[:2]
            assert all(row[s['id'] + '_pass'] in ['0', '1'] and (row[s['id'] + '_pass'] == '1') == (not row[s['id'] + '_exclusions']) for s in screens)
            records[key] = row
    pairs = {key[0] for key in records}; assert len(pairs) == expected_pairs and set(records) == {(p, m) for p in pairs for m in MASKS[:2]}
    return records


def project(node, records, screens):
    ends = [(node['model_id_' + side], node['version_' + side]) for side in ['a', 'b']]
    assert node['pair_key'] == hashlib.sha256(json.dumps(sorted(ends), separators=(',', ':')).encode()).hexdigest()
    assert int(node['same_model']) == int(ends[0] == ends[1])
    if ends[0] == ends[1]: return dict(disposition='identical_model_no_alignment', bits=[0, 0, 0])
    bits = []
    for mask in MASKS[:2]:
        row = records[node['pair_key'], mask]
        by_endpoint = {(row['model_' + s], int(row['version_' + s])): int(row['length_' + s]) for s in ['a', 'b']}
        assert set(by_endpoint) == set(ends) and all(by_endpoint[end] == node['length_' + side] for end, side in zip(ends, ['a', 'b']))
        bits.append(sum(1 << i for i, spec in enumerate(screens) if row[spec['id'] + '_pass'] == '1'))
    return dict(disposition='distinct_model_pair', bits=[*bits, bits[0] & bits[1]])


def extras(t, b, pt, pb):
    result = dict(guide=t['guide'], family=t['family'], focal_taxon=t['taxon_id'], gene_node=t['gene_node'], target_pair_key=t['pair_key'], background_pair_key=b['pair_key'], target_comparison_disposition=pt['disposition'], background_comparison_disposition=pb['disposition'], target_sequence_distance=t['sequence_distance'], background_sequence_distance=b['sequence_distance'])
    for i, mask in enumerate(MASKS):
        result.update({f'target_{mask}_pass_bits': pt['bits'][i], f'control_{mask}_pass_bits': pb['bits'][i], f'joint_{mask}_pass_bits': pt['bits'][i] & pb['bits'][i]})
    return result


def run(plan_path):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text()); source, bindings = load(plan, plan_path)
    tq, bq = [quality(source[key], plan['screens'], plan['expected'][count]) for key, count in [('target_quality', 'target_pairs'), ('background_quality', 'background_pairs')]]
    tp = {key: project(n, tq, plan['screens']) for key, n in source['targets'].items()}; bp = {key: project(n, bq, plan['screens']) for key, n in source['backgrounds'].items()}; del tq, bq
    dispositions = Counter(label + ':' + source[which][key]['guide'] + ':' + value['disposition'] for label, which, records in [('target', 'targets', tp), ('background', 'backgrounds', bp)] for key, value in records.items())
    out = Path(plan['output']); assert shutil.disk_usage(out.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30; out.mkdir(exist_ok=False)
    scenario_ids = [s['scenario_id'] for s in source['scenarios']]; all_scenarios = set(scenario_ids); seen = set(); baseline, baseline_pass, selected_counts, categories = Counter(), Counter(), Counter(), Counter(); selected = unmatched = 0
    with gzip.open(source['selection'] / 'selections.tsv.gz', 'rt') as original, (source['selection'] / 'target_policy_selection_status.tsv').open() as status, gzip.open(out / 'selected_pair_coverage.tsv.gz', 'wt', compresslevel=1) as selected_out, gzip.open(out / 'target_policy_coverage_status.tsv.gz', 'wt', compresslevel=1) as status_out:
        selections = csv.DictReader(original, delimiter='\t'); groups = iter(itertools.groupby(selections, key=lambda row: (row['target_id'], row['policy']))); current = next(groups, None)
        statuses = csv.DictReader(status, delimiter='\t'); sw = uw = None
        for index, raw in enumerate(statuses, 1):
            tid, policy = raw['target_id'], raw['policy']; key = tid, policy
            assert key not in seen and tid in tp and policy in plan['policies']; seen.add(key); t = source['targets'][tid]; p = tp[tid]; guide = t['guide']; assert guide in plan['guides']
            matched = raw['matched_scenarios'].split(',') if raw['matched_scenarios'] else []; missing = raw['unmatched_scenarios'].split(',') if raw['unmatched_scenarios'] else []
            assert len(set(matched)) == len(matched) and len(set(missing)) == len(missing) and not set(matched) & set(missing) and set(matched + missing) == all_scenarios
            baseline[guide, policy] += 1; unmatched += len(missing)
            for m, mask in enumerate(MASKS):
                for j, spec in enumerate(plan['screens']): baseline_pass[guide, policy, mask, spec['id']] += bool(p['bits'][m] & (1 << j))
            exported_status = dict(raw, guide=guide, family=t['family'], focal_taxon=t['taxon_id'], gene_node=t['gene_node'], target_pair_key=t['pair_key'], target_comparison_disposition=p['disposition'], **{f'target_{mask}_pass_bits': p['bits'][m] for m, mask in enumerate(MASKS)})
            if uw is None: uw = csv.DictWriter(status_out, list(exported_status), delimiter='\t', lineterminator='\n'); uw.writeheader()
            uw.writerow(exported_status)
            observed = set()
            if current is not None and current[0] == key:
                for chosen in current[1]:
                    sid, bid = chosen['scenario_id'], chosen['background_id']; assert sid in matched and sid not in observed and bid in bp and chosen['endpoint_order'] in ['0', '1']; observed.add(sid)
                    b, pb = source['backgrounds'][bid], bp[bid]; assert b['guide'] == guide
                    row = dict(chosen, **extras(t, b, p, pb))
                    if sw is None: sw = csv.DictWriter(selected_out, list(row), delimiter='\t', lineterminator='\n'); sw.writeheader()
                    sw.writerow(row); selected += 1; selected_counts[guide, policy, sid] += 1
                    for m, mask in enumerate(MASKS):
                        for j, spec in enumerate(plan['screens']):
                            category = int(bool(p['bits'][m] & (1 << j))) * 2 + int(bool(pb['bits'][m] & (1 << j)))
                            categories[guide, policy, sid, mask, spec['id'], category] += 1
                current = next(groups, None)
            assert observed == set(matched)
            if index % 100000 == 0: print('Full matched target/policy coverage records', index, '/', plan['expected']['target_policy_records'], flush=True)
        assert current is None
    assert seen == {(tid, p) for tid in tp for p in plan['policies']}
    assert len(seen) == plan['expected']['target_policy_records'] and selected == plan['expected']['selected_records'] and unmatched == plan['expected']['unmatched_decisions']
    # A full design without selected records still needs a well-defined empty export.
    if selected == 0:
        with gzip.open(out / 'selected_pair_coverage.tsv.gz', 'wt') as handle: handle.write('target_id\tpolicy\tscenario_id\tbackground_id\tendpoint_order\tscore\teligible_candidates\tequal_score_candidates\tscore_gap_to_second\t' + '\t'.join(extras(next(iter(source['targets'].values())), next(iter(source['backgrounds'].values())), next(iter(tp.values())), next(iter(bp.values())))) + '\n')
    rows = 0
    with (out / 'matched_attrition_counts.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, ATTRITION_FIELDS, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for guide, policy, sid, mask, spec in itertools.product(plan['guides'], plan['policies'], scenario_ids, MASKS, plan['screens']):
            n = baseline[guide, policy]; k = selected_counts[guide, policy, sid]; values = [categories[guide, policy, sid, mask, spec['id'], i] for i in range(4)]; assert sum(values) == k
            source_prefix = '|'.join([guide, policy, sid]); mr = source['matching_receipt']
            assert n == mr['counts'].get(source_prefix + '|targets', 0) and k == mr['counts'].get(source_prefix + '|matched', 0) and n - k == mr['counts'].get(source_prefix + '|unmatched', 0)
            target_all = baseline_pass[guide, policy, mask, spec['id']]; target_matched = values[2] + values[3]; assert 0 <= target_all - target_matched <= n - k
            writer.writerow(dict(guide=guide, policy=policy, scenario_id=sid, mask=mask, screen=spec['id'], all_target_records=n, matched_records=k, unmatched_records=n-k, target_pass_all_records=target_all, target_pass_matched_records=target_matched, target_pass_unmatched_records=target_all-target_matched, control_pass_matched_records=values[1]+values[3], joint_pass_matched_records=values[3], target_only_pass_matched_records=values[2], control_only_pass_matched_records=values[1], neither_pass_matched_records=values[0])); rows += 1
    verify(bindings)
    result = dict(status='complete_full_matched_coverage_pending_independent_readback', plan_sha256=sha(plan_path), target_nodes=len(tp), background_nodes=len(bp), target_policy_records=len(seen), scenarios=len(scenario_ids), scenario_decisions=len(seen)*len(scenario_ids), selected_records=selected, unmatched_decisions=unmatched, selection_screen_cells=selected*len(MASKS)*len(plan['screens']), full_scenario_screen_cells=len(seen)*len(scenario_ids)*len(MASKS)*len(plan['screens']), attrition_rows=rows, screens=plan['screens'], masks=MASKS, guides=plan['guides'], policies=plan['policies'], node_dispositions=dict(dispositions), source_hashes=bindings, artifacts={name:sha(out/name) for name in ['selected_pair_coverage.tsv.gz','target_policy_coverage_status.tsv.gz','matched_attrition_counts.tsv']}, scientific_eligibility=False, scope=plan['scope'])
    with (out/'receipt.json').open('x') as handle: handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2));return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);run(parser.parse_args().plan)
