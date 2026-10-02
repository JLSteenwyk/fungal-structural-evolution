#!/usr/bin/env python3
"""Index every frozen selection without collapsing genes onto physical models."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import itertools
import json
from pathlib import Path
import shutil
from full_matching_case_sources import load, identity, MASKS, BIT_FIELDS, CASE_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def metadata(target, background):
    assert target['guide'] == background['guide']
    for node in [target, background]:
        ends = sorted((node['model_id_' + s], node['version_' + s]) for s in ['a', 'b'])
        import hashlib
        assert node['pair_key'] == hashlib.sha256(json.dumps(ends, separators=(',', ':')).encode()).hexdigest()
        assert int(node['same_model']) == int(ends[0] == ends[1])
    return dict(case_id=identity('fixed-matched-logical-case-v1', target['node_id'], background['node_id']),
        physical_case_id=identity('fixed-matched-physical-case-v1', target['pair_key'], background['pair_key']),
        target_id=target['node_id'], background_id=background['node_id'], guide=target['guide'],
        target_family=target['family'], background_family=background['family'], focal_taxon=target['taxon_id'], gene_node=target['gene_node'],
        target_gene_a=target['gene_a'], target_gene_b=target['gene_b'], background_gene_a=background['gene_a'], background_gene_b=background['gene_b'],
        background_taxon_a=background['taxon_a'], background_taxon_b=background['taxon_b'],
        target_pair_key=target['pair_key'], background_pair_key=background['pair_key'],
        target_same_model=target['same_model'], background_same_model=background['same_model'],
        target_sequence_distance=target['sequence_distance'], background_sequence_distance=background['sequence_distance'])


def run(path):
    path = Path(path); plan = json.loads(path.read_text()); source, bindings = load(plan, path)
    out = Path(plan['output']); assert shutil.disk_usage(out.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    out.mkdir(exist_ok=True)
    lock = (out / 'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (out / 'receipt.json').exists(), 'Completed stage cannot be restarted'
    ph = sha(path); marker = out / 'stage_plan.json'
    state = dict(plan_sha256=ph, schema='fixed-matched-logical-case-v1', recovery='Full deterministic replay; no partial census is accepted.')
    if marker.exists(): assert json.loads(marker.read_text()) == state, 'Different plan cannot reuse output'
    else: marker.write_text(json.dumps(state, indent=2) + '\n')
    scenarios = {r['scenario_id']: i for i, r in enumerate(source['scenarios'])}
    policies = {p: i for i, p in enumerate(plan['policies'])}
    cases = {}; selected = unmatched = statuses = 0; seen = set(); retained = 0
    links = out / 'selection_case_links.tsv.gz.tmp'; casefile = out / 'case_index.tsv.gz.tmp'
    with gzip.open(source['root'] / 'selected_pair_coverage.tsv.gz', 'rt') as raw, gzip.open(source['root'] / 'target_policy_coverage_status.tsv.gz', 'rt') as st, gzip.open(links, 'wt', compresslevel=1) as export:
        original = csv.DictReader(raw, delimiter='\t'); grouped = iter(itertools.groupby(original, key=lambda r: (r['target_id'], r['policy']))); current = next(grouped, None)
        writer = csv.DictWriter(export, ['source_row_ordinal', 'case_id', 'physical_case_id'] + original.fieldnames, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for row in csv.DictReader(st, delimiter='\t'):
            tid, policy = row['target_id'], row['policy']; key = tid, policy
            assert key not in seen and policy in policies; seen.add(key); statuses += 1
            target = source['nodes']['target'][tid]; assert row['guide'] == target['guide']
            assert row['family'] == target['family'] and row['focal_taxon'] == target['taxon_id'] and row['gene_node'] == target['gene_node'] and row['target_pair_key'] == target['pair_key']
            matched = row['matched_scenarios'].split(',') if row['matched_scenarios'] else []
            missing = row['unmatched_scenarios'].split(',') if row['unmatched_scenarios'] else []
            assert len(set(matched)) == len(matched) and len(set(missing)) == len(missing)
            assert not set(matched) & set(missing) and set(matched + missing) == set(scenarios)
            unmatched += len(missing); observed = set()
            if current is not None and current[0] == key:
                for raw_selection in current[1]:
                    sid, bid = raw_selection['scenario_id'], raw_selection['background_id']
                    assert sid in matched and sid not in observed and raw_selection['endpoint_order'] in ['0', '1']; observed.add(sid)
                    background = source['nodes']['background'][bid]; casekey = tid, bid
                    if casekey not in cases:
                        cases[casekey] = dict(metadata(target, background), selection_records=0, endpoint_order_bits=0, policy_bits=0, scenario_bits=0, **{f:int(raw_selection[f]) for f in BIT_FIELDS})
                    case = cases[casekey]
                    assert all(case[f] == int(raw_selection[f]) for f in BIT_FIELDS)
                    assert raw_selection['guide'] == case['guide'] and raw_selection['target_pair_key'] == case['target_pair_key'] and raw_selection['background_pair_key'] == case['background_pair_key']
                    assert raw_selection['target_sequence_distance'] == str(target['sequence_distance']) and raw_selection['background_sequence_distance'] == str(background['sequence_distance'])
                    for mask in MASKS:
                        t, b, j = [case[f'{role}_{mask}_pass_bits'] for role in ['target', 'control', 'joint']]
                        assert 0 <= t < 64 and 0 <= b < 64 and j == t & b
                        assert t == int(row[f'target_{mask}_pass_bits'])
                        retained += j.bit_count()
                    for role in ['target', 'control', 'joint']:
                        assert case[f'{role}_both_pass_bits'] == case[f'{role}_full_pass_bits'] & case[f'{role}_plddt70_pass_bits']
                    for role, n in [('target', target), ('control', background)]:
                        if n['same_model']: assert all(case[f'{role}_{m}_pass_bits'] == 0 for m in MASKS)
                    case['selection_records'] += 1; case['endpoint_order_bits'] |= 1 << int(raw_selection['endpoint_order'])
                    case['policy_bits'] |= 1 << policies[policy]; case['scenario_bits'] |= 1 << scenarios[sid]
                    selected += 1
                    writer.writerow(dict(source_row_ordinal=selected, case_id=case['case_id'], physical_case_id=case['physical_case_id'], **raw_selection))
                current = next(grouped, None)
            assert observed == set(matched)
            if statuses % 100000 == 0: print('full_matching_case_statuses', statuses, 'selections', selected, 'logical_cases', len(cases), flush=True)
        assert current is None
    assert seen == {(tid,p) for tid in source['nodes']['target'] for p in policies}
    assert statuses == plan['expected']['target_policy_records'] and selected == plan['expected']['selected_records'] and unmatched == plan['expected']['unmatched_decisions']
    assert retained == sum(int(r['joint_pass_matched_records']) for r in source['attrition'].values())
    with gzip.open(casefile, 'wt', compresslevel=1) as handle:
        writer = csv.DictWriter(handle, CASE_FIELDS, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for key in sorted(cases): writer.writerow(cases[key])
    disposition = Counter(); guide_census = {}
    for guide in plan['guides']:
        group = [r for r in cases.values() if r['guide'] == guide]
        guide_census[guide] = dict(logical_cases=len(group), physical_cases=len({r['physical_case_id'] for r in group}), selected_records=sum(r['selection_records'] for r in group), target_nodes=len({r['target_id'] for r in group}), background_nodes=len({r['background_id'] for r in group}), target_families=len({r['target_family'] for r in group}), background_families=len({r['background_family'] for r in group}), focal_taxa=len({r['focal_taxon'] for r in group}))
    for r in cases.values(): disposition[f"target_same_model={r['target_same_model']},background_same_model={r['background_same_model']}"] += 1
    summary = dict(target_nodes=len(source['nodes']['target']), background_nodes=len(source['nodes']['background']), target_policy_records=statuses, scenarios=len(scenarios), scenario_decisions=statuses*len(scenarios), selected_records=selected, unmatched_decisions=unmatched,
        logical_cases=len(cases), physical_cases=len({r['physical_case_id'] for r in cases.values()}), unique_selected_targets=len({t for t,b in cases}), unique_selected_backgrounds=len({b for t,b in cases}), maximum_selection_reuse=max((r['selection_records'] for r in cases.values()),default=0), case_dispositions=dict(disposition), guide_census=guide_census,
        future_case_mask_rows=2*len(cases), future_order_pair_cells=8*len(cases), retained_selection_screen_cells=retained, screens=plan['screens'], masks=MASKS, guides=plan['guides'], policies=plan['policies'])
    verify(bindings)
    # Replaying after interruption replaces only this stage's unaccepted artifacts.
    links.replace(out / 'selection_case_links.tsv.gz'); casefile.replace(out / 'case_index.tsv.gz')
    receipt = dict(status='complete_full_matching_case_index_pending_independent_readback', plan_sha256=ph, **summary, source_hashes=bindings,
        artifacts={n:sha(out/n) for n in ['stage_plan.json', 'selection_case_links.tsv.gz', 'case_index.tsv.gz']},
        original_unmatched_status_path=str(source['root']/'target_policy_coverage_status.tsv.gz'), original_unmatched_status_sha256=sha(source['root']/'target_policy_coverage_status.tsv.gz'), scientific_eligibility=False, scope=plan['scope'])
    with (out/'receipt.json').open('x') as handle: handle.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='source_hashes'},indent=2),flush=True)
    return receipt


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
