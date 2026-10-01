"""Closed full matching exports for descriptive screening balance and reuse."""
import csv
import gzip
import json
from pathlib import Path
from background_measurement_union_sources import closed_source
from full_matched_coverage_sources import SUMMARY_FIELDS as MATCHED_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

KEY = ['guide', 'policy', 'scenario_id', 'mask', 'screen']
MASKS = ['full', 'plddt70', 'both']
REUSE_FIELDS = KEY + ['background_id', 'background_pair_key', 'retained_target_records', 'reciprocal_node_reuse_weight', 'retained_records_using_physical_pair', 'reciprocal_physical_pair_reuse_weight']
SUMMARY_FIELDS = ['target_nodes', 'background_nodes', 'selected_records', 'scenarios', 'strata', 'coverage_rows', 'balance_rows', 'reuse_rows', 'retained_selection_screen_cells', 'features', 'screens', 'masks', 'guides', 'policies']


def load(plan, plan_path):
    bindings = dict(plan['pins']); bind(bindings, plan_path)
    closed = closed_source(plan['matching_completion'], 'complete_verified_full_fixed_matched_coverage_attrition', 'complete_verified_full_matched_coverage_archive', 2, bindings)
    assert closed['scientific_eligibility'] is False and closed['source_plan'] == plan['matching_plan'] and closed['source_plan_sha256'] == sha(plan['matching_plan'])
    rp, ap = [Path(closed[k]) for k in ['producer_receipt', 'independent_readback']]
    r, a = [json.loads(p.read_text()) for p in [rp, ap]]
    assert bindings[str(rp)] == closed['producer_receipt_sha256'] and bindings[str(ap)] == closed['independent_readback_sha256']
    assert r['status'] == 'complete_full_matched_coverage_pending_independent_readback' and a['status'] == 'passed_full_matched_coverage_sql_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['matching_plan']) and a['producer_receipt_sha256'] == sha(rp)
    assert all(r[k] == a[k] == closed[k] for k in MATCHED_FIELDS)
    mp = json.loads(Path(plan['matching_plan']).read_text()); assert rp.parent == Path(mp['output'])
    for key in ['target_nodes', 'background_nodes', 'selected_records', 'scenarios']:
        assert r[key] == mp['expected'][key] == plan['expected'][key]
    assert r['screens'] == mp['screens'] == plan['screens'] and r['masks'] == MASKS and r['guides'] == plan['guides'] and r['policies'] == plan['policies']
    for name, digest in r['artifacts'].items(): bind(bindings, rp.parent / name, digest)
    graph = Path(mp['graph']); nodes = {}
    for kind, name in [('target', 'target_nodes.jsonl'), ('background', 'background_nodes.jsonl')]:
        path = graph / name; assert str(path) in bindings; index = {}
        with path.open() as handle:
            for line in handle:
                n = json.loads(line); assert n['node_id'] not in index; index[n['node_id']] = n
        assert len(index) == plan['expected'][kind + '_nodes']; nodes[kind] = index
    attrs = {}
    with (rp.parent / 'matched_attrition_counts.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = tuple(row[k] for k in KEY); assert key not in attrs; attrs[key] = row
    selection = Path(mp['selection']); scenarios = json.loads((selection / 'scenarios.json').read_text()); assert str(selection / 'scenarios.json') in bindings
    expected = {(g, p, s['scenario_id'], m, screen['id']) for g in plan['guides'] for p in plan['policies'] for s in scenarios for m in MASKS for screen in plan['screens']}
    assert len(scenarios) == plan['expected']['scenarios'] and set(attrs) == expected and len(attrs) == r['attrition_rows']
    verify(bindings)
    return dict(nodes=nodes, root=rp.parent, scenarios=scenarios, attrition=attrs, prior_receipt=rp), bindings


def target_flags(source, policies):
    """Read closed whole target/policy export; never infer eligibility from selection."""
    flags = {}; seen = set()
    with gzip.open(source['root'] / 'target_policy_coverage_status.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            tid, policy = row['target_id'], row['policy']; assert (tid, policy) not in seen and policy in policies; seen.add((tid, policy))
            node = source['nodes']['target'][tid]; values = [int(row['target_' + m + '_pass_bits']) for m in MASKS]
            assert all(0 <= b < 64 for b in values) and values[2] == values[0] & values[1]
            assert row['guide'] == node['guide'] and row['target_pair_key'] == node['pair_key']
            if node['same_model']: assert values == [0, 0, 0] and row['target_comparison_disposition'] == 'identical_model_no_alignment'
            else: assert row['target_comparison_disposition'] == 'distinct_model_pair'
            assert tid not in flags or flags[tid] == values; flags[tid] = values
    assert seen == {(tid, p) for tid in source['nodes']['target'] for p in policies}
    return flags
