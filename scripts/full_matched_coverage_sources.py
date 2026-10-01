"""Whole closed coverage and frozen matching inputs for full matched attrition."""
import json
from pathlib import Path
from background_measurement_union_sources import closed_source
from full_matched_graph_sources import load_graph
from full_background_coverage_sources import SUMMARY_FIELDS as COVERAGE_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

MASKS = ['full', 'plddt70', 'both']
ATTRITION_FIELDS = ['guide', 'policy', 'scenario_id', 'mask', 'screen', 'all_target_records', 'matched_records', 'unmatched_records', 'target_pass_all_records', 'target_pass_matched_records', 'target_pass_unmatched_records', 'control_pass_matched_records', 'joint_pass_matched_records', 'target_only_pass_matched_records', 'control_only_pass_matched_records', 'neither_pass_matched_records']
SUMMARY_FIELDS = ['target_nodes', 'background_nodes', 'target_policy_records', 'scenarios', 'scenario_decisions', 'selected_records', 'unmatched_decisions', 'selection_screen_cells', 'full_scenario_screen_cells', 'attrition_rows', 'screens', 'masks', 'guides', 'policies', 'node_dispositions']


def load(plan, plan_path):
    bindings = dict(plan['pins']); bind(bindings, plan_path)
    bc = closed_source(plan['background_completion'], 'complete_verified_full_background_original_length_coverage', 'complete_verified_full_background_coverage_archive', 2, bindings)
    assert bc['scientific_eligibility'] is False and bc['source_plan'] == plan['background_plan'] and bc['source_plan_sha256'] == sha(plan['background_plan'])
    brp, bap = [Path(bc[k]) for k in ['producer_receipt', 'independent_readback']]; br, ba = [json.loads(p.read_text()) for p in [brp, bap]]
    assert bindings[str(brp)] == bc['producer_receipt_sha256'] and bindings[str(bap)] == bc['independent_readback_sha256']
    assert br['status'] == 'complete_full_background_coverage_pending_independent_readback' and ba['status'] == 'passed_full_background_coverage_sql_decimal_readback'
    assert ba['producer_receipt_sha256'] == sha(brp) and br['plan_sha256'] == ba['plan_sha256'] == sha(plan['background_plan'])
    assert all(br[k] == ba[k] == bc[k] for k in COVERAGE_FIELDS) and bc['pairs'] == plan['expected']['background_pairs']
    bp = json.loads(Path(plan['background_plan']).read_text()); assert brp.parent == Path(bp['output'])
    tc = Path(plan['target_completion']); target = json.loads(tc.read_text()); bind(bindings, tc)
    assert target['status'] == 'complete_verified_full_expanded_pair_event_and_taxon_coverage_screens' and len(target['services']) == 2 and target['scientific_eligibility'] is False
    for path, digest in target['source_hashes'].items(): bind(bindings, path, digest)
    tp = json.loads(Path(plan['target_plan']).read_text()); trp = Path(tp['output']) / 'receipt.json'; tr = json.loads(trp.read_text()); assert str(trp) in bindings
    assert tr['status'] == 'complete_full_expanded_duplication_coverage_pending_independent_readback' and tr['plan_sha256'] == sha(plan['target_plan'])
    assert tr['pairs'] == target['summary']['pairs'] == plan['expected']['target_pairs']
    for name in ['pair_mask_rows', 'pair_pass_counts', 'pair_both_masks_pass_counts']: assert tr[name] == target['summary'][name]
    assert tr['screens'] == tp['screens'] == br['screens'] == plan['screens'] and br['target_pairs'] == tr['pairs']
    assert br['target_pass_counts'] == tr['pair_pass_counts'] and br['target_both_masks_pass_counts'] == tr['pair_both_masks_pass_counts']
    mc = Path(plan['matching_completion']); matching = json.loads(mc.read_text()); bind(bindings, mc)
    assert matching['status'] == 'complete_verified_expanded_background_fixed_matching' and len(matching['services']) == 2 and matching['scientific_eligibility'] is False
    for path, digest in matching['source_hashes'].items(): bind(bindings, path, digest)
    selection = Path(plan['selection']); mrp = selection / 'receipt.json'; mr = json.loads(mrp.read_text()); assert str(mrp) in bindings
    assert mr['status'] == 'complete_metadata_background_control_selection_pending_readback'
    for key in ['target_policy_records', 'scenarios', 'selected_records', 'source_edges']: assert mr[key] == matching['summary'][key] == plan['expected'][key]
    assert matching['summary']['scenario_decisions'] == plan['expected']['target_policy_records'] * plan['expected']['scenarios']
    assert matching['summary']['unmatched_decisions'] == plan['expected']['unmatched_decisions'] == matching['summary']['scenario_decisions'] - mr['selected_records']
    for root, record in [(brp.parent, br), (trp.parent, tr), (selection, mr)]:
        for name, digest in record['artifacts'].items(): bind(bindings, root / name, digest)
    graph = Path(plan['graph']); gr = load_graph(plan, bindings, mr)
    for name, digest in gr['artifacts'].items(): bind(bindings, graph / name, digest)
    cp = Path(plan['input_census_completion']); census_closed = json.loads(cp.read_text()); bind(bindings, cp)
    assert census_closed['status'] == 'complete_verified_full_fixed_matching_node_physical_identity_census'
    assert len(census_closed['services']) == 1 and census_closed['scientific_eligibility'] is False
    for path, digest in census_closed['source_hashes'].items(): bind(bindings, path, digest)
    crp = Path(plan['input_census']); assert str(crp) in census_closed['source_hashes']
    census = json.loads(crp.read_text())
    assert census['status'] == 'passed_complete_fixed_matching_node_physical_identity_census'
    assert census['plan_sha256'] == sha(plan['input_census_plan'])
    for field, expected in [('target_nodes', 'target_nodes'), ('background_nodes', 'background_nodes'), ('target_physical_pairs', 'target_pairs'), ('background_physical_pairs', 'background_pairs'), ('source_inventory_models', 'source_inventory_models')]:
        assert census[field] == census_closed['summary'][field] == plan['expected'][expected]
    assert census['used_target_pairs'] == plan['expected']['target_pairs'] and census['used_background_pairs'] == plan['expected']['background_pairs']
    assert not census['confidence_descriptor_mismatch_endpoint_occurrences'] and not census['maximum_confidence_descriptor_absolute_differences']
    for name in ['receipt.json', 'target_nodes.jsonl', 'background_nodes.jsonl']:
        path = str(graph / name); assert census['source_hashes'][path] == bindings[path]
    nodes = []
    for name, count in [('target_nodes.jsonl', 'target_nodes'), ('background_nodes.jsonl', 'background_nodes')]:
        records = {}
        with (graph / name).open() as handle:
            for line in handle:
                raw = json.loads(line); assert raw['node_id'] not in records; records[raw['node_id']] = raw
        assert len(records) == gr[count] == plan['expected'][count]; nodes.append(records)
    scenarios = json.loads((selection / 'scenarios.json').read_text()); assert len(scenarios) == mr['scenarios'] and len({s['scenario_id'] for s in scenarios}) == len(scenarios)
    assert plan['expected']['target_policy_records'] == len(nodes[0]) * len(plan['policies'])
    for path in [plan['background_plan'], plan['target_plan']]: bind(bindings, path)
    verify(bindings)
    return dict(targets=nodes[0], backgrounds=nodes[1], target_quality=trp.parent / 'pair_mask_coverage.tsv', background_quality=brp.parent / 'pair_mask_coverage.tsv', selection=selection, matching_receipt=mr, scenarios=scenarios), bindings
