"""Closed fixed selections and schemas for the complete logical-case handoff.

Only source I/O, schemas and identity conventions are shared with the reader.
Membership, counts and aggregation are reconstructed independently there.
"""
import hashlib
import json
from full_screened_balance_sources import load as load_closed_matching

MASKS = ['full', 'plddt70', 'both']
BIT_FIELDS = [f'{role}_{mask}_pass_bits' for mask in MASKS for role in ['target', 'control', 'joint']]
CASE_FIELDS = ['case_id', 'physical_case_id', 'target_id', 'background_id', 'guide',
    'target_family', 'background_family', 'focal_taxon', 'gene_node',
    'target_gene_a', 'target_gene_b', 'background_gene_a', 'background_gene_b',
    'background_taxon_a', 'background_taxon_b', 'target_pair_key', 'background_pair_key',
    'target_same_model', 'background_same_model', 'target_sequence_distance',
    'background_sequence_distance', 'selection_records', 'endpoint_order_bits',
    'policy_bits', 'scenario_bits'] + BIT_FIELDS
SUMMARY_FIELDS = ['target_nodes', 'background_nodes', 'target_policy_records', 'scenarios',
    'scenario_decisions', 'selected_records', 'unmatched_decisions', 'logical_cases',
    'physical_cases', 'unique_selected_targets', 'unique_selected_backgrounds',
    'maximum_selection_reuse', 'case_dispositions', 'guide_census',
    'future_case_mask_rows', 'future_order_pair_cells', 'retained_selection_screen_cells',
    'screens', 'masks', 'guides', 'policies']


def identity(namespace, *values):
    return hashlib.sha256(json.dumps([namespace, *values], separators=(',', ':')).encode()).hexdigest()


def load(plan, path):
    source, bindings = load_closed_matching(plan, path)
    mp = json.loads(open(plan['matching_plan']).read())
    expected = plan['expected']
    for k in ['target_policy_records', 'unmatched_decisions']:
        assert expected[k] == mp['expected'][k]
    assert expected['target_policy_records'] == expected['target_nodes'] * len(plan['policies'])
    assert expected['selected_records'] + expected['unmatched_decisions'] == expected['target_policy_records'] * expected['scenarios']
    return source, bindings
