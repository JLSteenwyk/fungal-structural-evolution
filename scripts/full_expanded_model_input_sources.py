"""Closed source I/O and schemas for complete expanded working-model inputs."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import pyarrow as pa
from background_measurement_union_sources import closed_source
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

ORDERS = ['00', '01', '10', '11', 'mean']
MASKS = ['full', 'plddt70']
GATES = ['mask', 'both_masks']
DISTANCE = ['gene_distance_delta', 'gene_distance_square_delta', 'gene_distance_cube_delta']
IDENTITY = ['aligned_identity_delta', 'aligned_identity_square_delta', 'aligned_identity_cube_delta']
NUISANCE = ['original_coverage_delta', 'log_aligned_length_ratio', 'aligned_plddt70_fraction_delta',
    'log_original_length_ratio', 'mean_ca_plddt_delta', 'fraction_ca_plddt_below50_delta']
NUMERIC = ['rmsd_delta', 'native_tm_dissimilarity_delta'] + DISTANCE + IDENTITY + NUISANCE
INTEGER = ['case_row', 'species_pattern_row', 'numerical_usable', 'joint_mask_pass_bits', 'joint_both_pass_bits', 'selection_records']
FIELDS = ['input_id', 'row_identity', 'case_id', 'case_row', 'physical_case_id', 'target_id', 'background_id',
    'guide', 'target_family', 'background_family', 'focal_taxon', 'gene_node',
    'target_gene_a', 'target_gene_b', 'background_gene_a', 'background_gene_b',
    'background_taxon_a', 'background_taxon_b', 'target_pair_key', 'background_pair_key',
    'family_component', 'species_pattern_id', 'species_pattern_row', 'mask', 'order_contrast',
    'numerical_usable', 'disposition', 'joint_mask_pass_bits', 'joint_both_pass_bits',
    'selection_records', 'target_state_keys', 'background_state_keys'] + NUMERIC
STRATUM = ['guide', 'policy', 'scenario_id', 'mask', 'eligibility_gate', 'screen', 'order_contrast']
COUNT_FIELDS = STRATUM + ['all_target_records', 'unmatched_records', 'selected_records', 'same_model_records', 'quality_retained_records',
    'quality_excluded_records', 'quality_retained_numeric_records', 'unique_targets', 'unique_background_nodes',
    'unique_background_physical_pairs', 'unique_target_families', 'unique_family_components',
    'unique_species_patterns', 'unique_focal_taxa', 'maximum_background_node_reuse',
    'maximum_background_physical_pair_reuse', 'sum_inverse_background_node_reuse_weights',
    'sum_inverse_background_physical_pair_reuse_weights', 'sum_inverse_family_component_weights']
SUMMARY_FIELDS = ['logical_cases', 'physical_cases', 'selected_records', 'unmatched_decisions',
    'case_mask_rows', 'model_input_rows', 'partitions', 'settings', 'valid_input_rows',
    'invalid_input_rows', 'quality_retained_setting_occurrences', 'model_specifications_per_setting',
    'future_model_setting_rows', 'working_tree_alternatives']
CATALOG_FIELDS = ['role', 'pair_key', 'mask', 'order', 'numerical_usable', 'native_status',
    'rmsd_recomputed', 'tm_left_native', 'tm_right_native', 'sequence_identity_exact',
    'joint_plddt70_fraction', 'original_coverage_a', 'original_coverage_b', 'aligned_length',
    'model_a', 'version_a', 'model_b', 'version_b', 'original_length_a', 'original_length_b']


def schema():
    return pa.schema([(f, pa.float64() if f in NUMERIC else pa.int64() if f in INTEGER else pa.string()) for f in FIELDS])


def identity(row_identity, order):
    return hashlib.sha256(json.dumps(['expanded-model-input-v1', row_identity, order], separators=(',', ':')).encode()).hexdigest()


def read_table(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def load(plan, path):
    bindings = dict(plan['pins']); bind(bindings, path); stages = {}
    definitions = [('cases', 'complete_verified_full_matching_logical_case_index'),
        ('measurements', 'complete_verified_full_expanded_matched_case_measurements'),
        ('catalog', 'complete_verified_full_expanded_directed_measurement_catalog'),
        ('covariance', 'complete_verified_full_expanded_covariance')]
    for name, status in definitions:
        c = closed_source(plan[name + '_completion'], status, status + '_archive', 2, bindings)
        assert c['scientific_eligibility'] is False
        config_path = plan[name + '_plan']; config = json.loads(Path(config_path).read_text()); bind(bindings, config_path)
        root = Path(config['output']); rp = root / 'receipt.json'; assert c['producer_receipt'] == str(rp)
        assert bindings[str(rp)] == c['producer_receipt_sha256'] == sha(rp)
        r = json.loads(rp.read_text()); assert r['plan_sha256'] == sha(config_path) and r['scientific_eligibility'] is False
        for n, d in r['artifacts'].items(): bind(bindings, root / n, d)
        stages[name] = dict(root=root, config=config, completion=c)
        del r
    assert stages['measurements']['config']['case_index_completion'] == plan['cases_completion']
    assert stages['measurements']['config']['catalog_completion'] == plan['catalog_completion']
    assert stages['covariance']['config']['case_completion'] == plan['cases_completion']
    for name in ['cases', 'measurements', 'covariance']:
        for k in ['logical_cases', 'physical_cases', 'selected_records']:
            assert stages[name]['completion'][k] == plan['expected'][k]
    assert stages['measurements']['completion']['unmatched_decisions'] == plan['expected']['unmatched_decisions']
    assert stages['catalog']['completion']['directed_states'] == plan['expected']['directed_states']
    cases = read_table(stages['cases']['root'] / 'case_index.tsv.gz')
    assert len(cases) == len({c['case_id'] for c in cases}) == plan['expected']['logical_cases']
    for key in ['guides', 'policies', 'screens']: assert stages['cases']['config'][key] == plan[key]
    cov = {r['case_id']: r for r in read_table(stages['covariance']['root'] / 'case_covariance_index.tsv.gz')}
    assert set(cov) == {c['case_id'] for c in cases}
    for case in cases:
        for f in ['physical_case_id', 'target_id', 'background_id', 'guide', 'target_family', 'background_family', 'selection_records']:
            assert cov[case['case_id']][f] == case[f]
    matching_path = stages['cases']['config']['matching_plan']; matching = json.loads(Path(matching_path).read_text())
    assert stages['catalog']['config']['matching_plan'] == matching_path
    nodes = {}; wanted = {r: {c[r + '_id'] for c in cases} for r in ['target', 'background']}
    for role in wanted:
        source = Path(matching['graph']) / (role + '_nodes.jsonl'); assert str(source) in bindings
        nodes[role] = {}
        with source.open() as handle:
            for line in handle:
                node = json.loads(line)
                if node['node_id'] in wanted[role]:
                    assert node['node_id'] not in nodes[role]; nodes[role][node['node_id']] = node
        assert set(nodes[role]) == wanted[role]
    wanted_pairs = {r: {c[r + '_pair_key'] for c in cases if c[r + '_same_model'] == '0'} for r in wanted}
    catalog = {}; count = 0
    with gzip.open(stages['catalog']['root'] / 'directed_measurements.tsv.gz', 'rt') as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            count += 1
            if r['pair_key'] in wanted_pairs[r['role']]:
                key = r['role'], r['pair_key'], r['mask'], int(r['order']); assert key not in catalog
                catalog[key] = {f: r[f] for f in CATALOG_FIELDS}
    assert count == plan['expected']['directed_states']
    assert set(catalog) == {(role, pair, m, o) for role, pairs in wanted_pairs.items() for pair in pairs for m in MASKS for o in [0, 1]}
    matching_root = Path(matching['output']); attrs = {}
    for r in read_table(matching_root / 'matched_attrition_counts.tsv'):
        key = tuple(r[f] for f in ['guide', 'policy', 'scenario_id', 'mask', 'screen']); assert key not in attrs; attrs[key] = r
    assert str(matching_root / 'matched_attrition_counts.tsv') in bindings
    scenarios_path = Path(matching['selection']) / 'scenarios.json'; assert str(scenarios_path) in bindings
    scenarios = json.loads(scenarios_path.read_text()); assert len(scenarios) == plan['expected']['scenarios']
    verify(bindings)
    return dict(stages=stages, cases=cases, cov=cov, nodes=nodes, catalog=catalog, attrition=attrs, scenarios=scenarios), bindings
