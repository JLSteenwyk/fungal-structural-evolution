"""Closed source I/O and explicit schemas for the expanded covariance handoff."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
from background_measurement_union_sources import closed_source
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

CASE_FIELDS = ['case_id', 'physical_case_id', 'target_id', 'background_id', 'guide',
    'target_family', 'background_family', 'family_component', 'species_pattern_id',
    'species_pattern_row', 'focal_index', 'background_index_a', 'background_index_b',
    'selection_records', 'target_same_model', 'background_same_model']
INCIDENCE_FIELDS = ['case_id', 'entity_kind', 'entity_id', 'entity_value', 'side',
    'endpoint', 'signed_loading', 'unsigned_loading']
PATTERN_FIELDS = ['species_pattern_id', 'row_index', 'taxa', 'twice_weights']
FAMILY_FIELDS = ['family', 'family_component', 'component_families']
SUMMARY_FIELDS = ['logical_cases', 'physical_cases', 'selected_records', 'species_columns',
    'patterns', 'zero_patterns', 'rank', 'families', 'family_components',
    'entity_occurrences', 'unique_entities', 'trees', 'covariance_entries_per_tree']


def identity(kind, value):
    return hashlib.sha256(json.dumps(['expanded-covariance-' + kind + '-v1', value],
        separators=(',', ':')).encode()).hexdigest()


def load(plan, path):
    bindings = dict(plan['pins']); bind(bindings, path)
    completion = closed_source(plan['case_completion'],
        'complete_verified_full_matching_logical_case_index',
        'complete_verified_full_matching_logical_case_index_archive', 2, bindings)
    config = json.loads(Path(plan['case_plan']).read_text())
    root = Path(config['output']); rp = root / 'receipt.json'
    receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_full_matching_case_index_pending_independent_readback'
    assert completion['scientific_eligibility'] is receipt['scientific_eligibility'] is False
    assert completion['producer_receipt'] == str(rp)
    assert completion['producer_receipt_sha256'] == bindings[str(rp)] == sha(rp)
    assert receipt['plan_sha256'] == sha(plan['case_plan'])
    for key in ['logical_cases', 'physical_cases', 'selected_records']:
        assert receipt[key] == completion[key] == plan['expected'][key]
    for name, digest in receipt['artifacts'].items(): bind(bindings, root / name, digest)
    matching = json.loads(Path(config['matching_plan']).read_text())
    graph = Path(matching['graph'])
    with gzip.open(root / 'case_index.tsv.gz', 'rt') as handle:
        cases = list(csv.DictReader(handle, delimiter='\t'))
    assert len(cases) == len({r['case_id'] for r in cases}) == receipt['logical_cases']
    wanted = {role: {r[role + '_id'] for r in cases} for role in ['target', 'background']}
    nodes = {}
    for role in wanted:
        source = graph / (role + '_nodes.jsonl'); assert str(source) in bindings
        nodes[role] = {}
        with source.open() as handle:
            for line in handle:
                row = json.loads(line)
                if row['node_id'] in wanted[role]:
                    assert row['node_id'] not in nodes[role]; nodes[role][row['node_id']] = row
        assert set(nodes[role]) == wanted[role]
    for case in cases:
        assert case['target_family'] == case['background_family']
        for role in nodes:
            node = nodes[role][case[role + '_id']]
            for field, source in [('family', 'family'), ('pair_key', 'pair_key'),
                ('gene_a', 'gene_a'), ('gene_b', 'gene_b'), ('same_model', 'same_model')]:
                assert case[role + '_' + field] == str(node[source])
            assert case['guide'] == node['guide']
        assert case['focal_taxon'] == nodes['target'][case['target_id']]['taxon_id']
        for end in ['a', 'b']:
            assert case['background_taxon_' + end] == nodes['background'][case['background_id']]['taxon_' + end]
    kp = Path(plan['kernel_plan']); kernel_plan = json.loads(kp.read_text())
    kernel = Path(plan['kernel_root']); krp = kernel / 'receipt.json'
    kr = json.loads(krp.read_text()); audit = json.loads(Path(plan['kernel_readback']).read_text())
    assert kr['status'] == 'complete_species_distance_kernels_pending_independent_readback'
    assert kr['plan_sha256'] == sha(kp)
    assert audit['status'] == 'passed_full_species_distance_kernel_readback'
    assert audit['source_receipt_sha256'] == sha(krp)
    assert set(kr['trees']) == set(kernel_plan['trees']) == set(plan['trees'])
    for p in [kp, krp, plan['kernel_readback']]: bind(bindings, p)
    for p, digest in kernel_plan['pins'].items(): bind(bindings, p, digest)
    for name, digest in kr['artifacts'].items(): bind(bindings, kernel / name, digest)
    for tree, spec in kernel_plan['trees'].items():
        bind(bindings, spec['tree'], kr['trees'][tree]['tree_sha256'])
        bind(bindings, spec['audit'], kr['trees'][tree]['audit_sha256'])
        assert json.loads(Path(spec['audit']).read_text())['status'] == spec['audit_status']
    taxa = json.loads((kernel / 'taxa.json').read_text())
    assert len(taxa) == len(set(taxa)) == plan['expected']['species_columns']
    assert all(r[k] in taxa for r in cases for k in ['focal_taxon', 'background_taxon_a', 'background_taxon_b'])
    verify(bindings)
    return dict(cases=cases, nodes=nodes, taxa=taxa, kernel=kernel,
        kernel_plan=kernel_plan, case_receipt=receipt), bindings
