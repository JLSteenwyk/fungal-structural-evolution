"""Closed full mapping/input source I/O for same-residue geometry stages."""
import gzip
import hashlib
import itertools
import json
from pathlib import Path
from full_triad_common_sources import load_design_inputs
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

DEFINITIONS = [('reference_common', 'reference_common_triples'),
               ('cycle_consistent', 'cycle_consistent_triples')]
METRICS = ['rmsd_ab', 'geometry_ab', 'relative_curvature_ab', 'sequence_identity_ab',
           'rmsd_ar', 'geometry_ar', 'relative_curvature_ar', 'sequence_identity_ar',
           'rmsd_br', 'geometry_br', 'relative_curvature_br', 'sequence_identity_br',
           'rmsd_ar_minus_br', 'mean_plddt_a', 'mean_plddt_b',
           'mean_plddt_reference', 'joint_plddt70_fraction', 'fit_status']
BASE_FIELDS = ['triad_id', 'model_a', 'version_a', 'model_b', 'version_b',
               'model_reference', 'version_reference', 'mask', 'order_ab',
               'order_ar', 'order_br', 'mapping_definition', 'triples_sha256',
               'common_residues', 'length_a', 'length_b', 'length_reference',
               'retained_a', 'retained_b', 'retained_reference', 'coverage_a',
               'coverage_b', 'coverage_reference', 'retained_coverage_a',
               'retained_coverage_b', 'retained_coverage_reference',
               'source_exclusions', 'mapping_disagreement_count']
SUMMARY_FIELDS = ['target_contexts', 'reference_tie_records',
                  'duplicate_reference_links', 'correspondence_work_triads',
                  'mapping_states', 'fit_rows', 'counts', 'core_screen_pass_counts',
                  'screen_pass_counts', 'triple_occurrences', 'fitted_pdb_inputs']


def fields(plan):
    result = BASE_FIELDS + METRICS
    for s in plan['screens']:
        result += [s['id'] + suffix for suffix in ['_core_pass', '_core_exclusions',
                   '_three_pair_pass', '_three_pair_exclusions', '_pass', '_exclusions']]
    return result


def triple_sha(triples):
    return hashlib.sha256(json.dumps(triples, separators=(',', ':')).encode()).hexdigest()


def load_sources(plan, plan_path):
    bindings = dict(plan['pins']); bind(bindings, plan_path)
    mp = json.loads(Path(plan['mapping_plan']).read_text())
    root = Path(mp['output']); rp = root / 'receipt.json'
    cp = Path(plan['mapping_completion']); c = json.loads(cp.read_text())
    assert c['status'] == 'complete_verified_full_triad_original_residue_mapping'
    assert c['exact_process_journals_checked'] == 2
    assert c['producer_receipt'] == str(rp) and c['producer_receipt_sha256'] == sha(rp)
    bind(bindings, cp); bind(bindings, c['full_hash_archive'], c['full_hash_archive_sha256'])
    archive = json.loads(Path(c['full_hash_archive']).read_text())
    assert archive['status'] == 'complete_verified_full_triad_original_residue_mapping_archive'
    assert len(archive['services']) == 2 and len(archive['source_hashes']) == c['bound_source_hashes']
    for path, digest in archive['source_hashes'].items(): bind(bindings, path, digest)
    r = json.loads(rp.read_text()); ap = Path(c['independent_readback']); a = json.loads(ap.read_text())
    assert a['status'] == 'passed_full_triad_original_residue_mapping_readback'
    assert r['status'] == 'complete_full_triad_original_residue_mapping_pending_independent_readback'
    assert a['producer_receipt_sha256'] == sha(rp) and sha(ap) == c['independent_readback_sha256']
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['mapping_plan'])
    for key in ['target_contexts', 'reference_tie_records', 'duplicate_reference_links',
                'correspondence_work_triads', 'mask_order_states']:
        assert r[key] == a[key] == c[key] == plan['expected'][key]
    assert r['artifacts']['common_residue_maps.jsonl.gz'] == sha(root / 'common_residue_maps.jsonl.gz')
    assert archive['source_hashes'][str(root / 'common_residue_maps.jsonl.gz')] == r['artifacts']['common_residue_maps.jsonl.gz']
    assert plan['screens'] == mp['screens']
    triads, inputs, design, counts, source_bindings = load_design_inputs(mp)
    for path, digest in source_bindings.items(): bind(bindings, path, digest)
    assert len(triads) == c['correspondence_work_triads']
    verify(bindings)
    return triads, inputs, r, root, bindings


def iterate_maps(root, triads, inputs):
    """Exhaust the closed full ordered source grid; do not preselect coverage."""
    with gzip.open(root / 'common_residue_maps.jsonl.gz', 'rt') as handle:
        for triad in triads:
            models = triad['models']
            for mask in ['full', 'plddt70']:
                for orders in itertools.product([0, 1], repeat=3):
                    line = next(handle, None); assert line is not None, 'Missing mapping state'
                    row = json.loads(line)
                    assert (row['triad_id'], row['mask'], row['orders']) == (triad['triad_id'], mask, list(orders))
                    assert row['models'] == models and row['role_order'] == ['a', 'b', 'reference']
                    assert row['original_lengths'] == [inputs[(*model, 'full')]['original_length'] for model in models]
                    assert row['retained_input_lengths'] == [inputs[(*model, mask)]['retained_residues'] for model in models]
                    source_excluded = any(e['mapping_exclusions'] for e in row['edge_provenance'])
                    for definition, field in DEFINITIONS:
                        triples = row[field]; n = len(triples)
                        assert all(len(t) == 3 and all(isinstance(p, int) for p in t) for t in triples)
                        for column, model in enumerate(models):
                            positions = [t[column] for t in triples]
                            assert len(set(positions)) == n and set(positions) <= set(inputs[(*model, mask)]['original_positions'])
                        expected_status = 'source_excluded' if source_excluded else ('fewer_than_three_common_residues' if n < 3 else 'pending_common_coordinate_geometry')
                        assert row['common_core_fit_input_status'][definition] == expected_status
                    assert len(row['reference_common_triples']) == row['common_reference_count']
                    assert len(row['cycle_consistent_triples']) == row['cycle_consistent_count']
                    assert {tuple(t) for t in row['cycle_consistent_triples']} <= {tuple(t) for t in row['reference_common_triples']}
                    yield row
        assert next(handle, None) is None, 'Extra mapping state'
