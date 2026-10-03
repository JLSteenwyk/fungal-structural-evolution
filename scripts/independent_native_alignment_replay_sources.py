"""Closed full recovery and independent inventory for saved-alignment replay."""
from collections import Counter
import json
from pathlib import Path

from independent_baliphy_scalar_sources import load as recovery_sources
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

SUMMARY_FIELDS = ['full_chains', 'checked_chains', 'failed_chains', 'full_groups',
    'complete_groups', 'unresolved_groups', 'saved_alignments', 'candidate_frames',
    'state_observations', 'candidate_anchor_coordinates', 'unanchored_residue_observations',
    'cutoff_count_cells', 'cutoff_counts', 'intact_chains_in_unresolved_groups']


def load(plan, path):
    source, bindings = recovery_sources(plan, path)
    inventory_path = Path(plan['native_inventory']); bind(bindings, inventory_path)
    inventory = json.loads(inventory_path.read_text())
    assert inventory['status'] == 'complete_full_independent_native_alignment_inventory_sample_replay_pending'
    assert inventory['scientific_eligibility'] is False
    completion_path = Path(plan['inventory_completion']); bind(bindings, completion_path)
    completion = json.loads(completion_path.read_text())
    assert completion['status'] == 'verified_original_native_alignment_inventory_completion'
    assert completion['inventory'] == str(inventory_path)
    assert completion['inventory_sha256'] == sha(inventory_path)
    assert completion['native_alignment_samples_replayed'] is False
    assert completion['scientific_eligibility'] is False
    assert completion['original_terminal_handle']['status'] == 'verified_original_terminal_success_with_bound_completed_artifacts'
    for key, digest in [('full_chain_details', 'full_chain_details_sha256'),
                        ('full_source_hash_archive', 'full_source_hash_archive_sha256')]:
        bind(bindings, inventory[key], inventory[digest])
    assert completion['details_sha256'] == inventory['full_chain_details_sha256']
    assert completion['source_hash_archive_sha256'] == inventory['full_source_hash_archive_sha256']
    archive = json.loads(Path(inventory['full_source_hash_archive']).read_text())
    assert archive['status'] == 'full_native_alignment_inventory_source_hash_archive'
    assert len(archive['source_hashes']) == inventory['source_hashes_rechecked'] == completion['bound_sources']
    for p, digest in archive['source_hashes'].items(): bind(bindings, p, digest)
    details = json.loads(Path(inventory['full_chain_details']).read_text())
    assert set(details) == set(source['native']) and len(details) == inventory['full_chains'] == 1620
    lookup = {str(Path(p).resolve()): d for p, d in bindings.items()}
    def doc(p, expected=None):
        p = Path(p); d = lookup[str(p.resolve())]
        assert expected is None or expected == d
        bind(bindings, p, d)
        return json.loads(p.read_text())
    producer = doc(source['complete']['producer_receipt'], source['complete']['producer_receipt_sha256'])
    assert set(producer['states']) == set(details)
    checked = failed = 0
    for cid, entry in details.items():
        native = source['native'][cid]; chain = native['chain']; selected = native['selected_disposition']
        assert entry['chain_id'] == cid and entry['seed'] == chain['seed']
        assert entry['model_input_identity'] == native['model_input_identity']
        assert entry['native_disposition'] == selected['status'] and entry['scientific_eligibility'] is False
        if selected['status'] != 'all_saved_alignments_and_candidate_nodes_checked':
            assert entry['status'] == 'unresolved_failed_native_chain_retained'
            failed += 1; continue
        checked += 1
        assert entry['status'] == 'full_native_source_metadata_and_runtime_labels_checked_sample_decode_pending'
        assert entry['sample_audit'] == selected['sample_audit']
        assert entry['sample_audit_sha256'] == selected['sample_audit_sha256']
        assert entry['input_alignment'] == chain['alignment']
        assert entry['input_alignment_sha256'] == chain['alignment_sha256']
        audit = doc(entry['sample_audit'], entry['sample_audit_sha256'])
        assert entry['selected_attempt_receipt'] == audit['attempt_receipt']
        assert entry['selected_attempt_receipt_sha256'] == audit['attempt_receipt_sha256']
        receipt = doc(audit['attempt_receipt'], audit['attempt_receipt_sha256'])
        assert receipt['exit_code'] == 0 and audit['iterations'] == 1000
        for name in ['C1.P1.fastas', 'runtime-tree.nwk']:
            keys = [p for p in receipt['artifacts'] if Path(p).name == name]; assert len(keys) == 1
            p = Path(audit['attempt_receipt']).parent / keys[0]; value = entry['native_files'][name]
            assert Path(value['path']).resolve() == p.resolve()
            assert value['sha256'] == receipt['artifacts'][keys[0]] == lookup[str(p.resolve())]
            assert p.stat().st_size == value['bytes_at_observation']
        state = producer['states'][cid]; report = doc(state['receipt'], state['receipt_sha256'])
        assert len(report['summaries']) == 1
        summary = report['summaries'][0]
        assert summary['samples'] == entry['expected_saved_alignments'] == 101
        assert summary['source_audit_sha256'] == entry['sample_audit_sha256']
        for key, name in [('original_state_array', 'states.npz'), ('coordinates', 'coordinates.json')]:
            matches = [p for p in summary['artifacts'] if Path(p).name == name]; assert len(matches) == 1
            value = entry[key]
            assert Path(value['path']).resolve() == Path(matches[0]).resolve()
            assert value['sha256'] == summary['artifacts'][matches[0]] == lookup[str(Path(matches[0]).resolve())]
        for key in ['candidate_anchor_coordinates', 'state_observations']:
            assert entry[key] == summary[key]
        assert entry['inherited_unanchored_residue_observations'] == summary['unanchored_residue_observations']
    assert checked == inventory['checked_chains'] == 1618 and failed == inventory['failed_chains'] == 2
    assert inventory['full_groups'] == len(source['groups']) == 405 and inventory['complete_groups'] == 403
    verify(bindings)
    source.update(details=details, inventory=inventory, doc=doc)
    return source, bindings


def summarize(results, source, plan):
    assert set(results) == set(source['details'])
    total = Counter(); cutoff = Counter(); failed = checked = intact_unresolved = 0
    for cid, row in results.items():
        entry = source['details'][cid]
        assert row['chain_id'] == cid and row['seed'] == entry['seed']
        assert row['model_input_identity'] == entry['model_input_identity']
        assert row['scientific_eligibility'] is False
        if entry['status'] == 'unresolved_failed_native_chain_retained':
            assert row['status'] == entry['status'] and row['cutoffs'] == {}
            assert row['saved_alignments'] == 0; failed += 1; continue
        assert row['status'] == 'all_native_saved_alignments_and_projections_replayed_not_posterior_qualification'
        checked += 1
        group = source['groups'][entry['model_input_identity']]
        intact_unresolved += group['status'] == 'unresolved_failed_native_chain_retained'
        for key in ['saved_alignments', 'candidate_frames', 'state_observations',
                    'candidate_anchor_coordinates', 'unanchored_residue_observations', 'cutoff_count_cells']:
            total[key] += row[key]
        assert row['saved_alignments'] == 101 and row['candidate_frames'] == 404
        assert row['state_observations'] == entry['state_observations']
        assert row['candidate_anchor_coordinates'] == entry['candidate_anchor_coordinates']
        assert row['unanchored_residue_observations'] == entry['inherited_unanchored_residue_observations']
        assert set(row['cutoffs']) == {'250', '500'}
        for cut, value in row['cutoffs'].items():
            assert value['retained_samples'] == (75 if cut == '250' else 50)
            assert value['count_cells'] == 22 * row['candidate_anchor_coordinates']
            cutoff[cut] += value['count_cells']
    assert checked == 1618 and failed == 2 and intact_unresolved == 6
    inventory = source['inventory']['totals']
    for key, old in [('saved_alignments', 'saved_alignments_to_decode'), ('candidate_frames', 'candidate_frames_to_decode'),
                     ('state_observations', 'state_observations'), ('candidate_anchor_coordinates', 'candidate_anchor_coordinates'),
                     ('unanchored_residue_observations', 'inherited_unanchored_residue_observations')]:
        assert total[key] == inventory[old]
    assert total['cutoff_count_cells'] == 44 * total['candidate_anchor_coordinates'] == sum(cutoff.values())
    result = dict(full_chains=1620, checked_chains=checked, failed_chains=failed, full_groups=405,
        complete_groups=403, unresolved_groups=2, **dict(total), cutoff_counts=dict(cutoff),
        intact_chains_in_unresolved_groups=intact_unresolved)
    assert all(result[k] == v for k, v in plan['expected'].items())
    return result
