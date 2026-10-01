"""Complete closed new-native and actual old-result sources for the background union."""
import csv
import json
from collections import Counter
from pathlib import Path
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

NUMERIC_INTS = ['aligned_length', 'joint_plddt70_pairs']
NUMERIC_FLOATS = ['rmsd_recomputed', 'rmsd_native', 'rmsd_rounding_error', 'sequence_identity_exact', 'tm_left_native', 'tm_right_native', 'coverage_left', 'coverage_right', 'joint_plddt70_fraction']
GEOMETRY_INTS = ['aligned_length', 'rank_left', 'rank_right', 'determinant_correction']
GEOMETRY_FLOATS = ['rms_width1_left', 'rms_width2_left', 'rms_width3_left', 'rms_width1_right', 'rms_width2_right', 'rms_width3_right', 'width2_to_width1_left', 'width2_to_width1_right', 'cross_s1', 'cross_s2', 'cross_s3', 'minimum_rotation_curvature', 'relative_rotation_curvature', 'relative_numeric_tolerance']
SUMMARY_FIELDS = ['full_pairs', 'directed_dispositions', 'source_dispositions', 'native_status_counts', 'numerical_counts']


def closed_source(path, status, archive_status, journals, bindings):
    path = Path(path); c = json.loads(path.read_text()); assert c['status'] == status and c['exact_process_journals_checked'] == journals
    bind(bindings, path); bind(bindings, c['full_hash_archive'], c['full_hash_archive_sha256'])
    archive = json.loads(Path(c['full_hash_archive']).read_text()); assert archive['status'] == archive_status and len(archive['services']) == journals
    assert len(archive['source_hashes']) == c['bound_source_hashes']
    assert all(c.get(key) == value for key, value in archive['summary'].items())
    for p, digest in archive['source_hashes'].items(): bind(bindings, p, digest)
    return c


def load_sources(plan, plan_path):
    bindings = dict(plan['pins']); bind(bindings, plan_path)
    new = closed_source(plan['native_completion'], 'complete_verified_full_new_background_native_numeric_geometry', 'complete_verified_full_new_background_measurement_archive', 3, bindings)
    old = closed_source(plan['reuse_completion'], 'complete_verified_full_background_input_checkpoint_reuse_with_dispositions', 'complete_verified_full_background_reuse_disposition_archive', 10, bindings)
    np, rp = Path(plan['native_plan']), Path(old['source_plan']); nc, rc = [json.loads(p.read_text()) for p in [np, rp]]
    assert rc['target_alignment_plan'] == str(np) and old['source_plan_sha256'] == sha(rp)
    assert new['new_pairs'] == nc['expected']['new_pairs'] == rc['expected']['new_pairs'] == plan['new_pairs']
    assert new['full_background_pairs'] == old['full_background_pairs'] == nc['expected']['full_pairs'] == rc['expected']['full_pairs'] == plan['full_pairs']
    assert old['counts'] == {'new_native_measurement_pending': 4 * plan['new_pairs'], 'verified_identical_input_checkpoint_and_retained_disposition': 4 * (plan['full_pairs'] - plan['new_pairs'])}, 'Incompatible reuse states require an additional closed native workload before full union'
    assert old['directed_dispositions'] == 4 * plan['full_pairs'] and new['directed_dispositions'] == 4 * plan['new_pairs']
    assert old['selected_source_dispositions'] == {name: 4 * count for name, count in plan['source_pair_counts'].items()}
    for closure, fields in [(new, [('native_receipt', 'native_receipt_sha256'), ('assessment_receipt', 'assessment_receipt_sha256'), ('independent_readback', 'independent_readback_sha256')]),
                            (old, [('producer_receipt', 'producer_receipt_sha256'), ('independent_readback', 'independent_readback_sha256')])]:
        for path_field, digest_field in fields:
            assert bindings.get(closure[path_field]) == closure[digest_field]
            bind(bindings, closure[path_field], closure[digest_field])
    nr = json.loads(Path(new['native_receipt']).read_text()); ar = json.loads(Path(new['assessment_receipt']).read_text()); qr = json.loads(Path(new['independent_readback']).read_text())
    rr = json.loads(Path(old['producer_receipt']).read_text()); ra = json.loads(Path(old['independent_readback']).read_text())
    assert nr['status'] == 'complete_expanded_background_alignment_dispositions_pending_numeric_geometry' and nr['plan_sha256'] == sha(np)
    assert ar['status'] == 'complete_expanded_background_numeric_geometry_pending_independent_readback' and ar['native_receipt_sha256'] == sha(new['native_receipt'])
    assert qr['status'] == 'passed_full_expanded_background_measurement_quaternion_readback' and qr['producer_receipt_sha256'] == sha(new['assessment_receipt'])
    assert rr['status'] == 'complete_full_expanded_background_reuse_qualification_pending_independent_readback' and rr['plan_sha256'] == sha(rp)
    assert ra['status'] == 'passed_full_expanded_background_input_checkpoint_numeric_quaternion_reuse_readback' and ra['producer_receipt_sha256'] == sha(old['producer_receipt'])
    assert rr['counts'] == ra['counts'] == old['counts'] and rr['directed_dispositions'] == ra['directed_dispositions'] == old['directed_dispositions']
    assert nr['directed_dispositions'] == ar['directed_dispositions'] == qr['directed_dispositions'] == new['directed_dispositions']
    assert nr['counts'] == ar['counts'] == qr['counts'] == new['counts']
    for name in ['numerical_counts', 'selected_source_dispositions']:
        assert rr[name] == ra[name] == old[name]
    assert new['scientific_eligibility'] is old['scientific_eligibility'] is False
    for record in [nr, ar, qr, rr, ra]: assert record['scientific_eligibility'] is False
    native_root, geometry_root, reuse_root = [Path(p).parent for p in [new['native_receipt'], new['assessment_receipt'], old['producer_receipt']]]
    assert native_root == Path(nc['output']) and reuse_root == Path(rc['output'])
    for root, record in [(native_root, nr), (geometry_root, ar), (reuse_root, rr)]:
        for name, digest in record['artifacts'].items(): bind(bindings, root / name, digest)
    ledger = native_root / 'full_background_work_partition.tsv'
    assert ledger.read_bytes() == (geometry_root / ledger.name).read_bytes() == (reuse_root / ledger.name).read_bytes()
    full, new_pairs = {}, set()
    with ledger.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            pair = row['pair_key']; assert pair not in full; full[pair] = [(row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))]
            if row['measurement_disposition'] == 'native_measurement_required': new_pairs.add(pair)
            else: assert row['measurement_disposition'] == 'pending_actual_input_result_and_numeric_reuse_checks'
    assert len(full) == plan['full_pairs'] and len(new_pairs) == plan['new_pairs']
    checkpoints = {}
    with (native_root / 'checkpoint_manifest.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            pair, mask, order = Path(row['path']).stem.rsplit('-', 2); key = pair, mask, int(order)
            assert key not in checkpoints and row['path'] == f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json'
            checkpoints[key] = dict(path=str(native_root / row['path']), sha256=row['sha256'], status=row['status']); bind(bindings, native_root / row['path'], row['sha256'])
    assert set(checkpoints) == {(pair, mask, order) for pair in new_pairs for mask in ['full', 'plddt70'] for order in [0, 1]}
    assert dict(Counter(key[1] + ':' + row['status'] for key, row in checkpoints.items())) == nr['counts']
    bind(bindings, np); bind(bindings, rp); verify(bindings)
    return dict(full=full, new_pairs=new_pairs, checkpoints=checkpoints, native_plan_sha256=sha(np), native_bundle=nr['input_bundle_sha256'],
                new_measurements=geometry_root / 'disposition_measurements.jsonl.gz', reuse_ledger=reuse_root / 'background_reuse_dispositions.jsonl.gz', ledger=ledger), bindings
