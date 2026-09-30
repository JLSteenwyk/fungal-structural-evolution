"""Source I/O and proof lineage for the full reference measurement union."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def bind(bindings, path, digest=None):
    path = str(path)
    digest = sha(path) if digest is None else digest
    assert path not in bindings or bindings[path] == digest, path
    bindings[path] = digest


def verify(bindings):
    for path, digest in bindings.items(): assert sha(path) == digest, path


def artifacts(bindings, root, receipt):
    bind(bindings, root / 'receipt.json')
    for name, digest in receipt['artifacts'].items(): bind(bindings, root / name, digest)


def keyed_table(path):
    rows = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['pair_key'], row['mask'], int(row['order'])
            assert key not in rows
            rows[key] = row
    return rows


def pair_ends(row):
    return [(row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))]


def native_bundle(plan, bindings):
    """Bind complete new-native output and prior numeric/quaternion proofs."""
    np = json.loads(Path(plan['native_plan']).read_text())
    dp = json.loads(Path(plan['native_diagnostic_plan']).read_text())
    gp = json.loads(Path(plan['native_geometry_plan']).read_text())
    qp = json.loads(Path(plan['native_geometry_readback_plan']).read_text())
    root, diagnostic, geometry = [Path(p['output']) for p in [np, dp, gp]]
    assert dp['source_plan'] == plan['native_plan'] and gp['diagnostic_plan'] == plan['native_diagnostic_plan']
    assert qp['source_plan'] == plan['native_geometry_plan'] and qp['output'] == plan['native_readback']
    nr, dr, gr = [json.loads((p / 'receipt.json').read_text()) for p in [root, diagnostic, geometry]]
    qr = json.loads(Path(plan['native_readback']).read_text())
    assert nr['status'] == 'complete_reference_alignment_dispositions_pending_readback'
    assert dr['status'] == 'complete_reference_alignment_rmsd_diagnostic_not_scientific_acceptance'
    assert gr['status'] == 'complete_reference_alignment_geometry_pending_independent_readback'
    assert qr['status'] == 'passed_full_reference_geometry_readback'
    for p, r, config in [(root, nr, plan['native_plan']), (diagnostic, dr, plan['native_diagnostic_plan']),
                         (geometry, gr, plan['native_geometry_plan'])]:
        assert r['plan_sha256'] == sha(config)
        artifacts(bindings, p, r)
    assert qr['plan_sha256'] == sha(plan['native_geometry_readback_plan'])
    assert dr['producer_receipt_sha256'] == sha(root / 'receipt.json')
    assert gr['diagnostic_receipt_sha256'] == sha(diagnostic / 'receipt.json')
    assert qr['producer_receipt_sha256'] == sha(geometry / 'receipt.json')
    assert qr['alignments_checked'] == gr['alignments'] == dr['numerically_checked_alignments']
    assert qr['counts'] == gr['counts'] and nr['counts'] == dr['counts']
    assert nr['directed_dispositions'] == dr['directed_dispositions'] == 4 * plan['new_pairs']
    assert nr['distinct_model_pairs'] == plan['new_pairs'] and nr['full_reference_pairs'] == plan['full_pairs']
    assert nr['existing_catalog_pairs_pending_reuse'] == dr['existing_catalog_pairs_pending_reuse'] == plan['full_pairs'] - plan['new_pairs']
    numeric = keyed_table(diagnostic / 'numeric_readback.tsv')
    geometries = keyed_table(geometry / 'alignment_geometry.tsv')
    assert len(numeric) == len(geometries) == gr['alignments'] and set(numeric) == set(geometries)
    full, new = {}, {}
    with (root / 'full_reference_work_partition.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            pair = row['pair_key']; ends = pair_ends(row)
            assert pair == hashlib.sha256(json.dumps(sorted(ends), separators=(',', ':')).encode()).hexdigest()
            assert pair not in full and ends[0] != ends[1]
            full[pair] = ends
            if row['measurement_disposition'] == 'native_measurement_required': new[pair] = ends
            else: assert row['measurement_disposition'] == 'pending_actual_input_result_and_numeric_reuse_checks'
    assert len(full) == plan['full_pairs'] and len(new) == plan['new_pairs']
    assert (diagnostic / 'full_reference_work_partition.tsv').read_bytes() == (root / 'full_reference_work_partition.tsv').read_bytes()
    checkpoints = {}
    with (root / 'checkpoint_manifest.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            pair, mask, order = Path(row['path']).stem.rsplit('-', 2)
            key = pair, mask, int(order)
            assert key not in checkpoints and pair in new and mask in ['full', 'plddt70'] and int(order) in [0, 1]
            assert row['path'] == f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json'
            checkpoints[key] = dict(path=str(root / row['path']), sha256=row['sha256'], status=row['status'])
            bind(bindings, root / row['path'], row['sha256'])
    assert set(checkpoints) == {(p, m, o) for p in new for m in ['full', 'plddt70'] for o in [0, 1]}
    assert dict(Counter(k[1] + ':' + v['status'] for k, v in checkpoints.items())) == nr['counts']
    assert set(numeric) == {k for k, v in checkpoints.items() if v['status'] == 'aligned'}
    for key in numeric:
        assert numeric[key]['aligned_length'] == geometries[key]['aligned_length']
        assert numeric[key]['rmsd_status'] == geometries[key]['rmsd_status']
    for config in [plan['native_plan'], plan['native_diagnostic_plan'], plan['native_geometry_plan'], plan['native_geometry_readback_plan']]:
        bind(bindings, config)
    for path, digest in np['pins'].items(): bind(bindings, path, digest)
    for path, digest in nr['upstream_bindings'].items(): bind(bindings, path, digest)
    bind(bindings, plan['native_readback'])
    return dict(full=full, new=new, numeric=numeric, geometry=geometries, checkpoints=checkpoints,
                plan=np, root=root, native_receipt=nr, diagnostic_receipt=dr, geometry_receipt=gr)


def load_union_sources(plan):
    bindings = dict(plan['pins'])
    reuse = json.loads(Path(plan['reuse_completion']).read_text())
    assert reuse['status'] == 'complete_verified_full_reference_input_checkpoint_numeric_reuse'
    assert reuse['full_reference_pairs'] == plan['full_pairs'] and reuse['directed_dispositions'] == 4 * plan['full_pairs']
    assert reuse['counts'] == {'new_native_measurement_pending': 4 * plan['new_pairs'],
                              'verified_identical_input_checkpoint_and_retained_disposition': 4 * (plan['full_pairs'] - plan['new_pairs'])}
    bind(bindings, reuse['full_hash_archive'], reuse['full_hash_archive_sha256'])
    archive = json.loads(Path(reuse['full_hash_archive']).read_text())
    assert archive['summary']['directed_dispositions'] == reuse['directed_dispositions'] and len(archive['services']) == 2
    bind(bindings, reuse['producer_receipt'], reuse['producer_receipt_sha256'])
    bind(bindings, reuse['independent_readback'], reuse['independent_readback_sha256'])
    rp = json.loads(Path(reuse['producer_receipt']).read_text())
    ap = json.loads(Path(reuse['independent_readback']).read_text())
    assert rp['status'] == 'complete_full_reference_reuse_qualification_pending_independent_readback'
    assert ap['status'] == 'passed_full_reference_input_checkpoint_numeric_reuse_readback'
    assert ap['producer_receipt_sha256'] == reuse['producer_receipt_sha256']
    assert ap['counts'] == rp['counts'] == reuse['counts']
    artifacts(bindings, Path(reuse['producer_receipt']).parent, rp)
    reuse_plan = json.loads(Path(reuse['source_plan']).read_text())
    assert reuse_plan['target_alignment_plan'] == plan['native_plan']
    bind(bindings, reuse['source_plan'], reuse['source_plan_sha256'])
    native_closure = json.loads(Path(plan['native_completion']).read_text())
    assert native_closure['status'] == 'complete_verified_new_reference_native_alignment_numeric_geometry'
    assert native_closure['new_pairs'] == plan['new_pairs'] and native_closure['full_pairs'] == plan['full_pairs']
    bind(bindings, native_closure['full_hash_archive'], native_closure['full_hash_archive_sha256'])
    closed = json.loads(Path(native_closure['full_hash_archive']).read_text())
    assert len(closed['services']) == 4 and closed['summary']['directed_dispositions'] == 4 * plan['new_pairs']
    native = native_bundle(plan, bindings)
    assert closed['source_hashes'][str(native['root'] / 'receipt.json')] == sha(native['root'] / 'receipt.json')
    assert closed['source_hashes'][plan['native_readback']] == sha(plan['native_readback'])
    for path in [plan['reuse_completion'], plan['native_completion']]: bind(bindings, path)
    reuse_path = Path(reuse['producer_receipt']).parent / 'reference_reuse_dispositions.jsonl'
    verify(bindings)
    return reuse_path, native, bindings
