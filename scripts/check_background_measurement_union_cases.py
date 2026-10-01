#!/usr/bin/env python3
"""Test complete union I/O with synthetic closed sources and rehashed false exports.

This is a software fixture, not native optimization or scientific/journal proof.
The native/numeric/geometry records and prior completion journals are synthetic.
Actual full-source byte binding, checkpoint roles, complete grids, normalization,
nullable states and independent SQL export reconstruction are exercised unchanged.
"""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from background_measurement_union_sources import NUMERIC_INTS, NUMERIC_FLOATS, GEOMETRY_INTS, GEOMETRY_FLOATS
from run_ortholog_pair_guide_comparison import sha


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(rows)


def zipped(path, rows):
    with gzip.open(path, 'wt') as handle:
        for row in rows: handle.write(json.dumps(row) + '\n')


def invoke(command, success=True):
    result = subprocess.run(command, capture_output=True, text=True)
    if success and result.returncode: raise RuntimeError(result.stderr + result.stdout)
    if not success: assert result.returncode != 0, 'False export or incompatible source was accepted'
    return result


def sources(root):
    native, assessment, reuse, old_root = [root / name for name in ['native', 'assessment', 'reuse', 'old']]
    for folder in [native, assessment, reuse, old_root]: folder.mkdir()
    np, sp = root / 'native-plan.json', root / 'reuse-plan.json'
    write(np, dict(output=str(native), expected=dict(full_pairs=3, new_pairs=1)))
    write(sp, dict(output=str(reuse), target_alignment_plan=str(np), expected=dict(full_pairs=3, new_pairs=1)))
    bundle = '1' * 64; pairs, checkpoints, new_rows, reuse_rows, all_checkpoints = [], [], [], [], []
    for index, (a, b, owner) in enumerate([('A', 'B', 'reference_old'), ('A', 'C', 'background_old'), ('B', 'C', '')]):
        ends = [(a, 1), (b, 2)]; pair = hashlib.sha256(json.dumps(sorted(ends), separators=(',', ':')).encode()).hexdigest()
        row = dict(pair_key=pair, model_a=a, version_a=1, model_b=b, version_b=2, measurement_disposition='native_measurement_required' if not owner else 'pending_actual_input_result_and_numeric_reuse_checks')
        pairs.append(row)
        for mask in ['full', 'plddt70']:
            for order in [0, 1]:
                directed = ends if order == 0 else ends[::-1]; source_order = 1 - order if index == 0 else order
                status = 'aligned'
                if index == 2 and (mask != 'full' or order): status = 'native_error' if mask == 'full' else ('parse_error' if order == 0 else 'timeout')
                if index == 1 and mask == 'plddt70': status = 'input_unavailable'
                path = (native if index == 2 else old_root) / f'pairs/{pair[:2]}/{pair}-{mask}-{source_order}.json'
                raw = dict(pair_key=pair, mask=mask, order=source_order, inputs=[dict(model_id=e[0], version=e[1]) for e in directed], status=status,
                           plan_sha256=sha(np) if index == 2 else '2' * 64, input_manifest_sha256=bundle, elapsed_seconds=0.25)
                numeric = geometry = None; reasons = [status]
                if status == 'aligned':
                    length = 2 if index == 0 and mask == 'plddt70' and order == 0 else 10
                    discrepant = index == 0 and mask == 'full' and order == 1
                    nonunique = index == 0 and order == 0
                    numeric = {name: 0.5 for name in NUMERIC_FLOATS}; numeric.update({name: length for name in NUMERIC_INTS})
                    numeric.update(rmsd_status='outside_printed_rounding' if discrepant else 'within_printed_rounding', rmsd_native=0.2, rmsd_recomputed=0.21 if discrepant else 0.2, rmsd_rounding_error=0.01 if discrepant else 0.)
                    geometry = {name: 0.25 for name in GEOMETRY_FLOATS}; geometry.update({name: 3 for name in GEOMETRY_INTS})
                    geometry.update(aligned_length=length, determinant_correction=1, geometry_status='degenerate_at_numeric_tolerance' if nonunique else 'unique_at_numeric_tolerance')
                    raw.update(returncode=0, metrics=dict(aligned_length=length, rmsd=0.2, tm_left=0.5, tm_right=0.5))
                    reasons = (['rmsd_discrepancy'] if discrepant else []) + (['fewer_than_three_pairs'] if length < 3 else []) + (['nonunique_rotation'] if nonunique else [])
                    if index != 2:
                        numeric = dict(pair_key=pair, mask=mask, order=str(source_order), **{k: str(v) for k, v in numeric.items()})
                        geometry = dict(pair_key=pair, mask=mask, order=str(source_order), rmsd_status=numeric['rmsd_status'], **{k: str(v) for k, v in geometry.items()})
                elif status == 'timeout': raw.update(timeout_seconds=30, error='synthetic timeout')
                elif status in ['native_error', 'parse_error']: raw.update(returncode=7 if status == 'native_error' else 0, error='synthetic ' + status)
                write(path, raw); all_checkpoints.append(path)
                current = {k: v for k, v in row.items() if k != 'measurement_disposition'}
                current.update(mask=mask, order=order, selected_source=owner, reuse_status='verified_identical_input_checkpoint_and_retained_disposition' if owner else 'new_native_measurement_pending')
                if owner:
                    current.update(source_checkpoint=str(path), source_checkpoint_sha256=sha(path), source_order=source_order, source_native_status=status,
                                   source_numeric=numeric, source_geometry=geometry, numerical_exclusion_reasons=reasons, numerical_usable=not reasons)
                else:
                    checkpoints.append(dict(path=str(path.relative_to(native)), sha256=sha(path), status=status))
                    new_rows.append(dict(pair_key=pair, mask=mask, order=order, model_left=directed[0][0], version_left=directed[0][1], model_right=directed[1][0], version_right=directed[1][1],
                                         checkpoint_path=str(path), checkpoint_sha256=sha(path), native_status=status, numerical=numeric, geometry=geometry,
                                         numerical_exclusion_reasons=reasons, numerical_usable=not reasons))
                reuse_rows.append(current)
    for folder in [native, assessment, reuse]: table(folder / 'full_background_work_partition.tsv', pairs)
    table(native / 'checkpoint_manifest.tsv', checkpoints); zipped(assessment / 'disposition_measurements.jsonl.gz', new_rows); zipped(reuse / 'background_reuse_dispositions.jsonl.gz', reuse_rows)
    counts = dict(Counter(r['mask'] + ':' + r['native_status'] for r in new_rows))
    nr, ar, qr, rr, ra = native / 'receipt.json', assessment / 'receipt.json', root / 'new-readback.json', reuse / 'receipt.json', reuse / 'readback.json'
    write(nr, dict(status='complete_expanded_background_alignment_dispositions_pending_numeric_geometry', plan_sha256=sha(np), input_bundle_sha256=bundle, counts=counts, directed_dispositions=4, scientific_eligibility=False,
                   artifacts={name: sha(native / name) for name in ['full_background_work_partition.tsv', 'checkpoint_manifest.tsv']}))
    write(ar, dict(status='complete_expanded_background_numeric_geometry_pending_independent_readback', native_receipt_sha256=sha(nr), counts=counts, directed_dispositions=4, scientific_eligibility=False,
                   artifacts={name: sha(assessment / name) for name in ['full_background_work_partition.tsv', 'disposition_measurements.jsonl.gz']}))
    write(qr, dict(status='passed_full_expanded_background_measurement_quaternion_readback', producer_receipt_sha256=sha(ar), counts=counts, directed_dispositions=4, scientific_eligibility=False))
    old_counts = dict(new_native_measurement_pending=4, verified_identical_input_checkpoint_and_retained_disposition=8)
    numerical_counts = dict(Counter('usable' if r['numerical_usable'] else ';'.join(r['numerical_exclusion_reasons']) for r in reuse_rows if r['selected_source']))
    selected = dict(reference_old=4, background_old=4, no_old_source=4)
    write(rr, dict(status='complete_full_expanded_background_reuse_qualification_pending_independent_readback', plan_sha256=sha(sp), counts=old_counts, directed_dispositions=12,
                   selected_source_dispositions=selected, numerical_counts=numerical_counts, scientific_eligibility=False,
                   artifacts={name: sha(reuse / name) for name in ['full_background_work_partition.tsv', 'background_reuse_dispositions.jsonl.gz']}))
    write(ra, dict(status='passed_full_expanded_background_input_checkpoint_numeric_quaternion_reuse_readback', producer_receipt_sha256=sha(rr), counts=old_counts, directed_dispositions=12,
                   selected_source_dispositions=selected, numerical_counts=numerical_counts, scientific_eligibility=False))
    bindings = {str(p): sha(p) for p in [np, sp, nr, ar, qr, rr, ra, *all_checkpoints]}
    nc, oc, na, oa = [root / name for name in ['new-closure.json', 'old-closure.json', 'new-archive.json', 'old-archive.json']]
    ns = dict(full_background_pairs=3, new_pairs=1, directed_dispositions=4, counts=counts)
    os = dict(full_background_pairs=3, directed_dispositions=12, counts=old_counts, selected_source_dispositions=selected, numerical_counts=numerical_counts)
    for path, status, summary, number in [(na, 'complete_verified_full_new_background_measurement_archive', ns, 3), (oa, 'complete_verified_full_background_reuse_disposition_archive', os, 10)]:
        write(path, dict(status=status, summary=summary, source_hashes=bindings, services=[dict(synthetic_fixture_only=True)] * number))
    write(nc, dict(status='complete_verified_full_new_background_native_numeric_geometry', **ns, exact_process_journals_checked=3, bound_source_hashes=len(bindings),
                   full_hash_archive=str(na), full_hash_archive_sha256=sha(na), native_receipt=str(nr), native_receipt_sha256=sha(nr), assessment_receipt=str(ar), assessment_receipt_sha256=sha(ar),
                   independent_readback=str(qr), independent_readback_sha256=sha(qr), scientific_eligibility=False))
    write(oc, dict(status='complete_verified_full_background_input_checkpoint_reuse_with_dispositions', **os, exact_process_journals_checked=10, bound_source_hashes=len(bindings),
                   full_hash_archive=str(oa), full_hash_archive_sha256=sha(oa), producer_receipt=str(rr), producer_receipt_sha256=sha(rr), independent_readback=str(ra), independent_readback_sha256=sha(ra),
                   source_plan=str(sp), source_plan_sha256=sha(sp), scientific_eligibility=False))
    pp = root / 'union-plan.json'; write(pp, dict(full_pairs=3, new_pairs=1, source_pair_counts=dict(reference_old=1, background_old=1, no_old_source=1), native_plan=str(np), native_completion=str(nc), reuse_completion=str(oc), output=str(root / 'union'), pins={}, scope=__doc__))
    return pp, oc, oa


def run():
    with tempfile.TemporaryDirectory(prefix='background-union-fixture-') as temp:
        root = Path(temp); plan, old_closure, old_archive = sources(root); out = root / 'union'
        producer = [sys.executable, 'scripts/union_background_measurements.py', '--plan', str(plan)]
        reader = [sys.executable, 'scripts/readback_background_measurement_union.py', '--plan', str(plan), '--output']
        invoke(producer); invoke(reader + [str(root / 'passed.json')])
        receipt = json.loads((out / 'receipt.json').read_text()); artifact = out / 'background_measurement_dispositions.jsonl.gz'
        artifact_bytes, receipt_bytes = artifact.read_bytes(), (out / 'receipt.json').read_bytes()
        with gzip.open(artifact, 'rt') as handle: original = [json.loads(line) for line in handle]
        assert len(original) == 12 and receipt['source_dispositions'] == dict(reference_old=4, background_old=4, background_new_native=4)
        assert receipt['numerical_counts'] == dict(usable=4, nonunique_rotation=1, rmsd_discrepancy=1, **{'fewer_than_three_pairs;nonunique_rotation': 1}, input_unavailable=2, native_error=1, parse_error=1, timeout=1)
        reversed_old = next(row for row in original if row['selected_source'] == 'reference_old'); assert reversed_old['source_order'] != reversed_old['order']
        labels = ['missing_state', 'duplicate_state', 'wrong_source', 'wrong_source_order', 'reversed_roles', 'wrong_checkpoint', 'original_numeric', 'normalized_numeric', 'original_geometry', 'normalized_geometry', 'cleared_exclusion', 'nullable_promotion', 'changed_native_metric', 'changed_native_error', 'changed_counts']
        rejected = []
        for label in labels:
            rows = copy.deepcopy(original); r = copy.deepcopy(receipt)
            aligned = next(row for row in rows if row['source_native_status'] == 'aligned')
            excluded = next(row for row in rows if row['source_native_status'] == 'aligned' and not row['numerical_usable'])
            missing = next(row for row in rows if row['source_native_status'] == 'input_unavailable')
            if label == 'missing_state': rows.pop()
            elif label == 'duplicate_state': rows.append(copy.deepcopy(rows[0]))
            elif label == 'wrong_source': rows[0]['selected_source'] = 'background_new_native'
            elif label == 'wrong_source_order': rows[0]['source_order'] = 1 - rows[0]['source_order']
            elif label == 'reversed_roles': rows[0]['directed_endpoints'].reverse()
            elif label == 'wrong_checkpoint': rows[0]['source_checkpoint_sha256'] = '0' * 64
            elif label == 'original_numeric': aligned['source_numeric']['rmsd_recomputed'] = 999
            elif label == 'normalized_numeric': aligned['numerical']['rmsd_recomputed'] = 999
            elif label == 'original_geometry': aligned['source_geometry']['relative_rotation_curvature'] = '999'
            elif label == 'normalized_geometry': aligned['geometry']['relative_rotation_curvature'] = 999
            elif label == 'cleared_exclusion': excluded.update(numerical_exclusion_reasons=[], numerical_usable=True)
            elif label == 'nullable_promotion': missing.update(numerical=aligned['numerical'], numerical_usable=True)
            elif label == 'changed_native_metric': aligned['native_metrics']['rmsd'] = 999
            elif label == 'changed_native_error': next(row for row in rows if row['source_native_status'] == 'native_error')['native_returncode'] = 0
            else: r['numerical_counts']['usable'] += 1
            zipped(artifact, rows); r['artifacts'][artifact.name] = sha(artifact); write(out / 'receipt.json', r)
            invoke(reader + [str(root / (label + '.json'))], success=False); rejected.append(label)
        artifact.write_bytes(artifact_bytes); (out / 'receipt.json').write_bytes(receipt_bytes)
        # A coherently rehashed closure with incompatible counts must fail before output creation.
        archive = json.loads(old_archive.read_text()); closure = json.loads(old_closure.read_text())
        closure['counts'] = dict(new_native_measurement_pending=4, verified_identical_input_checkpoint_and_retained_disposition=7, incompatible_inputs_require_new_measurement=1)
        archive['summary']['counts'] = closure['counts']; write(old_archive, archive); closure['full_hash_archive_sha256'] = sha(old_archive); write(old_closure, closure)
        conf = json.loads(plan.read_text()); conf['output'] = str(root / 'incompatible-union'); write(plan, conf)
        failed = invoke(producer, success=False); assert 'Incompatible reuse states' in failed.stderr and not Path(conf['output']).exists()
        return dict(status='passed_complete_background_measurement_union_software_cases', full_pairs=3, directed_dispositions=12, source_dispositions=receipt['source_dispositions'],
                    numerical_counts=receipt['numerical_counts'], rejected_rehashed_exports=rejected, incompatible_closed_source_rejected_before_output=True,
                    synthetic_native_numeric_geometry_and_prior_journals=True, production_journal_closure_tested=False, scope=__doc__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path)
    result = run(); result['script_sha256'] = sha(__file__)
    args = parser.parse_args()
    if args.output:
        with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
