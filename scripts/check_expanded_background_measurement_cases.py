#!/usr/bin/env python3
"""Real-USalign software checks for the four-collection background measurement path."""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import math
import tempfile
from collections import Counter
from pathlib import Path
from duplication_alignment_inputs import render_ca
from expanded_background_native_handoff import load_handoff
from run_duplication_alignments import run_job
import run_expanded_background_alignments as runner
import assess_expanded_background_measurements as assessment
import readback_expanded_background_measurements as reader
from run_ortholog_pair_guide_comparison import sha


def save(path, value): path.write_text(json.dumps(value, indent=2) + '\n')


def compressed(path, values):
    with gzip.open(path, 'wt') as handle:
        for value in values: handle.write(json.dumps(value) + '\n')


def table(path, rows):
    with path.open('w') as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)


def fingerprint(row): return hashlib.sha256(json.dumps(row, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True); args = parser.parse_args(); assert not args.output.exists()
    binary = '/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign'; rejected = []
    with tempfile.TemporaryDirectory(prefix='full-background-native-software-', dir='results') as directory:
        temp = Path(directory); source_root = temp / 'handoff'; source_root.mkdir(); sources = {}; all_rows = {}; seq = 'ACDEFGHIKLMNPQRSTVWY'
        for label, names in [('primary', 'AB'), ('reference', 'C'), ('background', 'CD'), ('legacy_reference', 'BE')]:
            folder = temp / label; folder.mkdir(); rows = []
            for name in names:
                xyz = [[3 * math.cos(i), 3 * math.sin(i), i] for i in range(len(seq))]
                if name == 'E': xyz = [[i, 0, 0] for i in range(len(seq))]
                confidence = [90 if name != 'D' or i < 2 else 50 for i in range(len(seq))]
                for mask in ['full', 'plddt70']:
                    blob, letters, positions = render_ca(dict(status='validated', sequence=seq, ca_xyz=xyz, ca_plddt=confidence), 70 if mask == 'plddt70' else None)
                    row = dict(model_id=name, version=1, mask=mask, status='ready' if len(positions) >= 3 else 'too_few_retained_residues',
                               source_sha256='synthetic-native-source-' + name, sequence=letters, original_positions=positions, retained_residues=len(positions), original_length=len(seq), coordinate_shard=label)
                    if row['status'] == 'ready':
                        p = folder / (name + '-' + mask + '.pdb'); p.write_bytes(blob); row.update(path=str(p), sha256=sha(p))
                    else: row['reason'] = 'fewer_than_three_retained_residues'
                    rows.append(row)
            manifest = folder / 'inputs.jsonl'; manifest.write_text(''.join(json.dumps(row) + '\n' for row in rows))
            sources[label] = dict(inputs=str(folder)); all_rows[label] = rows
        normalized = []; counts = Counter()
        for name in 'ABCDE':
            for mask in ['full', 'plddt70']:
                origins = []; originals = []
                for label, rows in all_rows.items():
                    for ordinal, row in enumerate(rows, 1):
                        if row['model_id'] == name and row['mask'] == mask:
                            origins.append(dict(collection=label, source_manifest=str(Path(sources[label]['inputs']) / 'inputs.jsonl'), source_row_number=ordinal, source_row_sha256=fingerprint(row)))
                            originals.append(row)
                semantic = {k: v for k, v in originals[0].items() if k not in ['path', 'coordinate_shard']}
                assert all({k: v for k, v in row.items() if k not in ['path', 'coordinate_shard']} == semantic for row in originals)
                selected = {k: v for k, v in originals[0].items() if k not in ['original_positions', 'coordinate_shard']}
                selected.update(semantic_sha256=fingerprint(semantic), collection_sources=origins); normalized.append(selected); counts[mask + ':' + selected['status']] += 1
        compressed(source_root / 'active_inputs.jsonl.gz', normalized)
        partition = []
        for a, b, old in [('A', 'B', True), ('A', 'C', False), ('C', 'D', False), ('A', 'E', False)]:
            ends = [(a, 1), (b, 1)]; key = hashlib.sha256(json.dumps(ends, separators=(',', ':')).encode()).hexdigest()
            partition.append(dict(pair_key=key, model_a=a, version_a=1, model_b=b, version_b=1, matching_old_sources=json.dumps(['background_old'] if old else []),
                                  measurement_disposition='pending_actual_input_result_and_numeric_reuse_checks' if old else 'native_measurement_required'))
        table(source_root / 'full_background_work_partition.tsv', partition)
        hp = temp / 'handoff-plan.json'; save(hp, dict(input_sources=sources)); summary = dict(active_models=5, active_model_mask_states=10, full_pairs=4, new_pairs=3, pending_catalog_reuse_pairs=1)
        rp, ap, archive, closure = [source_root / name for name in ['receipt.json', 'readback.json', 'archive.json', 'closed.json']]
        artifacts = {name: sha(source_root / name) for name in ['active_inputs.jsonl.gz', 'full_background_work_partition.tsv']}
        save(rp, dict(status='complete_full_expanded_background_native_input_handoff_pending_independent_readback', plan_sha256=sha(hp), **summary, artifacts=artifacts))
        save(ap, dict(status='passed_full_expanded_background_native_input_handoff_sql_readback', plan_sha256=sha(hp), producer_receipt_sha256=sha(rp), **summary))
        bindings = {str(hp): sha(hp), **{str(Path(spec['inputs']) / 'inputs.jsonl'): sha(Path(spec['inputs']) / 'inputs.jsonl') for spec in sources.values()}}
        for row in normalized:
            if row['status'] == 'ready': bindings[row['path']] = row['sha256']
        # Only closure provenance is synthetic here. This is not an actual process-journal proof.
        save(archive, dict(status='complete_verified_full_expanded_background_native_input_handoff_archive', source_hashes=bindings, services=['software_fixture_only'] * 2))
        save(closure, dict(status='complete_verified_full_expanded_background_native_input_handoff', **summary, exact_process_journals_checked=2,
                          full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive), bound_source_hashes=len(bindings), producer_receipt=str(rp), producer_receipt_sha256=sha(rp),
                          independent_readback=str(ap), independent_readback_sha256=sha(ap), active_input_counts=dict(counts)))
        np = temp / 'native-plan.json'; native_plan = dict(handoff_completion=str(closure), handoff_plan=str(hp), expected=summary, pins={binary: sha(binary)}, usalign=binary,
            output=str(temp / 'native'), options=['-mol', 'prot', '-mm', '0', '-outfmt', '0', '-ter', '2'], workers=2, per_pair_timeout_seconds=10,
            resources=dict(minimum_free_disk_gib=0), scope='Synthetic real-native software fixture; no biological pilot or actual closure proof.')
        save(np, native_plan)
        restored, new, _, bundle, _, _ = load_handoff(native_plan, True)
        assert restored['E', 1, 'full']['collection_sources'][0]['collection'] == 'legacy_reference' and restored['E', 1, 'full']['original_positions'] == list(range(1, 21))
        assert restored['D', 1, 'plddt70']['original_positions'] == [1, 2]
        nr = runner.run(np); assert nr['directed_dispositions'] == 12 and nr['counts'] == {'full:aligned': 6, 'plddt70:aligned': 4, 'plddt70:input_unavailable': 2}
        # An identical checkpoint replay must leave the saved bytes unchanged.
        first = new[0]; job = first['pair_key'], ('A', 1), ('C', 1), 'full', 0
        checkpoint = Path(native_plan['output']) / 'pairs' / job[0][:2] / (job[0] + '-full-0.json'); before = checkpoint.read_bytes()
        run_job(job, restored, native_plan, sha(np), bundle); assert checkpoint.read_bytes() == before
        dp = temp / 'assessment-plan.json'; save(dp, dict(native_plan=str(np), pins={}, output=str(temp / 'assessment'), scope=native_plan['scope']))
        dr = assessment.run(dp)
        qp = temp / 'reader-plan.json'; qp_config = dict(native_plan=str(np), assessment_plan=str(dp), pins={}, output=str(temp / 'reader.json'), scope=native_plan['scope']); save(qp, qp_config)
        qr = reader.run(qp); assert qr['directed_dispositions'] == 12 and qr['numerically_checked_alignments'] == 10
        assert sum(v for k, v in qr['geometry_counts'].items() if k.endswith(':degenerate_at_numeric_tolerance')) >= 4
        result_root = temp / 'assessment'; measurement = result_root / 'disposition_measurements.jsonl.gz'
        with gzip.open(measurement, 'rt') as handle: original_rows = [json.loads(line) for line in handle]
        original_receipt = json.loads((result_root / 'receipt.json').read_text())
        def aligned(rows): return next(row for row in rows if row['native_status'] == 'aligned')
        def unavailable(rows): return next(row for row in rows if row['native_status'] == 'input_unavailable')
        mutations = {
            'missing_directed_state': lambda rows, r: rows.pop(),
            'duplicated_directed_state': lambda rows, r: rows.append(copy.deepcopy(rows[0])),
            'false_model_role': lambda rows, r: aligned(rows).update(model_left='wrong'),
            'false_checkpoint_hash': lambda rows, r: aligned(rows).update(checkpoint_sha256='wrong'),
            'invented_unavailable_numeric': lambda rows, r: unavailable(rows).update(numerical={}),
            'promoted_unavailable': lambda rows, r: unavailable(rows).update(numerical_usable=True, numerical_exclusion_reasons=[]),
            'changed_native_status': lambda rows, r: aligned(rows).update(native_status='timeout'),
            'changed_numeric_identity': lambda rows, r: aligned(rows)['numerical'].update(sequence_identity_exact=.5),
            'changed_rmsd': lambda rows, r: aligned(rows)['numerical'].update(rmsd_recomputed=8.),
            'changed_rank': lambda rows, r: aligned(rows)['geometry'].update(rank_left=0),
            'removed_degenerate_flag': lambda rows, r: next(row for row in rows if 'nonunique_rotation' in row['numerical_exclusion_reasons']).update(numerical_usable=True, numerical_exclusion_reasons=[]),
            'false_total': lambda rows, r: r.update(directed_dispositions=11),
        }
        for name, mutate in mutations.items():
            rows = copy.deepcopy(original_rows); r = copy.deepcopy(original_receipt); mutate(rows, r); compressed(measurement, rows)
            r['artifacts'][measurement.name] = sha(measurement); save(result_root / 'receipt.json', r)
            bad = temp / (name + '-reader-plan.json'); save(bad, {**qp_config, 'output': str(temp / (name + '-readback.json'))})
            try: reader.run(bad)
            except (AssertionError, ValueError, KeyError): rejected.append(name)
            else: raise AssertionError('Rehashed false scientific export accepted: ' + name)
        # Synthetic native failure/timeout/parse and RMSD discrepancy remain explicit.
        native_root = Path(native_plan['output']); manifest = list(csv.DictReader((native_root / 'checkpoint_manifest.tsv').open(), delimiter='\t'))
        for item, status in zip(manifest[:4], ['native_error', 'parse_error', 'timeout', 'aligned']):
            path = native_root / item['path']; record = json.loads(path.read_text())
            if status == 'aligned':
                old_rmsd = record['metrics']['rmsd']; new_rmsd = old_rmsd + 1.
                text = record['stdout']; import re
                record['stdout'] = re.sub(r'(RMSD=\s*)[0-9.]+', lambda match: match[1] + str(new_rmsd), text, count=1); record['metrics']['rmsd'] = new_rmsd
            else:
                record.pop('metrics'); record['status'] = status
                if status == 'native_error': record['returncode'] = 9
                elif status == 'parse_error': record.update(returncode=0, stdout='synthetic unparseable text', error='ValueError: synthetic parsing failure')
                else: record.pop('returncode'); record['timeout_seconds'] = 10
            save(path, record); item.update(sha256=sha(path), status=record['status'])
        table(native_root / 'checkpoint_manifest.tsv', manifest)
        nr['counts'] = dict(Counter(Path(row['path']).stem.rsplit('-', 2)[1] + ':' + row['status'] for row in manifest)); nr['artifacts']['checkpoint_manifest.tsv'] = sha(native_root / 'checkpoint_manifest.tsv'); save(native_root / 'receipt.json', nr)
        failure_dp = temp / 'failure-assessment-plan.json'; save(failure_dp, dict(native_plan=str(np), pins={}, output=str(temp / 'failures'), scope=native_plan['scope'])); failure_dr = assessment.run(failure_dp)
        failure_qp = temp / 'failure-reader-plan.json'; save(failure_qp, {**qp_config, 'assessment_plan': str(failure_dp), 'output': str(temp / 'failure-readback.json')}); failure_qr = reader.run(failure_qp)
        assert failure_qr['directed_dispositions'] == 12
        assert any(k.endswith(':outside_printed_rounding') for k in failure_qr['rmsd_status_counts'])
        for status in ['native_error', 'parse_error', 'timeout']: assert any(k.endswith(':' + status) for k in failure_qr['counts'])
        # Source-row provenance restores positions independently of the compact export.
        chosen_manifest = Path(restored['E', 1, 'full']['collection_sources'][0]['source_manifest'])
        rows = [json.loads(line) for line in chosen_manifest.read_text().splitlines()]; rows[-2]['original_positions'].reverse(); chosen_manifest.write_text(''.join(json.dumps(row) + '\n' for row in rows))
        try: load_handoff(native_plan, True)
        except AssertionError: rejected.append('changed_original_residue_source')
        else: raise AssertionError('Changed original-position source accepted')
    result = dict(status='passed_expanded_background_native_numeric_quaternion_software_checks', directed_synthetic_native_states=12, numerically_checked_synthetic_alignments=10,
                  input_collections=4, pending_synthetic_catalog_pairs=1, unchanged_checkpoint_replay_passed=True, synthetic_native_failure_parse_timeout_and_rmsd_discrepancy_preserved=True,
                  rejected_rehashed_false_exports=rejected[:12], changed_original_position_source_rejected=True,
                  maximum_absolute_fixture_rmsd_difference=qr['maximum_absolute_rmsd_difference'], maximum_scaled_fixture_quaternion_curvature_error=qr['maximum_scaled_quaternion_curvature_error'],
                  script_hashes={str(p): sha(p) for p in [Path(__file__), Path(runner.__file__), Path(assessment.__file__), Path(reader.__file__), Path('scripts/expanded_background_native_handoff.py'), Path('scripts/expanded_background_measurement_sources.py')]},
                  usalign_sha256=sha(binary), scope='Five synthetic20-residue models/four collections including global overlaps, legacy-only model, short pLDDT70 and collinear fit. Actual USalign both masks/orders on3new pairs; pending old pair never computed. Compact-input original-position reconstruction, immutable checkpoint replay, full diagnostic/SVD and independent vectorized/quaternion baseline,12rehashed false exports and native failure/parse/timeout/RMSD discrepancy preservation. Only closure provenance is synthesized; no actual completion-journal audit, biological pilot, full source qualification or old-result acceptance claim.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__': main()
