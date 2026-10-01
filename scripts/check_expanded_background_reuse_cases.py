#!/usr/bin/env python3
"""Actual-native old-background I/O and full-grid reuse/corruption software checks."""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import math
import re
import tempfile
from collections import Counter
from pathlib import Path
import numpy as np
from assess_background_alignment_geometry import geometry
from background_alignment_handoff import load_handoff as background_handoff
from duplication_alignment_inputs import render_ca
from duplication_alignment_numeric_diagnostic import check_alignment
from duplication_alignment_numeric_readback import load_pdb
import expanded_background_reuse_sources as sources
import qualify_expanded_background_reuse as producer
import readback_expanded_background_reuse as reader
from run_duplication_alignments import run_job
from run_ortholog_pair_guide_comparison import sha


def save(path, value): path.write_text(json.dumps(value, indent=2) + '\n')
def lines(path, values): path.write_text(''.join(json.dumps(value) + '\n' for value in values))
def compressed(path, values):
    with gzip.open(path, 'wt') as handle:
        for value in values: handle.write(json.dumps(value) + '\n')
def table(path, rows):
    with path.open('w') as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
def pair(a, b, disposition='new_model_pair'):
    return dict(pair_key=hashlib.sha256(json.dumps(sorted([(a, 1), (b, 1)]), separators=(',', ':')).encode()).hexdigest(), model_a=a, version_a=1, model_b=b, version_b=1, work_disposition=disposition)


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True); args = parser.parse_args(); assert not args.output.exists()
    binary = '/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign'; rejected = []
    with tempfile.TemporaryDirectory(prefix='full-background-reuse-software-', dir='results') as directory:
        temp = Path(directory); seq = 'ACDEFGHIKLMNPQRSTVWY'; models = {}; original_inputs = {}; specs = {}
        for name in 'ABCD':
            raw = temp / (name + '.cif'); raw.write_text('Synthetic raw source contract fixture ' + name)
            models[name] = dict(model_id=name, version=1, path=str(raw), sha256=sha(raw), sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(), length=len(seq))
        for label, name in [('primary', 'A'), ('reference', 'B'), ('background', 'C')]:
            folder = temp / label; folder.mkdir(); inventory = temp / (label + '-models'); inventory.mkdir(); lines(inventory / 'models.jsonl', [models[name]])
            table(inventory / 'model_pairs.tsv', [pair('B', 'A')]); save(inventory / 'receipt.json', dict(artifacts={filename: sha(inventory / filename) for filename in ['models.jsonl', 'model_pairs.tsv']}))
            coord, audit = temp / (label + '-coords'), temp / (label + '-coordinate-reader'); coord.mkdir(); audit.mkdir(); save(coord / 'receipt.json', dict(status='synthetic_coordinate_contract'))
            coordinate_status = 'passed_background_exported_ca_readback' if label == 'background' else 'passed_duplication_exported_ca_readback'
            save(audit / 'receipt.json', dict(status=coordinate_status, producer_receipt_sha256=sha(coord / 'receipt.json')))
            ip = temp / (label + '-input-plan.json'); save(ip, dict(output=str(folder), coordinates=str(coord), readback=str(audit)))
            xyz = [[i, 0, 0] if name == 'B' else [3 * math.cos(i), 3 * math.sin(i), i] for i in range(len(seq))]
            confidence = [90 if name != 'C' or i < 2 else 50 for i in range(len(seq))]; rows = []
            for mask in ['full', 'plddt70']:
                blob, letters, positions = render_ca(dict(status='validated', sequence=seq, ca_xyz=xyz, ca_plddt=confidence), 70 if mask == 'plddt70' else None)
                row = dict(model_id=name, version=1, mask=mask, status='ready' if len(positions) >= 3 else 'too_few_retained_residues', source_sha256=models[name]['sha256'],
                           sequence=letters, original_positions=positions, retained_residues=len(positions), original_length=len(seq))
                if row['status'] == 'ready':
                    p = folder / (mask + '.pdb'); p.write_bytes(blob); row.update(path=str(p), sha256=sha(p))
                else: row['reason'] = 'fewer_than_three_retained_residues'
                rows.append(row); original_inputs[name, 1, mask] = row
            lines(folder / 'inputs.jsonl', rows)
            save(folder / 'receipt.json', dict(status='synthetic_input_contract', plan_sha256=sha(ip), inventory_receipt_sha256=sha(inventory / 'receipt.json'),
                coordinate_receipt_sha256=sha(coord / 'receipt.json'), readback_receipt_sha256=sha(audit / 'receipt.json'), models=1, input_dispositions=2,
                counts=dict(Counter(row['mask'] + ':' + row['status'] for row in rows)), artifacts={'inputs.jsonl': sha(folder / 'inputs.jsonl')}))
            specs[label] = dict(inputs=str(folder), input_plan=str(ip), model_inventory=str(inventory), model_file='models.jsonl', status='synthetic_input_contract', inventory_receipt_field='inventory_receipt_sha256', coordinate_readback_status=coordinate_status)
        reference_inventory = temp / 'reference-inventory'; reference_inventory.mkdir(); lines(reference_inventory / 'models.jsonl', [models[name] for name in 'AB'])
        table(reference_inventory / 'model_pairs.tsv', [pair('B', 'A')]); save(reference_inventory / 'receipt.json', dict(artifacts={filename: sha(reference_inventory / filename) for filename in ['models.jsonl', 'model_pairs.tsv']}))
        background_inventory = temp / 'background-inventory'; background_inventory.mkdir(); lines(background_inventory / 'active_models.jsonl', [models[name] for name in 'ABC'])
        bg_pairs = [pair('A', 'B', 'already_in_reference_queue'), pair('A', 'C'), pair('B', 'C')]; table(background_inventory / 'model_pairs.tsv', bg_pairs)
        save(background_inventory / 'receipt.json', dict(status='complete_background_measurement_inventory_pending_readback', distinct_eligible_model_pairs=3, new_model_pairs=2,
            artifacts={filename: sha(background_inventory / filename) for filename in ['active_models.jsonl', 'model_pairs.tsv']}))
        bg_proof = temp / 'background-inventory-proof.json'; save(bg_proof, dict(status='passed_full_background_measurement_inventory_readback', producer_receipt_sha256=sha(background_inventory / 'receipt.json'), active_models=3))
        plans = {}
        for label, inventory in [('reference_old', reference_inventory), ('background_old', background_inventory)]:
            p = dict(inventory=str(inventory), output=str(temp / (label + '-native')), usalign=binary, options=['-mol', 'prot', '-mm', '0', '-outfmt', '0', '-ter', '2'], per_pair_timeout_seconds=10, pins={binary: sha(binary)})
            if label == 'background_old': p.update(input_sources=specs, inventory_readback=str(bg_proof), existing_queues=[dict(path=str(reference_inventory), label='already_in_reference_queue')], resources=dict(new_model_pairs=2))
            pp = temp / (label + '-native-plan.json'); save(pp, p); plans[label] = pp
        def measure(label):
            pp = plans[label]; p = json.loads(pp.read_text()); root = Path(p['output']); root.mkdir()
            if label == 'background_old': _, pairs, mb, bundle = background_handoff(p)
            else:
                pairs = [pair('B', 'A')]; mb = {str(Path(specs[name]['inputs']) / filename): sha(Path(specs[name]['inputs']) / filename) for name in ['primary', 'reference'] for filename in ['receipt.json', 'inputs.jsonl']}
                bundle = hashlib.sha256(json.dumps(mb, sort_keys=True).encode()).hexdigest()
            checkpoints = []; numeric_rows = []; geometry_rows = []
            for row in pairs:
                ends = [(row['model_a'], 1), (row['model_b'], 1)]
                for mask in ['full', 'plddt70']:
                    for order in [0, 1]:
                        desired = ends if order == 0 else ends[::-1]; job = row['pair_key'], *desired, mask, order
                        path, status = run_job(job, original_inputs, p, sha(pp), bundle)
                        record = json.loads(path.read_text())
                        if label == 'background_old' and row['model_a'] == 'A' and mask == 'full' and order == 0:
                            record['stdout'] = re.sub(r'(RMSD=\s*)[0-9.]+', lambda m: m[1] + str(record['metrics']['rmsd'] + 1.), record['stdout'], count=1)
                            record['metrics']['rmsd'] += 1.; save(path, record)
                        checkpoints.append(dict(path=str(path.relative_to(root)), sha256=sha(path), status=status))
                        if status != 'aligned': continue
                        coords = [load_pdb(original_inputs[(*end, mask)]) for end in desired]; numeric = check_alignment(record, *coords)
                        numeric_rows.append(dict(pair_key=row['pair_key'], mask=mask, order=order, **numeric))
                        texts = [record['metrics'][name] for name in ['alignment_left', 'alignment_right']]; present = [np.array(list(text)) != '-' for text in texts]; common = present[0] & present[1]
                        indexes = [(np.cumsum(v) - 1)[common] for v in present]; g = geometry(*[c[1][ix] for c, ix in zip(coords, indexes)])
                        geometry_rows.append(dict(pair_key=row['pair_key'], mask=mask, order=order, rmsd_status=numeric['rmsd_status'], **g))
            table(root / 'checkpoint_manifest.tsv', checkpoints)
            save(root / 'receipt.json', dict(status='complete_background_alignment_dispositions_pending_readback' if label == 'background_old' else 'complete_reference_alignment_dispositions_pending_readback',
                 plan_sha256=sha(pp), input_bindings=mb, input_bundle_sha256=bundle, distinct_model_pairs=len(pairs), directed_dispositions=len(checkpoints),
                 counts=dict(Counter(Path(row['path']).stem.rsplit('-', 2)[1] + ':' + row['status'] for row in checkpoints)), artifacts={'checkpoint_manifest.tsv': sha(root / 'checkpoint_manifest.tsv')}))
            dp, gp, qp = [temp / (label + '-' + name + '-plan.json') for name in ['diagnostic', 'geometry', 'reader']]
            diagnostic, geometric = temp / (label + '-diagnostic'), temp / (label + '-geometry'); diagnostic.mkdir(); geometric.mkdir()
            save(dp, dict(source_plan=str(pp), output=str(diagnostic))); save(gp, dict(diagnostic_plan=str(dp), output=str(geometric)))
            qr_path = temp / (label + '-geometry-readback.json'); save(qp, dict(source_plan=str(gp), output=str(qr_path)))
            table(diagnostic / 'numeric_readback.tsv', numeric_rows); table(geometric / 'alignment_geometry.tsv', geometry_rows)
            statuses = ('complete_background_alignment_rmsd_diagnostic_not_scientific_acceptance', 'complete_background_alignment_geometry_pending_independent_readback', 'passed_full_background_geometry_readback') if label == 'background_old' else ('complete_reference_alignment_rmsd_diagnostic_not_scientific_acceptance', 'complete_reference_alignment_geometry_pending_independent_readback', 'passed_full_reference_geometry_readback')
            save(diagnostic / 'receipt.json', dict(status=statuses[0], plan_sha256=sha(dp), producer_receipt_sha256=sha(root / 'receipt.json'), directed_dispositions=len(checkpoints),
                 numerically_checked_alignments=len(numeric_rows), counts=json.loads((root / 'receipt.json').read_text())['counts'], artifacts={'numeric_readback.tsv': sha(diagnostic / 'numeric_readback.tsv')}))
            gc = dict(Counter(row['mask'] + ':' + row['geometry_status'] for row in geometry_rows))
            save(geometric / 'receipt.json', dict(status=statuses[1], plan_sha256=sha(gp), diagnostic_receipt_sha256=sha(diagnostic / 'receipt.json'), alignments=len(geometry_rows), counts=gc, artifacts={'alignment_geometry.tsv': sha(geometric / 'alignment_geometry.tsv')}))
            save(qr_path, dict(status=statuses[2], plan_sha256=sha(qp), producer_receipt_sha256=sha(geometric / 'receipt.json'), alignments_checked=len(geometry_rows), counts=gc))
            return dict(kind='background' if label == 'background_old' else 'reference', alignment_plan=str(pp), alignment_status=json.loads((root / 'receipt.json').read_text())['status'],
                input_sources=[specs[name] for name in ['primary', 'reference']], diagnostic=str(diagnostic), diagnostic_status=statuses[0], geometry=str(geometric), geometry_status=statuses[1],
                geometry_readback=str(qr_path), geometry_readback_status=statuses[2], diagnostic_plan=str(dp), geometry_plan=str(gp), geometry_readback_plan=str(qp), launches=[])
        old = {label: measure(label) for label in plans}
        current = copy.deepcopy(original_inputs)
        for mask in ['full', 'plddt70']:
            current['D', 1, mask] = {**copy.deepcopy(current['A', 1, mask]), 'model_id': 'D', 'source_sha256': models['D']['sha256']}
        current_root = temp / 'current'; current_root.mkdir(); compressed(current_root / 'active_models.jsonl.gz', list(models.values())); partition = []
        for a, b, matches in [('A', 'B', ['reference_old', 'background_old']), ('A', 'C', ['background_old']), ('B', 'C', ['background_old']), ('A', 'D', [])]:
            partition.append({**pair(a, b), 'matching_old_sources': json.dumps(matches), 'measurement_disposition': 'pending_actual_input_result_and_numeric_reuse_checks' if matches else 'native_measurement_required'})
        table(current_root / 'full_background_work_partition.tsv', partition)
        target = temp / 'target-native-plan.json'; save(target, dict(usalign=binary, options=['-mol', 'prot', '-mm', '0', '-outfmt', '0', '-ter', '2'], per_pair_timeout_seconds=10))
        def current_handoff(plan, include_positions=False):
            bindings = {str(current_root / name): sha(current_root / name) for name in ['active_models.jsonl.gz', 'full_background_work_partition.tsv']}
            return current, [row for row in partition if not json.loads(row['matching_old_sources'])], bindings, 'synthetic_bundle', partition, current_root
        sources.load_handoff = current_handoff
        config = dict(target_alignment_plan=str(target), sources=old, expected=dict(active_models=4, full_pairs=4, new_pairs=1, selected_source_pairs=dict(reference_old=1, background_old=2)), pins={},
                      output=str(temp / 'reuse'), scope='Synthetic current proof I/O and original source proof contracts; real native PDB/checkpoint/software geometry checks, not a biological pilot or original journal proof.')
        pp = temp / 'reuse-plan.json'; save(pp, config); result = producer.run(pp); audit = reader.run(pp, temp / 'readback.json')
        assert result['directed_dispositions'] == 16 and result['counts'] == {'verified_identical_input_checkpoint_and_retained_disposition': 12, 'new_native_measurement_pending': 4}
        assert audit['numerically_reconstructed_reuse_alignments'] == 8 and audit['numerical_counts'] == {'nonunique_rotation': 6, 'rmsd_discrepancy': 1, 'input_unavailable': 4, 'usable': 1}
        output = Path(config['output']); path = output / 'background_reuse_dispositions.jsonl.gz'
        with gzip.open(path, 'rt') as handle: baseline = [json.loads(line) for line in handle]
        original_receipt = json.loads((output / 'receipt.json').read_text())
        ab = [row for row in baseline if row['model_a'] == 'A' and row['model_b'] == 'B']; assert all(row['selected_source'] == 'reference_old' and row['source_order'] == 1 - row['order'] for row in ab)
        def reused(rows): return next(row for row in rows if row['source_checkpoint'] is not None)
        mutations = {
            'missing_state': lambda rows, r: rows.pop(), 'duplicated_state': lambda rows, r: rows.insert(0, copy.deepcopy(rows[0])),
            'wrong_preference': lambda rows, r: next(row for row in rows if row['selected_source'] == 'reference_old').update(selected_source='background_old'),
            'wrong_original_order': lambda rows, r: reused(rows).update(source_order=7),
            'wrong_model_role': lambda rows, r: reused(rows).update(model_a='wrong'),
            'invented_native_hash': lambda rows, r: reused(rows).update(source_checkpoint_sha256='wrong'),
            'cleared_degenerate': lambda rows, r: next(row for row in rows if 'nonunique_rotation' in row['numerical_exclusion_reasons']).update(numerical_usable=True, numerical_exclusion_reasons=[]),
            'cleared_discrepancy': lambda rows, r: next(row for row in rows if 'rmsd_discrepancy' in row['numerical_exclusion_reasons']).update(numerical_usable=True, numerical_exclusion_reasons=[]),
            'promoted_new_pair': lambda rows, r: next(row for row in rows if row['selected_source'] == '').update(reuse_status='verified_identical_input_checkpoint_and_retained_disposition'),
            'changed_numeric': lambda rows, r: next(row for row in rows if row['source_numeric'] is not None)['source_numeric'].update(sequence_identity_exact='0.5'),
            'changed_geometry': lambda rows, r: next(row for row in rows if row['source_geometry'] is not None)['source_geometry'].update(rank_left='0'),
            'false_total': lambda rows, r: r.update(directed_dispositions=15),
        }
        for name, mutate in mutations.items():
            rows = copy.deepcopy(baseline); r = copy.deepcopy(original_receipt); mutate(rows, r); compressed(path, rows); r['artifacts'][path.name] = sha(path); save(output / 'receipt.json', r)
            try: reader.run(pp, temp / (name + '-readback.json'))
            except (AssertionError, ValueError, KeyError, StopIteration): rejected.append(name)
            else: raise AssertionError('Rehashed false reuse export accepted: ' + name)
        # A current mask incompatibility stays pending without switching old source.
        current['B', 1, 'full']['reason'] = 'synthetic_declared_input_difference'
        ip = temp / 'incompatible-plan.json'; save(ip, {**config, 'output': str(temp / 'incompatible')}); ir = producer.run(ip); ia = reader.run(ip, temp / 'incompatible-readback.json')
        assert ia['counts'] == {'incompatible_inputs_require_new_measurement': 4, 'verified_identical_input_checkpoint_and_retained_disposition': 8, 'new_native_measurement_pending': 4}
        assert ia['numerically_reconstructed_reuse_alignments'] == 4
    result = dict(status='passed_full_expanded_background_actual_native_reuse_software_checks', synthetic_full_pairs=4, synthetic_directed_states=16, synthetic_reused_states=12,
        synthetic_reconstructed_alignments=8, old_background_three_partition_IO_checked=True, reversed_old_reference_orders_preserved=True, numerical_flags_and_unavailable_states_preserved=True,
        incompatible_mask_requires_new_measurement_without_source_substitution=True, rejected_rehashed_false_exports=rejected,
        maximum_absolute_fixture_reused_rmsd_difference=audit['maximum_absolute_reused_rmsd_difference'], maximum_scaled_fixture_reused_quaternion_curvature_error=audit['maximum_scaled_reused_quaternion_curvature_error'],
        script_hashes={str(p): sha(p) for p in [Path(__file__), Path(producer.__file__), Path(reader.__file__), Path(sources.__file__)]}, usalign_sha256=sha(binary),
        scope='Complete16-state synthetic grid with actual native original reference/background records and full three-partition old-background input/proof/source contracts. Current closed input-union I/O and original proof/journal provenance synthesized only here. Full actual-source input compatibility, fixed reference-before-background choice, reversed old orders, original flag retention, SVD/quaternion reconstruction,12rehashed false outputs and incompatible masks checked. No real full input-union/old-journal proof, accepted production reuse, biological pilot or biological effect claim.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__': main()
