#!/usr/bin/env python3
"""Complete synthetic native reuse handoff with missing, incompatible and excluded cases."""
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
from duplication_alignment_inputs import render_ca
from duplication_alignment_numeric_readback import load_pdb
from duplication_alignment_numeric_diagnostic import check_alignment
from assess_reference_alignment_geometry import geometry
from qualify_reference_alignment_reuse import numeric_exclusions
from run_duplication_alignments import run_job
from screen_duplication_alignment_reuse import sha


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def pair(a, b):
    return dict(pair_key=hashlib.sha256(json.dumps(sorted([a, b]), separators=(',', ':')).encode()).hexdigest(),
                model_a=a[0], version_a=a[1], model_b=b[0], version_b=b[1])


def main():
    assert numeric_exclusions(dict(status='aligned'), dict(rmsd_status='outside_printed_rounding', aligned_length='2'),
                              dict(geometry_status='degenerate_at_numeric_tolerance')) == ['rmsd_discrepancy', 'fewer_than_three_pairs', 'nonunique_rotation']
    assert numeric_exclusions(dict(status='timeout'), None, None) == ['timeout']
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        versions = dict(A=6, B=10, C=6, D=10, E=6, F=10, G=6)
        keys = {n: (n, v) for n, v in versions.items()}
        sequence = 'ACDEFGHIKLMNPQRSTVWY'
        models = {}
        coordinate_records = {}
        for n, key in keys.items():
            raw = root / (n + '.fixture-raw')
            raw.write_text('synthetic raw source binding ' + n)
            models[key] = dict(model_id=n, version=key[1], path=str(raw), sha256=sha(raw),
                               sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(), length=len(sequence))
            xyz = [[float(i), 0., 0.] for i in range(len(sequence))] if n == 'C' else [[3 * math.cos(i), 3 * math.sin(i), i] for i in range(len(sequence))]
            coordinate_records[key] = dict(status='validated', sequence=sequence, ca_xyz=xyz, ca_plddt=[90.] * len(sequence))
        def catalog(name, names, pairs, status):
            folder = root / name
            folder.mkdir()
            (folder / 'models.jsonl').write_text(''.join(json.dumps(models[keys[n]]) + '\n' for n in names))
            table(folder / 'model_pairs.tsv', pairs)
            write(folder / 'receipt.json', dict(status=status, unique_distinct_model_pairs=len(pairs),
                                                artifacts={p.name: sha(p) for p in folder.iterdir()}))
            return folder
        fullpairs = [pair(keys['A'], keys['B']), pair(keys['A'], keys['C']), pair(keys['D'], keys['E']), pair(keys['F'], keys['G'])]
        inventory = catalog('current-catalog', 'ABCDEFG', fullpairs, 'complete_provisional_reference_comparison_inventory')
        proof = root / 'inventory-proof.json'
        write(proof, dict(status='passed_full_reference_comparison_ledger_and_native_model_readback', producer_receipt_sha256=sha(inventory / 'receipt.json')))
        def materialize(name, names, current=False):
            folder = root / name
            folder.mkdir()
            rows = []
            for n in names:
                key = keys[n]
                for mask in ['full', 'plddt70']:
                    rec = coordinate_records[key].copy()
                    if n == 'C' and mask == 'plddt70':
                        rec['ca_plddt'] = [90., 90.] + [0.] * (len(sequence) - 2)
                    if current and n == 'F' and mask == 'plddt70':
                        rec['ca_plddt'] = [0.] + [90.] * (len(sequence) - 1)
                    blob, letters, positions = render_ca(rec, threshold=70 if mask == 'plddt70' else None)
                    row = dict(model_id=n, version=key[1], mask=mask, source_sha256=models[key]['sha256'],
                               original_positions=positions, sequence=letters, retained_residues=len(positions),
                               original_length=len(sequence))
                    if len(positions) < 3:
                        row['status'] = 'too_few_retained_residues'
                    else:
                        path = folder / (n + '-' + mask + '.pdb')
                        path.write_bytes(blob)
                        row.update(status='ready', path=str(path), sha256=sha(path))
                    rows.append(row)
            (folder / 'inputs.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in rows))
            coord = root / (name + '-coordinate')
            reader = root / (name + '-coordinate-reader')
            coord.mkdir(); reader.mkdir()
            write(coord / 'receipt.json', dict(status='fixture_coordinate_producer', models=len(names)))
            write(reader / 'receipt.json', dict(status='passed_duplication_exported_ca_readback', producer_receipt_sha256=sha(coord / 'receipt.json')))
            ip = root / (name + '-plan.json')
            write(ip, dict(output=str(folder), coordinates=str(coord), readback=str(reader)))
            write(folder / 'receipt.json', dict(status='complete_duplication_alignment_input_materialization',
                 plan_sha256=sha(ip), coordinate_receipt_sha256=sha(coord / 'receipt.json'),
                 readback_receipt_sha256=sha(reader / 'receipt.json'), models=len(names), input_dispositions=len(rows),
                 counts=dict(Counter(r['mask'] + ':' + r['status'] for r in rows)), artifacts={'inputs.jsonl': sha(folder / 'inputs.jsonl')}))
            spec = dict(inputs=str(folder), input_plan=str(ip), status='complete_duplication_alignment_input_materialization')
            if current:
                proof = root / (name + '-proof.json')
                write(proof, dict(status='passed_full_duplication_alignment_input_readback', source_receipt_sha256=sha(folder / 'receipt.json')))
                spec.update(readback=str(proof), readback_status='passed_full_duplication_alignment_input_readback', readback_receipt_field='source_receipt_sha256')
            return spec, {(r['model_id'], r['version'], r['mask']): r for r in rows}
        current_spec, current_inputs = materialize('current-inputs', 'ABCDEFG', current=True)
        oldprimary, op_inputs = materialize('old-primary-inputs', 'AB')
        oldextra, oe_inputs = materialize('old-additional-inputs', 'CFG')
        allold = {**op_inputs, **oe_inputs}
        usalign = '/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign'
        target = root / 'target-plan.json'
        write(target, dict(usalign=usalign, options=['-mol', 'prot', '-mm', '0', '-outfmt', '0', '-ter', '2']))
        native_specs = {}
        def completed_source(label, kind, names, pairs, inputspecs):
            q = catalog(label + '-catalog', names, pairs, 'fixture_old_catalog')
            folder = root / (label + '-native')
            folder.mkdir()
            sp = root / (label + '-alignment-plan.json')
            config = dict(output=str(folder), usalign=usalign, options=json.loads(target.read_text())['options'],
                          pins={usalign: sha(usalign)}, per_pair_timeout_seconds=10)
            config['queue' if kind == 'primary' else 'inventory'] = str(q)
            write(sp, config)
            mb = {str(Path(s['inputs']) / name): sha(Path(s['inputs']) / name) for s in inputspecs for name in ['receipt.json', 'inputs.jsonl']}
            bundle = sha(Path(inputspecs[0]['inputs']) / 'inputs.jsonl') if kind == 'primary' else hashlib.sha256(json.dumps(mb, sort_keys=True).encode()).hexdigest()
            checkpoint_rows = []
            numeric_rows, geometry_rows = [], []
            for r in pairs:
                ends = [(r['model_a'], r['version_a']), (r['model_b'], r['version_b'])]
                for mask in ['full', 'plddt70']:
                    for order in [0, 1]:
                        directed = ends if order == 0 else ends[::-1]
                        path, state = run_job((r['pair_key'], *directed, mask, order), allold, config, sha(sp), bundle)
                        checkpoint_rows.append(dict(path=str(path.relative_to(folder)), sha256=sha(path), status=state))
                        if state != 'aligned':
                            continue
                        record = json.loads(path.read_text())
                        coords = [load_pdb(allold[(*e, mask)]) for e in directed]
                        numeric = check_alignment(record, *coords)
                        numeric_rows.append(dict(pair_key=r['pair_key'], mask=mask, order=order, **numeric))
                        strings = [record['metrics'][f'alignment_{s}'] for s in ['left', 'right']]
                        nongap = [np.array(list(s)) != '-' for s in strings]
                        paired = nongap[0] & nongap[1]
                        indices = [(np.cumsum(x) - 1)[paired] for x in nongap]
                        geometry_rows.append(dict(pair_key=r['pair_key'], mask=mask, order=order,
                              rmsd_status=numeric['rmsd_status'], **geometry(*[c[1][i] for c, i in zip(coords, indices)])))
            table(folder / 'checkpoint_manifest.tsv', checkpoint_rows)
            status = 'complete_duplication_alignment_dispositions_pending_readback' if kind == 'primary' else 'complete_reference_alignment_dispositions_pending_readback'
            receipt = dict(status=status, plan_sha256=sha(sp), directed_dispositions=len(checkpoint_rows), artifacts={'checkpoint_manifest.tsv': sha(folder / 'checkpoint_manifest.tsv')})
            if kind == 'primary':
                receipt.update(input_manifest_sha256=bundle, input_receipt_sha256=mb[str(Path(inputspecs[0]['inputs']) / 'receipt.json')])
            else:
                receipt.update(input_bundle_sha256=bundle, input_bindings=mb)
            write(folder / 'receipt.json', receipt)
            diag, geo = root / (label + '-numeric'), root / (label + '-geometry')
            diag.mkdir(); geo.mkdir()
            table(diag / 'numeric_readback.tsv', numeric_rows)
            table(geo / 'alignment_geometry.tsv', geometry_rows)
            write(diag / 'receipt.json', dict(status='fixture_complete_numeric', producer_receipt_sha256=sha(folder / 'receipt.json'), numerically_checked_alignments=len(numeric_rows), artifacts={'numeric_readback.tsv': sha(diag / 'numeric_readback.tsv')}))
            write(geo / 'receipt.json', dict(status='fixture_complete_geometry', diagnostic_receipt_sha256=sha(diag / 'receipt.json'), alignments=len(geometry_rows), artifacts={'alignment_geometry.tsv': sha(geo / 'alignment_geometry.tsv')}))
            proof = root / (label + '-geometry-proof.json')
            write(proof, dict(status='fixture_full_geometry_readback', producer_receipt_sha256=sha(geo / 'receipt.json'), alignments_checked=len(geometry_rows)))
            native_specs[label] = dict(kind=kind, alignment_plan=str(sp), alignment_status=status, input_sources=inputspecs,
                                       diagnostic=str(diag), diagnostic_status='fixture_complete_numeric', geometry=str(geo),
                                       geometry_status='fixture_complete_geometry', geometry_readback=str(proof), geometry_readback_status='fixture_full_geometry_readback')
        # The original primary direction is reversed relative to the current pair.
        completed_source('primary_expanded_completed', 'primary', 'AB', [pair(keys['B'], keys['A'])], [oldprimary])
        completed_source('reference_old', 'reference', 'ABCFG', [fullpairs[1], fullpairs[3]], [oldprimary, oldextra])
        union = root / 'union'
        union.mkdir()
        choices = [['primary_expanded_completed', 'reference_old'], ['reference_old'], [], ['reference_old']]
        candidate_rows = [{**r, 'new_sources': json.dumps(['reference_expanded']), 'matching_old_sources': json.dumps(s),
                           'changed_old_sources': '[]', 'disposition': 'matching_catalog_sources_pending_input_and_result_checks' if s else 'new_pair_requires_alignment'} for r, s in zip(fullpairs, choices)]
        table(union / 'pair_reuse_candidates.tsv', candidate_rows)
        write(union / 'receipt.json', dict(status='complete_full_pair_union_catalog_reuse_screen_not_authorization',
               source_hashes={str(inventory / n): sha(inventory / n) for n in ['receipt.json', 'models.jsonl', 'model_pairs.tsv']},
               per_new_source_counts={'reference_expanded:matching_catalog_sources_pending_input_and_result_checks': 3},
               artifacts={'pair_reuse_candidates.tsv': sha(union / 'pair_reuse_candidates.tsv')}))
        unionproof = root / 'union-proof.json'
        write(unionproof, dict(status='passed_full_pair_union_reuse_candidate_sql_readback', producer_receipt_sha256=sha(union / 'receipt.json')))
        pp = root / 'plan.json'
        out = root / 'qualified'
        write(pp, dict(inventory=str(inventory), inventory_readback=str(proof), reuse_candidates=str(union),
                       reuse_readback=str(unionproof), current_inputs=[current_spec], sources=native_specs,
                       target_alignment_plan=str(target), output=str(out), pins={}))
        subprocess.run([sys.executable, 'scripts/qualify_reference_alignment_reuse.py', '--plan', str(pp)], check=True, capture_output=True)
        reader = [sys.executable, 'scripts/readback_reference_alignment_reuse.py', '--plan', str(pp)]
        subprocess.run([*reader, '--output', str(root / 'passed.json')], check=True, capture_output=True)
        rows = [json.loads(line) for line in (out / 'reference_reuse_dispositions.jsonl').read_text().splitlines()]
        assert len(rows) == 16
        assert Counter(r['reuse_status'] for r in rows) == {'verified_identical_input_checkpoint_and_retained_disposition': 10, 'new_native_measurement_pending': 4, 'incompatible_inputs_require_new_measurement': 2}
        assert any(r['source_order'] == 1 and r['order'] == 0 for r in rows if r['pair_key'] == fullpairs[0]['pair_key'])
        excluded = next(i for i, r in enumerate(rows) if r['source_native_status'] == 'aligned' and not r['numerical_usable'])
        assert 'nonunique_rotation' in rows[excluded]['numerical_exclusion_reasons']
        assert any(r['source_native_status'] == 'input_unavailable' for r in rows)
        path = out / 'reference_reuse_dispositions.jsonl'
        initial = path.read_bytes()
        rp = out / 'receipt.json'
        original_receipt = rp.read_bytes()
        for i, (index, field, value) in enumerate([(excluded, 'numerical_usable', True),
                                                  (excluded, 'numerical_exclusion_reasons', []),
                                                  (0, 'selected_source', 'wrong'), (0, 'order', 7),
                                                  (0, 'source_checkpoint_sha256', 'wrong')]):
            modified = [json.loads(line) for line in initial.decode().splitlines()]
            modified[index][field] = value
            path.write_text(''.join(json.dumps(r) + '\n' for r in modified))
            receipt = json.loads(original_receipt)
            receipt['artifacts'][path.name] = sha(path)
            write(rp, receipt)
            result = subprocess.run([*reader, '--output', str(root / f'false-{i}.json')], capture_output=True)
            assert result.returncode != 0 and not (root / f'false-{i}.json').exists()
        path.write_bytes(initial); rp.write_bytes(original_receipt)
        path = out / 'model_input_identity_checks.jsonl'
        items = [json.loads(line) for line in path.read_text().splitlines()]
        items[0]['current_projection_sha256'] = 'wrong'
        path.write_text(''.join(json.dumps(r) + '\n' for r in items))
        receipt = json.loads(original_receipt); receipt['artifacts'][path.name] = sha(path); write(rp, receipt)
        result = subprocess.run([*reader, '--output', str(root / 'false-input.json')], capture_output=True)
        assert result.returncode != 0 and not (root / 'false-input.json').exists()
    print('Passed complete16-state native synthetic handoff:versions6/10, directed reversal, both masks, missing and incompatible inputs, retained degenerate/short states;6 rehashed false exports rejected.')


if __name__ == '__main__':
    main()
