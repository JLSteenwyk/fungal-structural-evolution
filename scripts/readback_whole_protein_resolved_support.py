#!/usr/bin/env python3
"""Verify every support link using saved NPZ matrices and source certificates."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from check_joint_support_certificate import check_certificate
from run_after_verified_dependencies import live, terminal_state
from screen_duplication_alignment_reuse import sha


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.plan.read_text())
    launch = json.loads(Path(config['launch']).read_text())
    while True:
        if live(launch):
            time.sleep(30); continue
        state = terminal_state(launch['unit'])
        if state['ActiveState'] in ['active', 'activating', 'deactivating', 'reloading']:
            time.sleep(30); continue
        assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0'), state
        break
    assert sha(launch['plan']) == launch['plan_sha256']
    plan = json.loads(Path(launch['plan']).read_text())
    root = Path(plan['output']); rp = root/'receipt.json'
    receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_resolved_whole_protein_support_pending_readback'
    assert receipt['plan_sha256'] == sha(launch['plan'])
    bindings = {str(args.plan): sha(args.plan), **config['pins'], **receipt['source_hashes'], str(rp): sha(rp)}
    bindings.update({str(root/name): h for name, h in receipt['artifacts'].items()})
    def verify():
        for path, expected in bindings.items():
            assert sha(path) == expected, path
    verify()
    def indexed(path, field):
        records = {}
        for line in path.open():
            row = json.loads(line)
            assert row[field] not in records
            records[row[field]] = row
        return records
    inputs = Path(plan['inputs']); support = Path(plan['support']); witnesses = Path(plan['witnesses'])
    recipes = indexed(inputs/'input_manifest.jsonl', 'fit_input_id')
    original_map = indexed(support/'input_support_map.jsonl', 'fit_input_id')
    certificates = indexed(support/'joint_support.jsonl', 'geometry_id')
    repairs = indexed(witnesses/'witnesses.jsonl', 'geometry_id')
    unresolved = {k for k, row in certificates.items() if row['result']['classification'].startswith('unresolved')}
    assert set(repairs) == unresolved and len(repairs) == 6
    assert set(recipes) == set(original_map) and len(recipes) == 75070
    proof = json.loads(Path(plan['support_proof']).read_text())
    assert proof['status'] == 'passed_full_whole_protein_joint_support_readback'
    assert proof['source_receipt_sha256'] == sha(support/'receipt.json')
    repaired_proof = json.loads(Path(plan['witness_proof']).read_text())
    assert repaired_proof['status'] == 'completed_serialized_whole_protein_support_witness_checks'
    assert repaired_proof['receipt_sha256'] == sha(witnesses/'receipt.json')
    checked = set(); seen = set(); input_counts = Counter(); changed_inputs = 0; certificate_checks = 0
    verified_map = {}
    for line in (root/'resolved_input_support.jsonl').open():
        row = json.loads(line); identifier = row['fit_input_id']
        assert identifier not in seen and identifier in recipes; seen.add(identifier)
        recipe = recipes[identifier]; spec = recipe['recipe']['specification']; path = inputs/recipe['path']
        assert sha(path) == recipe['sha256'] == row['input_sha256']
        assert digest(spec) == identifier
        with np.load(path, allow_pickle=False) as a:
            matrix = a['matrix']; ordered = a['row_identity']
        assert matrix.dtype == np.dtype('<f8') and np.isfinite(matrix).all()
        assert matrix.shape == (spec['records'], len(spec['columns'])) and np.all(matrix[:, 1] == 1)
        assert hashlib.sha256(matrix.tobytes()).hexdigest() == spec['values_sha256']
        identity = hashlib.sha256(ordered.tobytes()).hexdigest()
        assert identity == spec['ordered_identity_sha256']
        x = np.ascontiguousarray(matrix[:, 2:], dtype='<f8'); x[x == 0] = 0.
        geometry_spec = dict(records=len(x), columns=spec['columns'][2:],
                             values_sha256=hashlib.sha256(x.tobytes()).hexdigest(), ordered_identity_sha256=identity)
        geometry_id = digest(geometry_spec)
        assert geometry_id == row['geometry_id'] == original_map[identifier]['geometry_id']
        original = certificates[geometry_id]
        assert original['specification'] == geometry_spec
        assert original_map[identifier]['classification'] == original['result']['classification']
        if geometry_id in repairs:
            repair = repairs[geometry_id]
            assert repair['accepted'] and repair['original'] == original['result']
            assert repair['specification'] == geometry_spec and repair['representative_input'] == original['representative_input']
            certificate = repair['candidate']; kind = 'verified_repaired_witness'; certificate_file = witnesses/'witnesses.jsonl'
        else:
            certificate = original['result']; kind = 'original_verified_certificate'; certificate_file = support/'joint_support.jsonl'
        expected = dict(fit_input_id=identifier, input_sha256=recipe['sha256'], geometry_id=geometry_id,
                        original_classification=original['result']['classification'], classification=certificate['classification'],
                        source_kind=kind, certificate_file=str(certificate_file), certificate_file_sha256=bindings[str(certificate_file)])
        assert row == expected
        assert certificate['classification'] == 'zero_supported_to_numeric_tolerance' and certificate['certificate_valid']
        # Original certificates once per reconstructed geometry; repaired witnesses
        # against every associated input, preserving the earlier 30-input check scope.
        if geometry_id not in checked or geometry_id in repairs:
            check_certificate(x, certificate); certificate_checks += 1
        checked.add(geometry_id); verified_map[identifier] = row
        input_counts[row['classification']] += 1
        changed_inputs += row['original_classification'] != row['classification']
        if len(seen) % 1000 == 0:
            print('Reconstructed resolved support inputs', len(seen), '/75070', flush=True)
    assert seen == set(recipes) and checked == set(certificates) and len(checked) == 15270
    original_settings = (Path(plan['inventory'])/'setting_input_map.jsonl').open()
    resolved_settings = (root/'resolved_setting_support.jsonl').open()
    from itertools import zip_longest
    settings = changed_settings = 0; setting_counts = Counter(); keys = set()
    for old_line, new_line in zip_longest(original_settings, resolved_settings):
        assert old_line is not None and new_line is not None
        old = json.loads(old_line); new = json.loads(new_line)
        row = verified_map[old['fit_input_id']]
        expected = dict(**old, geometry_id=row['geometry_id'], original_support_classification=row['original_classification'],
                        support_classification=row['classification'], certificate_source_kind=row['source_kind'])
        assert new == expected
        key = tuple(old[k] for k in ['guide', 'policy', 'scenario_id', 'mask', 'cohort', 'screen', 'target_order', 'background_order', 'variant', 'outcome'])
        assert key not in keys; keys.add(key)
        settings += 1; changed_settings += row['original_classification'] != row['classification']
        setting_counts[row['classification']] += 1
    geometry_counts = Counter()
    for identifier, original in certificates.items():
        cert = repairs[identifier]['candidate'] if identifier in repairs else original['result']
        geometry_counts[cert['classification']] += 1
    assert settings == receipt['settings'] == 414720
    assert changed_settings == receipt['changed_settings'] == 112
    assert changed_inputs == receipt['changed_inputs'] == 30 and receipt['changed_geometries'] == 6
    assert receipt['inputs'] == len(seen) and receipt['geometries'] == len(checked)
    assert receipt['input_classification_counts'] == dict(input_counts)
    assert receipt['geometry_classification_counts'] == dict(geometry_counts)
    assert receipt['setting_classification_counts'] == dict(setting_counts)
    assert certificate_checks == 15264+30
    verify()
    anchors = [args.plan, Path(config['launch']), Path(launch['plan']), rp,
               root/'resolved_input_support.jsonl', root/'resolved_setting_support.jsonl',
               support/'joint_support.jsonl', witnesses/'witnesses.jsonl', Path(__file__), Path('scripts/check_joint_support_certificate.py')]
    result = dict(status='passed_full_resolved_whole_protein_support_readback', source_receipt_sha256=sha(rp),
                  source_hashes={str(p):sha(p) for p in anchors}, producer_terminal_state=state,
                  inputs=len(seen), geometries=len(checked), settings=settings, changed_inputs=changed_inputs,
                  changed_settings=changed_settings, certificate_checks=certificate_checks,
                  input_classification_counts=dict(input_counts), setting_classification_counts=dict(setting_counts),
                  scope='Every input matrix, ordered identity, covariate geometry and certificate choice reconstructed, all settings compared field by field. Original certificates checked per exact geometry and repaired witnesses against all associated inputs using compensated sums. No optimizer called. Support registry is verified; numerical hull inclusion is not interior overlap, model adequacy, calibrated uncertainty or biological acceptance.')
    with Path(config['output']).open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'source_hashes'}), flush=True)


if __name__ == '__main__':
    main()
