#!/usr/bin/env python3
"""Write and independently check new witnesses for every unresolved geometry."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from check_joint_support_certificate import check_certificate
from repair_joint_support_weights import repair_weights
from screen_duplication_alignment_reuse import sha


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    support, inputs = Path(plan['support']), Path(plan['inputs'])
    for root, proof_path, launch_path in [
        (support, plan['support_proof'], plan['support_check_launch']),
        (inputs, plan['input_proof'], plan['input_check_launch']),
    ]:
        proof = json.loads(Path(proof_path).read_text())
        assert proof['status'].startswith('passed_full_whole_protein_')
        assert proof['source_receipt_sha256'] == sha(root / 'receipt.json')
        launch = json.loads(Path(launch_path).read_text())
        state = dict(line.split('=', 1) for line in subprocess.check_output([
            'systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState',
            '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
        receipt = json.loads((root / 'receipt.json').read_text())
        bindings[str(root / 'receipt.json')] = sha(root / 'receipt.json')
        bindings.update({str(root / name): h for name, h in receipt['artifacts'].items()})

    def verify_bindings():
        for path, expected in bindings.items():
            assert sha(path) == expected, path

    verify_bindings()
    unresolved = {}
    for line in (support / 'joint_support.jsonl').open():
        row = json.loads(line)
        if row['result']['classification'].startswith('unresolved'):
            assert row['geometry_id'] not in unresolved
            unresolved[row['geometry_id']] = row
    assert len(unresolved) == plan['expected_geometries']
    selected = {}
    for line in (support / 'input_support_map.jsonl').open():
        row = json.loads(line)
        if row['geometry_id'] in unresolved:
            assert row['fit_input_id'] not in selected
            selected[row['fit_input_id']] = row['geometry_id']
    assert len(selected) == plan['expected_inputs']
    recipes = {}
    for line in (inputs / 'input_manifest.jsonl').open():
        row = json.loads(line)
        if row['fit_input_id'] in selected:
            assert row['fit_input_id'] not in recipes
            recipes[row['fit_input_id']] = row
    assert set(recipes) == set(selected)

    def load(identifier):
        item = recipes[identifier]
        path = inputs / item['path']
        assert sha(path) == item['sha256']
        bindings[str(path)] = item['sha256']
        spec = item['recipe']['specification']
        assert digest(spec) == identifier
        with np.load(path, allow_pickle=False) as arrays:
            matrix, identity = arrays['matrix'], arrays['row_identity']
        assert hashlib.sha256(matrix.tobytes()).hexdigest() == spec['values_sha256']
        identity_sha = hashlib.sha256(identity.tobytes()).hexdigest()
        assert identity_sha == spec['ordered_identity_sha256']
        x = np.ascontiguousarray(matrix[:, 2:], dtype='<f8')
        geometry = dict(records=len(x), columns=spec['columns'][2:],
                        values_sha256=hashlib.sha256(x.tobytes()).hexdigest(),
                        ordered_identity_sha256=identity_sha)
        assert digest(geometry) == selected[identifier]
        assert geometry == unresolved[selected[identifier]]['specification']
        return x

    output = Path(plan['output'])
    output.mkdir(parents=True, exist_ok=False)
    witnesses = []
    for geometry_id, original in sorted(unresolved.items()):
        identifier = original['representative_input']
        assert selected[identifier] == geometry_id
        result = repair_weights(load(identifier), original['result'])
        witnesses.append(dict(geometry_id=geometry_id, representative_input=identifier,
                              specification=original['specification'], **result))
    path = output / 'witnesses.jsonl'
    path.write_text(''.join(json.dumps(r, sort_keys=True, allow_nan=False)+'\n' for r in witnesses))

    # Read serialized results and independently reconstruct their barycenters
    # against every associated input, not only the representative input.
    saved = {r['geometry_id']: r for r in map(json.loads, path.read_text().splitlines())}
    assert set(saved) == set(unresolved) and len(saved) == len(witnesses)
    accepted_inputs = 0
    for identifier, geometry_id in selected.items():
        row = saved[geometry_id]
        assert row['original'] == unresolved[geometry_id]['result']
        check_certificate(load(identifier), row['candidate'])
        assert row['accepted'] == (row['candidate']['classification'] == 'zero_supported_to_numeric_tolerance')
        accepted_inputs += row['accepted']
    verify_bindings()
    result = dict(status='completed_serialized_whole_protein_support_witness_checks',
                  geometries=len(saved), accepted_geometries=sum(r['accepted'] for r in saved.values()),
                  associated_inputs_checked=len(selected), accepted_inputs=accepted_inputs,
                  source_hashes=bindings, artifacts={'witnesses.jsonl': sha(path)},
                  scope='Separate witnesses checked with independent compensated summation against every associated audited input. Original certificates unchanged. No interior overlap, calibrated inference or model adequacy established. Downstream support overlay integration remains separate.')
    (output / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes', 'artifacts']}))


if __name__ == '__main__':
    main()
