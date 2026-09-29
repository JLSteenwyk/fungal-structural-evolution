#!/usr/bin/env python3
"""Resolve full input/setting support with separately verified witness provenance."""
import argparse
from collections import Counter
import json
from pathlib import Path
from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    roots = {key: Path(plan[key]) for key in ['support', 'witnesses', 'inputs', 'inventory']}
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    receipts = {}
    for key, root in roots.items():
        receipts[key] = json.loads((root/'receipt.json').read_text())
        bindings[str(root/'receipt.json')] = sha(root/'receipt.json')
        bindings.update({str(root/name): digest for name, digest in receipts[key]['artifacts'].items()})
    proof = json.loads(Path(plan['support_proof']).read_text())
    assert proof['status'] == 'passed_full_whole_protein_joint_support_readback'
    assert proof['source_receipt_sha256'] == sha(roots['support']/'receipt.json')
    repair = json.loads(Path(plan['witness_proof']).read_text())
    assert repair['status'] == 'completed_serialized_whole_protein_support_witness_checks'
    assert repair['receipt_sha256'] == sha(roots['witnesses']/'receipt.json')
    assert repair['artifacts'] == receipts['witnesses']['artifacts']
    for kind, key in [('input', 'inputs'), ('inventory', 'inventory')]:
        p = json.loads(Path(plan[kind+'_proof']).read_text())
        assert p['status'] == ('passed_full_whole_protein_materialized_input_readback' if kind == 'input' else 'passed_full_whole_protein_input_inventory_readback')
        assert p['source_receipt_sha256'] == sha(roots[key]/'receipt.json')
    bindings.update(receipts['witnesses']['source_hashes'])
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    def read_unique(path, key):
        result = {}
        for line in path.open():
            row = json.loads(line)
            assert row[key] not in result
            result[row[key]] = row
        return result
    certificates = read_unique(roots['support']/'joint_support.jsonl', 'geometry_id')
    witnesses = read_unique(roots['witnesses']/'witnesses.jsonl', 'geometry_id')
    mapping = read_unique(roots['support']/'input_support_map.jsonl', 'fit_input_id')
    inputs = read_unique(roots['inputs']/'input_manifest.jsonl', 'fit_input_id')
    assert set(mapping) == set(inputs) and len(inputs) == proof['unique_inputs'] == 75070
    assert len(certificates) == proof['unique_geometries'] == 15270
    unresolved = {k for k, r in certificates.items() if r['result']['classification'].startswith('unresolved')}
    assert set(witnesses) == unresolved and len(witnesses) == repair['accepted_geometries'] == 6
    selected = {}; geometry_counts = Counter()
    for identifier, row in certificates.items():
        if identifier in witnesses:
            witness = witnesses[identifier]
            assert witness['original'] == row['result'] and witness['specification'] == row['specification']
            assert witness['representative_input'] == row['representative_input'] and witness['accepted']
            result = witness['candidate']
            kind = 'verified_repaired_witness'
            path = roots['witnesses']/'witnesses.jsonl'
        else:
            result = row['result']; kind = 'original_verified_certificate'
            path = roots['support']/'joint_support.jsonl'
        assert result['classification'] == 'zero_supported_to_numeric_tolerance' and result['certificate_valid']
        selected[identifier] = dict(classification=result['classification'], source_kind=kind,
                                    certificate_file=str(path), certificate_file_sha256=bindings[str(path)])
        geometry_counts[result['classification']] += 1
    output = Path(plan['output']); output.mkdir(parents=True, exist_ok=False)
    resolved = {}; input_counts = Counter(); changed_inputs = 0
    with (output/'resolved_input_support.jsonl').open('w') as f:
        for identifier, old in sorted(mapping.items()):
            geometry_id = old['geometry_id']; original = certificates[geometry_id]
            assert old['classification'] == original['result']['classification']
            item = inputs[identifier]; spec = item['recipe']['specification']
            assert original['specification']['records'] == spec['records']
            assert original['specification']['columns'] == spec['columns'][2:]
            assert original['specification']['ordered_identity_sha256'] == spec['ordered_identity_sha256']
            new = selected[geometry_id]
            row = dict(fit_input_id=identifier, input_sha256=item['sha256'], geometry_id=geometry_id,
                       original_classification=old['classification'], **new)
            resolved[identifier] = row
            input_counts[new['classification']] += 1
            changed_inputs += new['classification'] != old['classification']
            f.write(json.dumps(row, sort_keys=True)+'\n')
    counts = Counter(); settings = changed_settings = 0; seen = set()
    with (output/'resolved_setting_support.jsonl').open('w') as f:
        for line in (roots['inventory']/'setting_input_map.jsonl').open():
            original = json.loads(line)
            key = tuple(original[k] for k in ['guide', 'policy', 'scenario_id', 'mask', 'cohort', 'screen', 'target_order', 'background_order', 'variant', 'outcome'])
            assert key not in seen; seen.add(key)
            identifier = original['fit_input_id']
            assert identifier is not None and identifier in resolved
            support = resolved[identifier]
            row = dict(**original, geometry_id=support['geometry_id'],
                       original_support_classification=support['original_classification'],
                       support_classification=support['classification'],
                       certificate_source_kind=support['source_kind'])
            f.write(json.dumps(row, sort_keys=True)+'\n')
            counts[support['classification']] += 1; settings += 1
            changed_settings += support['original_classification'] != support['classification']
    assert settings == receipts['inventory']['settings'] == 414720
    assert changed_inputs == repair['accepted_inputs'] == 30 and changed_settings == 112
    verify()
    receipt = dict(status='complete_resolved_whole_protein_support_pending_readback',
                   source_hashes=bindings, plan_sha256=sha(args.plan), inputs=len(resolved), geometries=len(selected), settings=settings,
                   changed_geometries=len(witnesses), changed_inputs=changed_inputs, changed_settings=changed_settings,
                   geometry_classification_counts=dict(geometry_counts), input_classification_counts=dict(input_counts),
                   setting_classification_counts=dict(counts),
                   artifacts={name:sha(output/name) for name in ['resolved_input_support.jsonl','resolved_setting_support.jsonl']},
                   scope='Full covariate-support registry retaining original and corrected certificate provenance. All 75070 inputs and 414720 settings mapped. Independent matrix/identity/certificate readback required before consumption. Numerical zero-reference hull inclusion does not establish interior overlap, dense support, causal exchangeability, adequate models or calibrated inference.')
    (output/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_hashes','artifacts']}), flush=True)


if __name__ == '__main__':
    main()
