#!/usr/bin/env python3
"""Exercise new readback on all effective inputs using frozen short-run fixtures."""
import json
from pathlib import Path
import tempfile
from ancestral_chain_attempt import sha, write_json
from readback_independent_baliphy_chain import readback


def main():
    inputs = Path('results/ancestral/baliphy-independent-chain-inputs-20260927-v1/chain_inputs.json')
    rows = json.loads(inputs.read_text())
    selected = {}
    for row in rows:
        selected.setdefault(row['effective_input_group'], row)
    assert len(selected) == 135
    proof_path = Path('metadata/baliphy_sample_mapping_final_readback_completed_20260927.json')
    proof = json.loads(proof_path.read_text())
    ap = Path(proof['audit_receipt'])
    assert sha(ap) == proof['audit_receipt_sha256']
    old_audit = json.loads(ap.read_text())
    mapping = 'results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv'
    tested, samples, receipts = [], 0, {}
    with tempfile.TemporaryDirectory() as temporary:
        for group, row in selected.items():
            source = Path('results/ancestral/baliphy-sample-mapping-20260927-v1') / row['original_configuration_ids'][0]
            rp = source / 'receipt.json'
            assert sha(rp) == old_audit['pins'][str(rp)]
            receipts[str(rp)] = sha(rp)
            old = json.loads(rp.read_text())
            for name, digest in old['artifacts'].items():
                assert sha(source / name) == digest
            folder = Path(temporary) / group
            target = folder / 'independent-chain-1'
            target.mkdir(parents=True)
            original_logs = list(source.glob('mapping-check-*'))
            assert len(original_logs) == 1
            for name, origin in [('runtime-tree.nwk', source / 'runtime-tree.nwk'),
                                 ('C1.log', original_logs[0] / 'C1.log'),
                                 ('C1.P1.fastas', original_logs[0] / 'C1.P1.fastas')]:
                (target / name).symlink_to(origin.resolve())
            fixture = folder / 'receipt.json'
            write_json(fixture, dict(exit_code=0, artifacts={str(p.relative_to(folder)): sha(p) for p in target.iterdir()}))
            result = readback(row, fixture, 20, mapping)
            assert result['status'] == 'all_saved_alignments_and_candidate_nodes_checked'
            samples += len(result['candidate_samples'])
            tested.append(group)
    result = dict(status='passed_all_effective_input_short_fixture_readback',
        groups=len(tested), candidate_samples=samples,
        pins={str(inputs): sha(inputs), str(ap): sha(ap), mapping: sha(mapping),
              __file__: sha(__file__), 'scripts/readback_independent_baliphy_chain.py': sha('scripts/readback_independent_baliphy_chain.py')},
        source_receipts=receipts,
        scope='Frozen fixed-parameter short-run fixtures adapted in temporary directories. '
              'No new posterior sampling or proof of convergence.')
    output = Path('metadata/baliphy_independent_readback_fixture_check_20260927.json')
    assert not output.exists()
    write_json(output, result)
    print('Passed', len(tested), 'effective-input fixtures;', samples, 'candidate samples')


if __name__ == '__main__':
    main()
