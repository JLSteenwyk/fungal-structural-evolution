#!/usr/bin/env python3
"""Prepare full independent-seed model inputs without launching sampling."""
import argparse
import collections
import json
from pathlib import Path
from ancestral_chain_attempt import sha, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    pp = Path('metadata/baliphy_prior_initialization_plan_20260927.json')
    proof_path = Path('metadata/baliphy_prior_initialization_final_readback_completed_20260927.json')
    plan = json.loads(pp.read_text())
    proof = json.loads(proof_path.read_text())
    audit_path = Path(proof['audit_receipt'])
    assert sha(audit_path) == proof['audit_receipt_sha256']
    audit = json.loads(audit_path.read_text())
    assert audit['dispositions'] == {'initial_model_semantics_and_finite_scores_checked': 405}
    for name, digest in audit['pins'].items():
        assert sha(name) == digest, name
    for name, digest in audit['artifacts'].items():
        assert sha(audit_path.parent / name) == digest, name
    assert len(plan['jobs']) == 405
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / 'models').mkdir()
    rows = []
    pins = {str(p): sha(p) for p in [pp, proof_path, audit_path, Path(__file__),
            Path('scripts/ancestral_chain_attempt.py')]}
    for job in sorted(plan['jobs'], key=lambda j: j['job_id']):
        folder = Path(plan['output']) / job['job_id']
        receipt_path = folder / 'receipt.json'
        assert sha(receipt_path) == audit['pins'][str(receipt_path)]
        source = json.loads(receipt_path.read_text())
        assert source['job'] == job and source['exit_code'] == 0
        for name, digest in source['artifacts'].items():
            assert sha(folder / name) == digest
        pins[str(receipt_path)] = sha(receipt_path)
        original_path = folder / 'BAliPhy.Main.hs'
        original = original_path.read_text()
        marker = ';mcmcState <- makeMCMCState'
        assert original.count(marker) == 1
        insertion = ';T.writeFile (outputDirectory </> "runtime-tree.nwk") (writeNewick_rooted (addInternalLabels tree))\n'
        mapped = original.replace(marker, insertion + marker)
        program = out / 'models' / (job['job_id'] + '.hs')
        program.write_text(mapped)
        for chain in range(1, 5):
            rows.append(dict(chain_id=job['job_id'] + '-chain' + str(chain),
                seed=20261001 + len(rows), chain=chain,
                effective_input_group=job['effective_input_group'], prior_label=job['prior_label'],
                original_configuration_ids=job['original_configuration_ids'],
                family=job['family'], proteins=job['proteins'],
                alignment=job['alignment'], alignment_sha256=job['alignment_sha256'],
                tree=job['tree'], tree_sha256=job['tree_sha256'],
                source_receipt=str(receipt_path), source_receipt_sha256=sha(receipt_path),
                original_program_sha256=sha(original_path),
                program=str(program), program_sha256=sha(program)))
    assert len(rows) == len({r['seed'] for r in rows}) == len({r['chain_id'] for r in rows}) == 1620
    groups = collections.Counter((r['effective_input_group'], r['prior_label']) for r in rows)
    assert len(groups) == 405 and set(groups.values()) == {4}
    assert len({alias for r in rows for alias in r['original_configuration_ids']}) == 324
    write_json(out / 'chain_inputs.json', rows)
    receipt = dict(status='prepared_full_independent_chain_inputs_not_launched',
        chains=len(rows), effective_inputs=135, priors=3, chains_per_prior_input=4,
        original_configuration_aliases=324, pins=pins,
        artifacts={str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file()},
        required_before_launch=['full short-run audit and resource summary',
            'costed iteration horizon and runtime/memory/output limits',
            'controller with attempt recovery and sample mapping audit',
            'independent-chain convergence and extension policy'],
        scope='Prepared model sources and unique seeds only. All aliases retained; '
              'no posterior samples, convergence, root sensitivity or project completion claim.')
    write_json(out / 'receipt.json', receipt)
    print(json.dumps({k: v for k, v in receipt.items() if k not in ['pins', 'artifacts']}))


if __name__ == '__main__':
    main()
