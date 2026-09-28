#!/usr/bin/env python3
"""Prepare four-chain candidate-length traces using audited source-node identities."""
import argparse
import csv
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha, write_json
from advance_independent_chain_diagnostics import quartet_manifest


def length_rows(samples, expected):
    indexed, identities = {}, {}
    for row in samples:
        iteration, level, length = int(row['iteration']), int(row['level']), int(row['ungapped_length'])
        assert length == float(row['ungapped_length']) and length >= 0
        key = iteration, level
        assert key not in indexed and level in range(4) and iteration in expected
        indexed[key] = length
        if level in identities:
            assert identities[level] == row['source_node']
        identities[level] = row['source_node']
    assert set(indexed) == {(iteration, level) for iteration in expected for level in range(4)}
    rows = [dict(iter=iteration, **{'length_level%d' % level: indexed[iteration, level] for level in range(4)})
            for iteration in expected]
    return rows, identities


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer-plan', type=Path, required=True)
    parser.add_argument('--group', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--python', required=True)
    args = parser.parse_args()
    plan = json.loads(args.producer_plan.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest
    jobs = [job for job in json.loads(Path(plan['jobs']).read_text())
            if job['config']['model_input_identity'] == args.group]
    assert len(jobs) == 4
    ready = quartet_manifest(jobs, Path(plan['output']), plan['iterations'], sha(plan['mapping']))
    assert ready is not None, 'Four fully checked chains required'
    scalar_manifest, evidence = ready
    expected = list(range(0, plan['iterations'] + 1, 10))
    prepared, identities = [], None
    for job, original in zip(jobs, scalar_manifest['chains']):
        paths = list((Path(plan['output']) / job['chain']['chain_id']).glob('attempt-*-sample-audit.json'))
        assert len(paths) == 1 and sha(paths[0]) == evidence[str(paths[0])]
        audit = json.loads(paths[0].read_text())
        rows, nodes = length_rows(audit['candidate_samples'], expected)
        if identities is not None:
            assert identities == nodes, 'Source-node identities differ between chains'
        identities = nodes
        prepared.append((original, rows))
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    chains = []
    for original, rows in prepared:
        target = out / (original['chain_id'] + '.tsv')
        with target.open('w') as handle:
            writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        chains.append(dict(original, log=str(target), log_sha256=sha(target)))
    outputs = {}
    for fraction in plan['diagnostics']['burn_in_fractions']:
        cutoff = int(plan['iterations'] * fraction)
        manifest = dict(chains=chains, expected_iterations=expected,
            discard_through_iteration=cutoff, variables=['length_level%d' % level for level in range(4)])
        mp = out / ('manifest-%d.json' % cutoff)
        write_json(mp, manifest)
        destination = out / ('discard-%d' % cutoff)
        subprocess.run([args.python, 'scripts/ancestral_chain_diagnostics.py', '--manifest', str(mp),
                        '--output', str(destination)], check=True)
        dp = destination / 'diagnostics.json'
        outputs[str(cutoff)] = dict(path=str(dp), sha256=sha(dp),
            retained_samples_per_chain=sum(i > cutoff for i in expected))
    write_json(out / 'receipt.json', dict(status='candidate_length_screens_complete_not_posterior_qualification',
        model_input_identity=args.group, source_nodes_by_level=identities, evidence=evidence,
        producer_plan_sha256=sha(args.producer_plan), script_sha256=sha(__file__), outputs=outputs,
        artifacts={str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file()},
        scope='Audited candidate lengths at saved iterations, using stable source nodes across chains. '
              'No site-state, homology/alignment mixing, or posterior qualification; constant lengths remain review flags.'))


if __name__ == '__main__':
    main()
