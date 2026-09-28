#!/usr/bin/env python3
"""Paired six-run regression of MP ascertainment caching, with independent replay."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess
import numpy as np
from Bio import Phylo, SeqIO
from scipy.special import gammainc
from scipy.stats import gamma
from ancestral_chain_attempt import sha, write_json
from replay_fastml_indel_probabilities import infer, analytic_check


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    build_path = Path('results/software/fastml-precision-cache-build-20260928-v1/receipt.json')
    build = json.loads(build_path.read_text())
    plan_path = Path('metadata/ancestral_fastml_indel_plan_20260927.json')
    plan = json.loads(plan_path.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest
    for row in build['builds']:
        assert sha(row['binary']) == row['binary_sha256'] and sha(row['log']) == row['log_sha256']
    jobs = [j for j in plan['jobs'] if j['job_id'] in {
        'OG0000294-envelope-famsa-terminal_gap', 'OG0000294-envelope-mafft-terminal_gap',
        'OG0000294-envelope-mafft-terminal_unknown'}]
    assert len(jobs) == 3 and all(j['character_count'] == 1 for j in jobs)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    write_json(out / 'plan.json', dict(build_receipt_sha256=sha(build_path), source_plan_sha256=sha(plan_path),
        jobs=[j['job_id'] for j in jobs], variants=[r['label'] for r in build['builds']],
        cpus=1, memory_limit_gib=2, per_run_timeout_seconds=120, output_allowance_mib=100,
        planning_minutes=[1, 20], paid_resources=False, gpu=False,
        qualification=dict(precision_only_minimum_absolute_likelihood_discrepancy=0.4,
                           refreshed_maximum_absolute_likelihood_discrepancy=1e-6,
                           maximum_probability_discrepancy=1e-6), script_sha256=sha(__file__)))
    analytic_check()
    rows = []
    cpu = str(min(os.sched_getaffinity(0)))
    for job in jobs:
        original = Path(plan['output']) / job['job_id']
        for variant in build['builds']:
            folder = out / (job['job_id'] + '__' + variant['label'])
            folder.mkdir()
            settings = dict(plan['settings'], _seqFile=str(Path(job['characters']).resolve()),
                            _treeFile=str((original / 'input_tree.nwk').resolve()), _outDir=str(folder / 'RESULTS'))
            parameters = folder / 'parameters.txt'
            parameters.write_text(''.join(f'{k} {v}\n' for k, v in settings.items()))
            command = ['taskset', '-c', cpu, 'prlimit', '--as=2147483648',
                       str(Path(variant['binary']).resolve()), str(parameters)]
            with (folder / 'stdout.log').open('w') as handle:
                run = subprocess.run(command, cwd=folder, stdout=handle, stderr=subprocess.STDOUT, timeout=120)
            assert run.returncode == 0
            results = folder / 'RESULTS'
            text = (results / 'EstimatedParameters.txt').read_text()
            parameters_fit = {k: float(v) for k, v in re.findall(r'^(_\w+)\s+([\deE.+-]+)', text, re.M)}
            alpha, gain, loss = [parameters_fit[k] for k in ['_userAlphaRate', '_userGain', '_userLoss']]
            rates = 4 * np.diff(gammainc(alpha + 1, alpha * gamma.ppf([0, .25, .5, .75, 1], a=alpha, scale=1/alpha)))
            tree = Phylo.read(results / 'TheTree.INodes.ph', 'newick')
            assert tree.root.comment == 'N1'
            tree.root.name = tree.root.comment
            sequences = {r.id: str(r.seq) + '0' for r in SeqIO.parse(job['characters'], 'fasta')}
            post, ll = infer(tree, sequences, gain, loss, rates)
            replay = float((ll[:-1] - np.log(-np.expm1(ll[-1]))).sum())
            reported = float(re.search(r'Log-likelihood=\s*([\deE.+-]+)', text).group(1))
            maximum = 0.
            seen = set()
            with (results / 'AncestralReconstructPosterior.txt').open() as handle:
                for r in csv.DictReader(handle, delimiter='\t'):
                    key = r['Node'], int(r['POS'])
                    assert key not in seen and key[1] == 1
                    seen.add(key)
                    maximum = max(maximum, abs(float(r['Prob']) - post[key[0]][0]))
            assert len(seen) == len(post)
            row = dict(job_id=job['job_id'], variant=variant['label'], command=command,
                       reported_log_likelihood=reported, replay_log_likelihood=replay,
                       likelihood_difference=replay-reported, maximum_probability_difference=maximum,
                       parameters=dict(alpha=alpha, gain=gain, loss=loss),
                       artifacts={str(p.relative_to(folder)): sha(p) for p in folder.rglob('*') if p.is_file()})
            write_json(folder / 'receipt.json', row)
            rows.append(row)
            print(job['job_id'], variant['label'], replay-reported, maximum, flush=True)
    passed = all(r['maximum_probability_difference'] < 1e-6 and
                 (abs(r['likelihood_difference']) > .4 if r['variant'] == 'precision-only'
                  else abs(r['likelihood_difference']) < 1e-6) for r in rows)
    write_json(out / 'receipt.json', dict(status='passed_paired_cache_regression' if passed else 'failed_paired_cache_regression',
        runs=rows, plan_sha256=sha(out / 'plan.json'), analytic_enumeration='passed',
        scope='Three diagnosed one-character cases, two isolated variants. Does not qualify other '
              'inputs, optimizer routes, model adequacy or ancestral structure ensembles.'))
    assert passed


if __name__ == '__main__':
    main()
