#!/usr/bin/env python3
"""Replay the complete terminal FastML grid without overwriting the earlier snapshot."""
import argparse
import csv
import json
from pathlib import Path
import re
import numpy as np
from Bio import Phylo, SeqIO
from scipy.special import gammainc
from scipy.stats import gamma
from ancestral_chain_attempt import sha, write_json
from replay_fastml_indel_probabilities import infer, analytic_check


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    analytic_check()
    source = json.loads(Path(plan['source_plan']).read_text())
    proof = json.loads(Path(plan['terminal_readback']).read_text())
    source_jobs = {r['job_id']: r for r in source['jobs']}
    checked = {r['job_id']: r for r in proof['jobs']}
    assert set(source_jobs) == set(checked) and len(checked) == 156
    rows, empty = [], []
    for job_id, job in sorted(source_jobs.items()):
        folder = Path(source['output']) / job_id
        rp = folder / 'receipt.json'
        assert sha(rp) == checked[job_id]['source_receipt_sha256']
        receipt = json.loads(rp.read_text())
        assert receipt['job'] == job
        for rel, digest in receipt['artifacts'].items():
            assert sha(folder / rel) == digest
        if not job['character_count']:
            empty.append(job_id)
            continue
        assert receipt['exit_code'] == 0
        text = (folder / 'RESULTS/EstimatedParameters.txt').read_text()
        pars = {k: float(v) for k, v in re.findall(r'^(_\w+)\s+([\deE.+-]+)', text, re.M)}
        alpha, gain, loss = [pars[k] for k in ['_userAlphaRate', '_userGain', '_userLoss']]
        assert np.isfinite([alpha, gain, loss]).all() and min(alpha, gain, loss) > 0
        rates = 4 * np.diff(gammainc(alpha + 1, alpha * gamma.ppf(
            [0, .25, .5, .75, 1], a=alpha, scale=1/alpha)))
        assert np.isfinite(rates).all() and min(rates) > 0
        tree = Phylo.read(folder / 'RESULTS/TheTree.INodes.ph', 'newick')
        assert tree.root.name is None and tree.root.comment == 'N1'
        tree.root.name = tree.root.comment
        seq = {r.id: str(r.seq) + '0' for r in SeqIO.parse(job['characters'], 'fasta')}
        post, likelihood = infer(tree, seq, gain, loss, rates)
        assert np.isfinite(likelihood).all() and likelihood[-1] < 0
        assert all(np.isfinite(p).all() and ((p >= 0) & (p <= 1 + 1e-12)).all() for p in post.values())
        maximum, count = 0., 0
        with (folder / 'RESULTS/AncestralReconstructPosterior.txt').open() as handle:
            for r in csv.DictReader(handle, delimiter='\t'):
                maximum = max(maximum, abs(float(r['Prob']) - post[r['Node']][int(r['POS']) - 1]))
                count += 1
        assert count == checked[job_id]['probability_rows']
        corrected = float((likelihood[:-1] - np.log(-np.expm1(likelihood[-1]))).sum())
        reported = float(re.search(r'Log-likelihood=\s*([\deE.+-]+)', text).group(1))
        diagnostic = {}
        if abs(corrected - reported) > .01:
            initial = Phylo.read(folder / 'input_tree.nwk', 'newick')
            for index, node in enumerate(initial.get_nonterminals()):
                node.name = f'initial_{index}'
            _, initial_ll = infer(initial, seq, gain, loss, rates)
            stale = float((likelihood[:-1] - np.log(-np.expm1(initial_ll[-1]))).sum())
            diagnostic = dict(initial_tree_denominator_log_likelihood=stale,
                              initial_tree_denominator_difference=stale - reported)
        row = dict(job_id=job_id, probability_rows=count,
                   maximum_probability_difference=maximum,
                   corrected_log_likelihood=corrected, reported_log_likelihood=reported,
                   likelihood_difference=corrected - reported,
                   alpha=alpha, gain=gain, loss=loss, gamma_category_rates=rates.tolist(),
                   source_receipt_sha256=sha(rp), **diagnostic)
        write_json(out / (job_id + '.json'), row)
        rows.append(row)
        print(job_id, maximum, corrected - reported, flush=True)
    assert len(rows) == 153 and len(empty) == 3
    assert sum(r['probability_rows'] for r in rows) == proof['probability_rows']
    for path, digest in plan['pins'].items():
        assert sha(path) == digest
    write_json(out / 'receipt.json', dict(status='complete_serialized_parameter_replay_not_fit_qualification',
        plan_sha256=sha(args.plan), analytic_enumeration='passed', jobs=rows, empty_inputs=empty,
        probability_rows=sum(r['probability_rows'] for r in rows),
        artifacts={p.name: sha(p) for p in out.glob('*.json')},
        scope='Independent inside/outside replay with rounded parameters and branches for all '
              '153 nonempty inputs. All-zero exclusion uses final tree; initial-tree denominator '
              'is a diagnostic when absolute likelihood discrepancy exceeds 0.01. This cutoff '
              'does not certify smaller differences. No global optimization or model adequacy claim.'))


if __name__ == '__main__':
    main()
