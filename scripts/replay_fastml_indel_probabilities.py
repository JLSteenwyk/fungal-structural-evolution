#!/usr/bin/env python3
"""Independent binary CTMC inside/outside replay of serialized FastML results."""
import csv
import json
import math
import re
from io import StringIO
from pathlib import Path
import numpy as np
from Bio import Phylo, SeqIO
from scipy.linalg import expm
from scipy.special import gammainc, logsumexp
from scipy.stats import gamma
from prepare_case_ancestral_neighborhoods import sha


def infer(tree, sequences, gain, loss, rates):
    """Return all node marginals and unconditioned site log likelihoods."""
    width = len(next(iter(sequences.values())))
    order = list(tree.find_clades(order='preorder'))
    pi = np.array([loss, gain])/(gain+loss)
    q = np.array([[-gain, gain], [loss, -loss]])
    logs, marginals = [], []
    for rate in rates:
        transition = {n:expm(q*rate*n.branch_length) for n in order if n is not tree.root}
        inside, scales = {}, {}
        for node in reversed(order):
            part = np.ones((2, width))
            if node.is_terminal():
                letters = np.array(list(sequences[node.name]))
                part = np.array([(letters == '?') | (letters == str(s)) for s in [0,1]], dtype=float)
            scale = np.zeros(width)
            for child in node.clades:
                part *= transition[child] @ inside[child]
                scale += scales[child]
                norm = part.max(axis=0)
                assert np.all(norm > 0)
                part /= norm
                scale += np.log(norm)
            inside[node], scales[node] = part, scale
        logs.append(np.log(pi @ inside[tree.root]) + scales[tree.root])
        outside = {tree.root:np.repeat(pi[:,None], width, axis=1)}
        conditional = {}
        for node in order:
            joint = outside[node]*inside[node]
            conditional[node.name] = joint[1]/joint.sum(axis=0)
            for child in node.clades:
                part = outside[node].copy()
                for sibling in node.clades:
                    if sibling is child:
                        continue
                    part *= transition[sibling] @ inside[sibling]
                    part /= part.max(axis=0)
                message = transition[child].T @ part
                outside[child] = message/message.max(axis=0)
        marginals.append(conditional)
    logs = np.array(logs)
    site = logsumexp(logs, axis=0)-math.log(len(rates))
    weights = np.exp(logs-logsumexp(logs, axis=0))
    return {n.name:sum(weights[k]*m[n.name] for k,m in enumerate(marginals)) for n in order}, site


def analytic_check():
    # Independently enumerate all internal states for a two-internal-node tree.
    from itertools import product
    tree = Phylo.read(StringIO('(a:0.2,(b:0.3,c:0.4)I:0.1)R;'), 'newick')
    seq = {'a':'01?', 'b':'1?0', 'c':'001'}
    g,l = .4,.9
    pi = np.array([l,g])/(g+l)
    q = np.array([[-g,g],[l,-l]])
    rates = [.5,1.5]
    observed, logs = infer(tree,seq,g,l,rates)
    nodes = list(tree.find_clades())
    for pos in range(3):
        total = 0.
        ones = {n.name:0. for n in nodes}
        for rate in rates:
            for values in product([0,1], repeat=len(nodes)):
                assignment = dict(zip(nodes, values))
                if any(seq[n.name][pos]!='?' and int(seq[n.name][pos])!=assignment[n] for n in tree.get_terminals()):
                    continue
                prob = pi[assignment[tree.root]]/len(rates)
                for n in nodes:
                    for child in n.clades:
                        prob *= expm(q*rate*child.branch_length)[assignment[n],assignment[child]]
                total += prob
                for n in nodes:
                    ones[n.name] += prob*assignment[n]
        assert abs(logs[pos]-math.log(total)) < 1e-12
        assert max(abs(observed[n][pos]-v/total) for n,v in ones.items()) < 1e-12


def main():
    analytic_check()
    pp = Path('metadata/ancestral_fastml_indel_plan_20260927.json')
    plan = json.loads(pp.read_text())
    rows = []
    for rp in sorted(Path(plan['output']).glob('*/receipt.json')):
        receipt = json.loads(rp.read_text())
        job = receipt['job']
        if not job['character_count']:
            continue
        assert receipt['exit_code'] == 0 and receipt['plan_sha256'] == sha(pp)
        for rel,digest in receipt['artifacts'].items():
            assert sha(rp.parent/rel) == digest
        folder = rp.parent/'RESULTS'
        text = (folder/'EstimatedParameters.txt').read_text()
        pars = {k:float(v) for k,v in re.findall(r'^(_\w+)\s+([\deE.+-]+)',text,re.M)}
        alpha = pars['_userAlphaRate']
        rates = 4*np.diff(gammainc(alpha+1, alpha*gamma.ppf([0,.25,.5,.75,1],a=alpha,scale=1/alpha)))
        tree = Phylo.read(folder/'TheTree.INodes.ph','newick')
        assert tree.root.comment == 'N1'
        tree.root.name = tree.root.comment
        seq = {r.id:str(r.seq)+'0' for r in SeqIO.parse(job['characters'],'fasta')}
        post, likelihood = infer(tree,seq,pars['_userGain'],pars['_userLoss'],rates)
        maximum = 0.
        with (folder/'AncestralReconstructPosterior.txt').open() as handle:
            for r in csv.DictReader(handle,delimiter='\t'):
                maximum = max(maximum, abs(float(r['Prob'])-post[r['Node']][int(r['POS'])-1]))
        conditioned = likelihood[:-1]-np.log(-np.expm1(likelihood[-1]))
        reported = float(re.search(r'Log-likelihood=\s*([\deE.+-]+)',text).group(1))
        diagnostic = {}
        if abs(conditioned.sum()-reported) > .01:
            initial = Phylo.read(rp.parent/'input_tree.nwk','newick')
            for index,node in enumerate(initial.get_nonterminals()):
                node.name = f'initial_{index}'
            _,initial_likelihood = infer(initial,seq,pars['_userGain'],pars['_userLoss'],rates)
            stale = float((likelihood[:-1]-np.log(-np.expm1(initial_likelihood[-1]))).sum())
            diagnostic = dict(initial_tree_denominator_log_likelihood=stale,
                              initial_tree_denominator_difference=stale-reported,
                              disposition='likelihood_discrepancy_requires_resolution')
        rows.append(dict(job_id=job['job_id'],maximum_probability_difference=maximum,
            corrected_log_likelihood=float(conditioned.sum()),reported_log_likelihood=reported,
            likelihood_difference=float(conditioned.sum()-reported),alpha=alpha,
            gain=pars['_userGain'],loss=pars['_userLoss'],gamma_category_rates=rates.tolist(),
            source_receipt_sha256=sha(rp), **diagnostic))
        print(job['job_id'],maximum,conditioned.sum()-reported,flush=True)
    output = Path('metadata/fastml_indel_probability_replay_snapshot_20260927.json')
    output.write_text(json.dumps(dict(status='serialized_parameter_replay_snapshot',
        analytic_enumeration='passed',jobs=rows,script_sha256=sha(__file__),plan_sha256=sha(pp),
        scope='Independent inside/outside replay from rounded serialized parameters and tree. Differences reported without asserting fit convergence, a full-precision match, model adequacy or complete production.'),indent=2)+'\n')


if __name__ == '__main__':
    main()
