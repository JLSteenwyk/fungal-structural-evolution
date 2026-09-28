#!/usr/bin/env python3
"""Replay raw probabilities with explicit binary64 boundary-roundoff reporting."""
import csv
import re
from pathlib import Path
import numpy as np
from Bio import Phylo, SeqIO
from scipy.special import gammainc
from scipy.stats import gamma
from replay_fastml_indel_probabilities import infer


BOUNDARY_TOLERANCE = 64 * np.finfo(np.float64).eps


def boundary_error(value):
    """Return raw boundary excess without clipping; reject larger/nonfinite errors.

    64 binary64 epsilon is an engineering screen, not a proof of a forward
    error bound. Counts and maximum excess remain visible in every readback.
    """
    if not np.isfinite(value):
        raise ValueError('Nonfinite probability')
    excess = max(0., -value, value - 1.)
    if excess > BOUNDARY_TOLERANCE:
        raise ValueError('Probability exceeds explicit boundary tolerance: ' + repr(value))
    return excess


def readback(folder, job):
    folder = Path(folder) / 'RESULTS'
    text = (folder / 'EstimatedParameters.txt').read_text()
    pars = {k: float(v) for k, v in re.findall(r'^(_\w+)\s+([\deE.+-]+)', text, re.M)}
    alpha, gain, loss = [pars[k] for k in ['_userAlphaRate', '_userGain', '_userLoss']]
    assert np.isfinite([alpha, gain, loss]).all() and min(alpha, gain, loss) > 0
    rates = 4 * np.diff(gammainc(alpha + 1, alpha * gamma.ppf([0, .25, .5, .75, 1], a=alpha, scale=1/alpha)))
    assert np.isfinite(rates).all() and min(rates) > 0
    tree = Phylo.read(folder / 'TheTree.INodes.ph', 'newick')
    assert tree.root.name is None and tree.root.comment == 'N1'
    tree.root.name = tree.root.comment
    nodes = [n.name for n in tree.find_clades()]
    assert None not in nodes and len(nodes) == len(set(nodes))
    observed = {r.id: str(r.seq) for r in SeqIO.parse(job['characters'], 'fasta')}
    assert len(observed) == job['proteins']
    assert set(observed) == {n.name for n in tree.get_terminals()}
    assert {len(s) for s in observed.values()} == {job['character_count']}
    seq = {k: v + '0' for k, v in observed.items()}
    post, ll = infer(tree, seq, gain, loss, rates)
    assert np.isfinite(ll).all() and ll[-1] < 0
    assert all(np.isfinite(v).all() and ((v >= 0) & (v <= 1 + 1e-12)).all() for v in post.values())
    seen, max_difference, max_tip_error = set(), 0., 0.
    boundary_count, maximum_boundary_error = 0, 0.
    with (folder / 'AncestralReconstructPosterior.txt').open() as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            node, pos, value = r['Node'], int(r['POS']), float(r['Prob'])
            assert node in post and 1 <= pos <= job['character_count']
            assert r['State'] == '1'
            excess = boundary_error(value)
            boundary_count += int(excess > 0)
            maximum_boundary_error = max(maximum_boundary_error, excess)
            assert (node, pos) not in seen
            seen.add((node, pos))
            max_difference = max(max_difference, abs(value - post[node][pos - 1]))
            if node in observed and observed[node][pos - 1] != '?':
                max_tip_error = max(max_tip_error, abs(value - int(observed[node][pos - 1])))
    assert len(seen) == len(nodes) * job['character_count'] and max_tip_error < 1e-6
    corrected = float((ll[:-1] - np.log(-np.expm1(ll[-1]))).sum())
    reported = float(re.search(r'Log-likelihood=\s*([\deE.+-]+)', text).group(1))
    assert np.isfinite([corrected, reported]).all()
    return dict(status='integrity_and_independent_replay_complete_not_fit_qualification',
        probability_rows=len(seen), maximum_probability_difference=max_difference,
        raw_boundary_excursions=boundary_count, maximum_raw_boundary_error=maximum_boundary_error,
        probability_boundary_tolerance=BOUNDARY_TOLERANCE, raw_values_clipped=False,
        maximum_known_tip_error=max_tip_error, reported_log_likelihood=reported,
        replay_log_likelihood=corrected, likelihood_difference=corrected-reported,
        fitted_parameters=dict(alpha=alpha, gain=gain, loss=loss),
        gamma_category_rates=rates.tolist(), total_branch_length=tree.total_branch_length(),
        scope='Raw boundary excursions within 64 binary64 epsilon recorded, never clipped. Fixed fitted-parameter replay; no optimization qualification, mask-conditional '
              'ascertainment model, or homologous joint insertion/deletion model.')
