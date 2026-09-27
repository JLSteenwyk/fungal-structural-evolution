#!/usr/bin/env python3
"""Stable binary inside/outside probabilities, conditional on supplied model parameters."""
import math
import numpy as np
from scipy.special import logsumexp


def transition_matrix(gain,loss,time):
    pi0,pi1=loss/(gain+loss),gain/(gain+loss)
    changed=-np.expm1(-(gain+loss)*time)
    return np.array([[1-pi1*changed,pi1*changed],[pi0*changed,1-pi0*changed]])


def infer(tree, sequences, gain, loss, rates):
    """Return all node marginals and unconditioned site log likelihoods."""
    width = len(next(iter(sequences.values())))
    order = list(tree.find_clades(order='preorder'))
    pi = np.array([loss, gain])/(gain+loss)
    assert len({n.name for n in order}) == len(order), "unique node names required"
    logs, marginals = [], []
    for rate in rates:
        transition = {n:transition_matrix(gain, loss, rate*n.branch_length) for n in order if n is not tree.root}
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

