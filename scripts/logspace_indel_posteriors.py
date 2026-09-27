#!/usr/bin/env python3
"""Independent reversible binary tree messages entirely in log probability."""
import math
import numpy as np
from scipy.special import logsumexp


def infer(tree,sequences,gain,loss,rates):
    nodes=list(tree.find_clades(order='preorder'))
    assert len({n.name for n in nodes})==len(nodes)
    width=len(next(iter(sequences.values())))
    graph={n:[] for n in nodes};parent={tree.root:None}
    for n in nodes:
        for c in n.clades:
            graph[n].append((c,c.branch_length));graph[c].append((n,c.branch_length));parent[c]=n
    pi=np.array([loss,gain])/(gain+loss);logpi=np.log(pi)
    evidence={}
    for node in nodes:
        e=np.zeros((2,width))
        if node.is_terminal():
            letters=np.array(list(sequences[node.name]))
            for state in [0,1]:
                e[state,(letters!='?') & (letters!=str(state))]=-np.inf
        evidence[node]=e
    joints={n:[] for n in nodes}
    for rate in rates:
        transitions={}
        for node in nodes:
            for other,length in graph[node]:
                t=(gain+loss)*rate*length
                with np.errstate(divide='ignore'):
                    logchange=np.log(-np.expm1(-t))
                transitions[node,other]=np.array([
                    [np.logaddexp(logpi[0],logpi[1]-t),logpi[1]+logchange],
                    [logpi[0]+logchange,np.logaddexp(logpi[1],logpi[0]-t)]])
        messages={}
        def send(source,destination):
            local=evidence[source].copy()
            for neighbor,_ in graph[source]:
                if neighbor is not destination:
                    local+=messages[neighbor,source]
            p=transitions[destination,source]
            messages[source,destination]=np.array([
                np.logaddexp(p[s,0]+local[0],p[s,1]+local[1]) for s in [0,1]])
        for node in reversed(nodes[1:]):
            send(node,parent[node])
        for node in nodes:
            for child in node.clades:
                send(node,child)
        for node in nodes:
            local=evidence[node].copy()+logpi[:,None]
            for neighbor,_ in graph[node]:
                local+=messages[neighbor,node]
            joints[node].append(local-math.log(len(rates)))
    posterior={};site=None
    for node in nodes:
        joint=logsumexp(joints[node],axis=0)
        total=np.logaddexp(joint[0],joint[1])
        if site is None:
            site=total
        assert np.max(abs(total-site))<1e-8
        posterior[node.name]=np.exp(joint[1]-total)
    return posterior,site
