#!/usr/bin/env python3
"""Exact matching-quartet score from independent node-component intersections.

A resolved quartet has two endpoint centers in each unrooted tree. Summing
compatible 2+1+1 selections across every species/gene node pair counts each
matching quartet four times. Component intersections either place the cherry
pair in the same species/gene component, or in disjoint paired components.
Gene multifurcations contribute only resolved 2+1+1 selections. Stars and
quartets missing any of their taxa therefore contribute zero.
"""
import numpy as np
from numba import njit

from coalescent_quartet_audit_v2 import component_tips, graph


def species_tripartitions(tree,tips):
    adjacency = graph(tree)
    index = {t:i for i,t in enumerate(tips)}
    rows = []
    for node,neighbors in adjacency.items():
        if len(neighbors) != 3: continue
        groups = sorted((component_tips(adjacency,n,node) for n in neighbors),
                        key=lambda g:tuple(sorted(g)))
        assert set.union(*groups) == set(tips) and sum(map(len,groups)) == len(tips)
        colors = np.full(len(tips),-1,dtype=np.int64)
        for color,g in enumerate(groups):
            for taxon in g: colors[index[taxon]]=color
        rows.append(colors)
    assert len(rows) == len(tips)-2
    return np.array(rows,dtype=np.int64)


@njit(cache=True)
def matching_quartets(tripartitions,children,offsets,leaves):
    nodes = len(leaves)
    counts = np.zeros((nodes,3),dtype=np.int64)
    incident = np.zeros((nodes+1,3),dtype=np.int64)
    total_intersections = np.int64(0)
    for colors in tripartitions:
        counts[:] = 0
        for node in range(nodes):
            if leaves[node] >= 0:
                counts[node,colors[leaves[node]]] = 1
            else:
                for k in range(offsets[node],offsets[node+1]):
                    child = children[k]
                    for color in range(3): counts[node,color] += counts[child,color]
        total = counts[nodes-1]
        for node in range(nodes):
            if leaves[node] >= 0: continue
            degree = 0
            for k in range(offsets[node],offsets[node+1]):
                child = children[k]
                incident[degree] = counts[child]
                degree += 1
            if node != nodes-1:
                incident[degree] = total-counts[node]
                degree += 1
            if degree < 3: continue
            squares = np.zeros((3,3),dtype=np.int64)
            for k in range(degree):
                for a in range(3):
                    for b in range(3): squares[a,b] += incident[k,a]*incident[k,b]
            for paired in range(3):
                other1,other2 = (paired+1)%3,(paired+2)%3
                for k in range(degree):
                    x = incident[k]
                    outside = total-x
                    # The two paired leaves occupy the same component in both trees.
                    distinct_outsiders = outside[other1]*outside[other2]-(squares[other1,other2]-x[other1]*x[other2])
                    total_intersections += (x[paired]*(x[paired]-1)//2)*distinct_outsiders
                    # Species and gene paired components contain complementary pairs.
                    distinct_pair = (outside[paired]*outside[paired]-(squares[paired,paired]-x[paired]*x[paired]))//2
                    total_intersections += distinct_pair*x[other1]*x[other2]
    assert total_intersections >= 0 and total_intersections%4 == 0
    return total_intersections//4
