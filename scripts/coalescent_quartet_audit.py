#!/usr/bin/env python3
"""Independent colored-tree quartet DP and coalescent numerical checks.

No ASTRAL libraries or native quartet algorithm are used. A rooted postorder
DP counts color-distinct quartets by their distribution among LCA children:
3+1, 2+2, and 2+1+1 are resolved; 1+1+1+1 is unresolved. Rooted triples
carry their cherry pair, making the 3+1 contribution independent of rooting.
Integer arithmetic precedes gene-specific normalization. The complement
calculation for the global denominator counts unresolved quartet centers.
"""
import hashlib
import itertools
import json
import math
import re

import dendropy
import numpy as np
from numba import njit
from scipy.special import betaln, betaincc, logsumexp


PARTITIONS = np.array([(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)], dtype=np.int64)
TRIPLES = np.array([(a, b, c) for a in range(4) for b in range(a+1, 4)
                    for c in range(4) if c not in (a, b)], dtype=np.int64)
TRIPLE_INDEX = np.full((4, 4, 4), -1, dtype=np.int64)
for i, (a, b, c) in enumerate(TRIPLES):
    TRIPLE_INDEX[a, b, c] = TRIPLE_INDEX[b, a, c] = i


def read_tree(path=None, data=None):
    return dendropy.Tree.get(path=str(path) if path is not None else None,
                            data=data, schema='newick', rooting='force-unrooted',
                            preserve_underscores=True)


def packed(tree, tips):
    """Postorder children and leaf indices; multifurcations are retained."""
    nodes = list(tree.postorder_node_iter())
    position = {n: i for i, n in enumerate(nodes)}
    tip_index = {t: i for i, t in enumerate(tips)}
    children, offsets, leaves = [], [0], []
    found = []
    for node in nodes:
        kids = list(node.child_node_iter())
        children.extend(position[c] for c in kids)
        offsets.append(len(children))
        if not kids:
            label = node.taxon.label
            assert label in tip_index, label
            leaves.append(tip_index[label])
            found.append(label)
        else:
            assert node.taxon is None
            leaves.append(-1)
    assert len(found) == len(set(found)) >= 4
    return (np.array(children, dtype=np.int64), np.array(offsets, dtype=np.int64),
            np.array(leaves, dtype=np.int64))


@njit(cache=True)
def colored_quartets(children, offsets, leaves, colors):
    size = len(leaves)
    counts = np.zeros((size, 4), dtype=np.int64)
    triples = np.zeros((size, 12), dtype=np.int64)
    quartets = np.zeros((size, 3), dtype=np.int64)
    for node in range(size):
        if leaves[node] >= 0:
            counts[node, colors[leaves[node]]] = 1
            continue
        begin, end = offsets[node], offsets[node+1]
        same_child_pairs = np.zeros((4, 4), dtype=np.int64)
        for k in range(begin, end):
            child = children[k]
            for a in range(4):
                counts[node, a] += counts[child, a]
                for b in range(4):
                    same_child_pairs[a, b] += counts[child, a]*counts[child, b]
            for j in range(12):
                triples[node, j] += triples[child, j]
            for j in range(3):
                quartets[node, j] += quartets[child, j]
        for k in range(begin, end):
            child = children[k]
            c = counts[child]
            outside = counts[node]-c
            for j in range(12):
                a, b, other = TRIPLES[j]
                triples[node, j] += c[a]*c[b]*outside[other]
            for j in range(3):
                a, b, d, e = PARTITIONS[j]
                # 3+1: the resolved cherry of a rooted triple and its outsider.
                quartets[node, j] += triples[child, TRIPLE_INDEX[a,b,d]]*outside[e]
                quartets[node, j] += triples[child, TRIPLE_INDEX[a,b,e]]*outside[d]
                quartets[node, j] += triples[child, TRIPLE_INDEX[d,e,a]]*outside[b]
                quartets[node, j] += triples[child, TRIPLE_INDEX[d,e,b]]*outside[a]
                # 2+2: color a,b in this child; d,e in another child.
                quartets[node, j] += c[a]*c[b]*(same_child_pairs[d,e]-c[d]*c[e])
                # 2+1+1: the two outsiders must come from different children.
                quartets[node, j] += c[a]*c[b]*(outside[d]*outside[e]-(same_child_pairs[d,e]-c[d]*c[e]))
                quartets[node, j] += c[d]*c[e]*(outside[a]*outside[b]-(same_child_pairs[a,b]-c[a]*c[b]))
    available = counts[-1, 0]*counts[-1, 1]*counts[-1, 2]*counts[-1, 3]
    return quartets[-1], available


def resolved_quartet_denominator(tree):
    """All resolved induced quartets, including partial and polytomous genes."""
    nodes = list(tree.postorder_node_iter())
    descendants = {}
    for node in nodes:
        descendants[node] = (1 if node.is_leaf() else
                             sum(descendants[c] for c in node.child_node_iter()))
    n = descendants[tree.seed_node]
    unresolved = 0
    for node in nodes:
        components = [descendants[c] for c in node.child_node_iter()]
        if node.parent_node is not None:
            components.append(n-descendants[node])
        coefficients = [1, 0, 0, 0, 0]
        for count in components:
            for j in range(4, 0, -1):
                coefficients[j] += count*coefficients[j-1]
        unresolved += coefficients[4]
    assert 0 <= unresolved <= math.comb(n, 4)
    return math.comb(n, 4)-unresolved


def posterior(counts, prior=0.5):
    counts = np.array(counts, dtype=float)
    assert np.isfinite(counts).all() and min(counts) >= 0 and prior > 0
    total = float(counts.sum())
    a, b = counts+1, total-counts+2*prior
    logs = counts*math.log(2)+betaln(a,b)+np.log(betaincc(a,b,1/3))
    result = np.exp(logs-logsumexp(logs))
    assert np.isfinite(result).all()
    return result


def map_length(concordant, total, prior=0.5):
    assert 0 <= concordant <= total+1e-8 and prior > 0
    return max(0.0, -math.log(1.5*(1-concordant/(total+2*prior))))


def graph(tree):
    nodes = list(tree.preorder_node_iter())
    adjacency = {n: set(n.child_nodes()) for n in nodes}
    for n in nodes:
        if n.parent_node is not None: adjacency[n].add(n.parent_node)
    root = tree.seed_node
    if len(adjacency[root]) == 2:
        a, b = adjacency.pop(root)
        adjacency[a].remove(root); adjacency[b].remove(root)
        adjacency[a].add(b); adjacency[b].add(a)
    assert all(len(neighbors) in (1,3) for neighbors in adjacency.values()), 'Nonbinary species tree'
    return adjacency


def component_tips(adjacency, start, blocked):
    result, stack = set(), [(start, blocked)]
    while stack:
        node, previous = stack.pop()
        if node.taxon is not None: result.add(node.taxon.label)
        stack.extend((other,node) for other in adjacency[node] if other is not previous)
    return result


def annotated_branches(tree, expected_tips):
    """Species internal branches, canonical partitions and original annotations."""
    tips = sorted(expected_tips)
    found = [n.taxon.label for n in tree.leaf_node_iter()]
    assert len(found) == len(set(found)) == len(tips) and set(found) == set(tips)
    index = {t: i for i,t in enumerate(tips)}
    adjacency = graph(tree)
    records = []
    seen = set()
    for node in tree.preorder_node_iter():
        if node.label is None:
            assert node.edge_length is None or node.edge_length == 0, 'Unexpected unannotated length'
            continue
        assert node is not tree.seed_node and node.taxon is None
        match = re.fullmatch(r'\[(.*)\]', node.label)
        assert match, node.label
        items = [field.split('=') for field in match.group(1).split(';')]
        assert len(items) == 11 and len({p[0] for p in items}) == 11
        values = {k:float(v) for k,v in items}
        assert set(values) == {'q1','q2','q3','f1','f2','f3','pp1','pp2','pp3','QC','EN'}
        other = node.parent_node
        if other not in adjacency:
            neighbors = tree.seed_node.child_nodes()
            assert len(neighbors) == 2
            other = next(n for n in neighbors if n is not node)
        assert other in adjacency[node] and len(adjacency[node]) == len(adjacency[other]) == 3
        groups = []
        for center, opposite in [(node,other),(other,node)]:
            pair = [component_tips(adjacency,n,center) for n in adjacency[center] if n is not opposite]
            groups.extend(sorted(pair,key=lambda s:tuple(sorted(s))))
        assert len(groups) == 4 and set.union(*groups) == set(tips)
        assert sum(map(len,groups)) == len(tips)
        colors = np.full(len(tips),-1,dtype=np.int64)
        for color, group in enumerate(groups):
            for taxon in group: colors[index[taxon]]=color
        side = sorted(groups[0]|groups[1])
        complement = sorted(set(tips)-set(side))
        side = min(side,complement,key=lambda s:(len(s),tuple(s)))
        split = hashlib.sha256(json.dumps(side,separators=(',',':')).encode()).hexdigest()
        assert split not in seen; seen.add(split)
        qc = math.prod(map(len,groups))
        assert values['QC'] == qc and qc > 0
        assert node.edge_length is not None and math.isfinite(node.edge_length) and node.edge_length >= 0
        records.append(dict(split=split,canonical_side=side,groups=[sorted(g) for g in groups],
                            colors=colors,native=values,length=node.edge_length))
    assert len(records) == len(tips)-3
    return sorted(records,key=lambda r:r['split'])


def close(actual, expected, tolerance=2e-8):
    assert math.isfinite(actual) and math.isfinite(expected)
    assert math.isclose(actual,expected,rel_tol=tolerance,abs_tol=tolerance), (actual,expected)


def numerical_check(values, length, independent_counts, tolerance=2e-8):
    """NNI alternatives compared as an unordered pair; no invented polarity."""
    counts = np.array([values['f1'],values['f2'],values['f3']])
    independent_counts = np.array(independent_counts,dtype=float)
    assert np.isfinite(counts).all() and min(counts) >= 0
    close(counts[0],independent_counts[0],tolerance)
    for actual,expected in zip(sorted(counts[1:]),sorted(independent_counts[1:])):
        close(actual,expected,tolerance)
    total = float(independent_counts.sum())
    close(values['EN'],total,tolerance)
    expected_pp = posterior(counts)
    for i in range(3):
        # Native zero-evidence q=NaN is an explicit unavailable fraction.
        if total == 0: assert math.isnan(values['q'+str(i+1)])
        else: close(values['q'+str(i+1)],counts[i]/total,tolerance)
        close(values['pp'+str(i+1)],expected_pp[i],tolerance)
    close(length,map_length(independent_counts[0],total),tolerance)
    close(sum(values['pp'+str(i+1)] for i in range(3)),1,tolerance)
    return dict(independent_effective_genes=total,
                native_and_independent_counts=counts.tolist(),
                independent_counts_canonical_groups=independent_counts.tolist(),
                independent_posterior=expected_pp.tolist(),
                independent_map_length=map_length(independent_counts[0],total),
                q_available=total>0)
