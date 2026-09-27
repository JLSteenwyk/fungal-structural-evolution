#!/usr/bin/env python3
"""Exact per-edge binary assignments among all minimum-cost tree labelings."""

def edge_costs(tree,states):
    nodes=list(tree.find_clades(order='preorder'));inside={};outside={tree.root:(0,0)};inf=len(nodes)+1
    for node in reversed(nodes):
        if node.is_terminal():
            s=states.get(node.name);inside[node]=tuple(0 if s is None or s==v else inf for v in (0,1))
        else:inside[node]=tuple(sum(min(inside[c][t]+(s!=t) for t in (0,1)) for c in node.clades) for s in (0,1))
    optimum=min(inside[tree.root]);edges=[]
    for parent in nodes:
        for child in parent.clades:
            rest=tuple(outside[parent][s]+sum(min(inside[c][t]+(s!=t) for t in (0,1)) for c in parent.clades if c is not child) for s in (0,1))
            outside[child]=tuple(min(rest[s]+(s!=t) for s in (0,1)) for t in (0,1))
            costs=tuple(rest[s]+(s!=t)+inside[child][t] for s,t in [(0,0),(0,1),(1,0),(1,1)])
            assert min(costs)==optimum
            edges.append((parent,child,costs))
    return optimum,edges
