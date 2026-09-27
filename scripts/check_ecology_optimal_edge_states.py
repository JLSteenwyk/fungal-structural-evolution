#!/usr/bin/env python3
"""Exhaustively enumerate small-tree assignments to check every edge constraint."""
from io import StringIO
from itertools import product
from Bio import Phylo
from ecology_optimal_edge_states import edge_costs
cases=0
for nwk in ['((a,b),(c,d));','(a,b,c,d);','(a,(b,(c,d)));']:
 tree=Phylo.read(StringIO(nwk),'newick');nodes=list(tree.find_clades());index={n:i for i,n in enumerate(nodes)}
 for labels in product([None,0,1],repeat=4):
  states=dict(zip('abcd',labels));score,edges=edge_costs(tree,states);best=10**9;expected=[[10**9]*4 for _ in edges]
  for bits in product([0,1],repeat=len(nodes)):
   if any(states[n.name] is not None and bits[index[n]]!=states[n.name] for n in nodes if n.is_terminal()):continue
   cost=sum(bits[index[p]]!=bits[index[c]] for p,c,_ in edges);best=min(best,cost)
   for j,(p,c,_) in enumerate(edges):
    k=2*bits[index[p]]+bits[index[c]];expected[j][k]=min(expected[j][k],cost)
  assert score==best
  for (_,_,actual),wanted in zip(edges,expected):
   # Impossible leaf constraints have a finite sentinel; feasible costs exact.
   assert [v==score for v in actual]==[v==best for v in wanted]
   for x,y in zip(actual,wanted):
    if y<10**9:assert x==y
  cases+=1
print('Passed',cases,'exhaustive state patterns, including unresolved tips and polytomies')
