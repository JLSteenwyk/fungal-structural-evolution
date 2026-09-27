#!/usr/bin/env python3
"""Independent all-edge network-flow verification for one saved bootstrap tree."""
import csv,gzip,hashlib,json
from collections import Counter
from io import StringIO
from pathlib import Path
import numpy as np
from Bio import Phylo
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def audit_tree(job):
    label,number,newick,states,taxa,folder,expected_hash=job
    path=Path(folder)/f'{label}-{number:04d}.json.gz';assert sha(path)==expected_hash
    with gzip.open(path,'rt') as f:data=json.load(f)
    tree=Phylo.read(StringIO(newick),'newick');nodes=list(tree.find_clades(order='preorder'));idx={v:i for i,v in enumerate(nodes)};n=len(nodes);inf=n+1
    assert len(tree.get_terminals())==len(taxa) and {v.name for v in tree.get_terminals()}==set(taxa)
    endpoints={(idx[v],idx[c]):c for v in nodes for c in v.clades};keys={};sides={}
    for pair,node in endpoints.items():
        side=sorted(v.name for v in node.get_terminals());canonical=min([side,sorted(set(taxa)-set(side))],key=lambda x:(len(x),x));key=hashlib.sha256(json.dumps(canonical,separators=(',',':')).encode()).hexdigest();keys[pair]=key;sides[key]=';'.join(canonical)
    assert len(set(keys.values()))==len(endpoints)
    assert len(data)==2 and {d['coding'] for d in data}=={'original32_coding','ramaria_unknown'}
    summaries=[];contributions=[];checked=0
    for d in data:
        current={t:s for t,s in states.items() if d['coding']!='ramaria_unknown' or t!='F113071'};rr=[];cc=[];vv=[]
        for parent,child in endpoints:rr.extend([parent,child]);cc.extend([child,parent]);vv.extend([1,1])
        def constraint(node,state):return (n,node) if state==0 else (node,n+1)
        for node in tree.get_terminals():
            if node.name in current:
                x,y=constraint(idx[node],current[node.name]);rr.append(x);cc.append(y);vv.append(inf)
        def cost(extra):
            x=list(rr);y=list(cc);v=list(vv)
            for node,state in extra:
                i,j=constraint(node,state);x.append(i);y.append(j);v.append(inf)
            matrix=coo_matrix((np.array(v,dtype=np.int64),(x,y)),shape=(n+2,n+2)).tocsr();value=int(maximum_flow(matrix,n,n+1).flow_value)
            return value if value<inf else None
        baseline=cost([]);assert d['minimum_changes']==baseline;seen=set();counts=Counter()
        for key,parent,child,status,stored in d['edges']:
            pair=parent,child;assert pair not in seen and keys[pair]==key;seen.add(pair)
            expected=[cost([(parent,s),(child,t)]) for s in (0,1) for t in (0,1)]
            assert stored==expected,(label,number,d['coding'],pair,stored,expected);checked+=4
            same=any(expected[i]==baseline for i in [0,3]);change=any(expected[i]==baseline for i in [1,2])
            correct='optional_change' if same and change else 'required_change' if change else 'no_change_in_any_optimum';assert status==correct
            counts[status]+=1;contributions.append([d['coding'],key,status,sides[key]])
        assert seen==set(endpoints)
        summaries.append(dict(tree_source=label,tree_index=number,coding=d['coding'],minimum_changes=baseline,edges=len(endpoints),**{k:counts[k] for k in ['required_change','optional_change','no_change_in_any_optimum']}))
    assert sha(path)==expected_hash
    return dict(tree_source=label,tree_index=number,shard_sha256=expected_hash,constrained_costs=checked,summaries=summaries,contributions=contributions)
