#!/usr/bin/env python3
"""Verify all endpoint-constrained costs by an independent network-flow calculation."""
import argparse,csv,json,hashlib
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow
from Bio import Phylo

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text());root=Path(p['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
 assert r['plan_sha256']==sha(a.plan)
 for path,h in p['pins'].items():assert sha(path)==h
 for name,h in r['artifacts'].items():assert sha(root/name)==h
 with open(p['evidence']) as f:evidence=list(csv.DictReader(f,delimiter='\t'))
 states={x['taxon_id']:int(x['state']=='brown_rot') for x in evidence if x['state'] in ['white_rot','brown_rot']};assert len(states)==6
 with (root/'edge_states.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 seen=set();checked=0
 for source in p['trees']:
  tree=Phylo.read(source['tree'],'newick');nodes=list(tree.find_clades(order='preorder'));idx={n:i for i,n in enumerate(nodes)};n=len(nodes);inf=n+1;alltips={x.name for x in tree.get_terminals()};endpoints={(idx[parent],idx[child]):child for parent in nodes for child in parent.clades}
  for scenario in p['scenarios']:
   mode=scenario['label'];current=dict(states,**scenario['uncertain_assignments']);rr=[];cc=[];vv=[]
   for parent,child in endpoints:rr.extend([parent,child]);cc.extend([child,parent]);vv.extend([1,1])
   def constraint(node,state):return (n,node) if state==0 else (node,n+1)
   for node in tree.get_terminals():
    if node.name in current:
     x,y=constraint(idx[node],current[node.name]);rr.append(x);cc.append(y);vv.append(inf)
   def cost(extra):
    x=list(rr);y=list(cc);v=list(vv)
    for node,state in extra:
     i,j=constraint(node,state);x.append(i);y.append(j);v.append(inf)
    matrix=coo_matrix((np.array(v,dtype=np.int64),(x,y)),shape=(n+2,n+2)).tocsr()
    result=maximum_flow(matrix,n,n+1).flow_value
    return int(result) if result<inf else None
   baseline=cost([]);selected=[row for row in rows if row['tree_source']==source['label'] and row['coding']==mode];assert len(selected)==len(endpoints)
   counts=Counter()
   for row in selected:
    parent,child=int(row['parent_index']),int(row['child_index']);key=(source['label'],mode,parent,child);assert key not in seen;seen.add(key);node=endpoints[parent,child]
    side=sorted(t.name for t in node.get_terminals());assert row['child_taxa']==';'.join(side)
    canonical=min([side,sorted(alltips-set(side))],key=lambda v:(len(v),v));assert row['split_id']==hashlib.sha256(json.dumps(canonical,separators=(',',':')).encode()).hexdigest()
    expected={s+t:cost([(parent,int(s)),(child,int(t))]) for s in '01' for t in '01'}
    assert int(row['minimum_changes'])==baseline
    for k,v in expected.items():assert row['cost_'+k]==('' if v is None else str(v)),(key,k,v,row['cost_'+k]);checked+=1
    optimal=[k for k,v in expected.items() if v==baseline];assert row['optimal_assignments']==';'.join(optimal)
    same=any(k[0]==k[1] for k in optimal);change=any(k[0]!=k[1] for k in optimal)
    status='optional_change' if same and change else 'required_change' if change else 'no_change_in_any_optimum';assert row['status']==status;counts[status]+=1
   summary=next(x for x in r['summaries'] if x['tree_source']==source['label'] and x['coding']==mode)
   assert summary['counts']==dict(counts) and summary['minimum_changes']==baseline and summary['edges']==len(endpoints) and summary['coded_tips']==len(current) and summary['unknown_tips']==len(alltips)-len(current)
   print(source['label'],mode,'checked',len(selected),'edges',flush=True)
 assert len(seen)==len(rows)==r['rows'] and sha(rp)==rh
 for path,h in p['pins'].items():assert sha(path)==h
 result=dict(status='passed_full_wood_decay_optimal_edge_network_flow_readback',producer_receipt_sha256=rh,rows=len(rows),constrained_costs=checked,script_sha256=sha(__file__),scope='All four endpoint-constrained costs on every edge reconstructed by integer network maximum flow, independently of inside/outside recurrence; split identities, assignment sets and summaries verified. Conditional equal-cost diagnostic, not rooted origins or ecological effects.')
 a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
