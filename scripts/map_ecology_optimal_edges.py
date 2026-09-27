#!/usr/bin/env python3
"""Map possible and required changes under an explicitly frozen trait coding."""
import argparse,csv,json,hashlib
from collections import Counter
from pathlib import Path
from Bio import Phylo
from ecology_optimal_edge_states import edge_costs

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text())
 for path,h in p['pins'].items():assert sha(path)==h
 with open(p['evidence']) as f:evidence=list(csv.DictReader(f,delimiter='\t'))
 states={r['taxon_id']:1 if r['state']=='ectomycorrhizal' else 0 for r in evidence if r['state'] in ['ectomycorrhizal','asymbiotic','saprotrophic']}
 with open(p['manifest']) as f:taxa={r['taxon_id'] for r in csv.DictReader(f,delimiter='\t')}
 out=Path(p['output']);out.mkdir(parents=True,exist_ok=False);summaries=[];rows=[]
 for source in p['trees']:
  receipt=json.loads(Path(source['receipt']).read_text());proof=json.loads(Path(source['readback']).read_text())
  assert proof['source_receipt_sha256']==sha(source['receipt']) and proof['status']=='passed_pmsf_profile_tree_and_bootstrap_readback'
  assert sha(source['tree'])==receipt['artifacts'][Path(source['tree']).name]
  tree=Phylo.read(source['tree'],'newick');assert len(tree.get_terminals())==len(taxa) and {n.name for n in tree.get_terminals()}==taxa
  for mode,omit in [('original32_coding',''),('ramaria_unknown','F113071')]:
   current={t:s for t,s in states.items() if t!=omit};score,edges=edge_costs(tree,current);counts=Counter();nodes=list(tree.find_clades(order='preorder'));idx={n:i for i,n in enumerate(nodes)};inf=len(nodes)+1
   for parent,child,cost in edges:
    side=sorted(n.name for n in child.get_terminals());other=sorted(taxa-set(side));canonical=min([side,other],key=lambda v:(len(v),v));split=hashlib.sha256(json.dumps(canonical,separators=(',',':')).encode()).hexdigest()
    unchanged=any(cost[i]==score for i in [0,3]);changed=any(cost[i]==score for i in [1,2]);status='optional_change' if unchanged and changed else 'required_change' if changed else 'no_change_in_any_optimum';counts[status]+=1
    row=dict(tree_source=source['label'],coding=mode,split_id=split,parent_index=idx[parent],child_index=idx[child],child_taxa=';'.join(side),minimum_changes=score,status=status,optimal_assignments=';'.join(k for k,v in zip(['00','01','10','11'],cost) if v==score))
    row.update({f'cost_{k}':v if v<inf else '' for k,v in zip(['00','01','10','11'],cost)});rows.append(row)
   summaries.append(dict(tree_source=source['label'],coding=mode,minimum_changes=score,coded_tips=len(current),unknown_tips=len(taxa)-len(current),edges=len(edges),counts=dict(counts)))
 with (out/'edge_states.tsv').open('w') as f:
  w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 for path,h in p['pins'].items():assert sha(path)==h
 r=dict(status='complete_ecology_optimal_edge_diagnostic_pending_readback',plan_sha256=sha(a.plan),rows=len(rows),summaries=summaries,artifacts={'edge_states.tsv':sha(out/'edge_states.tsv')},scope='All edges on two completed ML trees, original 32-species conditional binary coding and Ramaria omission. No new 45-species labels coded as ECM absence. Root-free equal-cost optimum; endpoint order follows stored Newick traversal and is NOT ancestral gain/loss direction. Possible assignments are jointly constrained across edges, not independent changes or probabilities. No bootstrap branch support, ecological replication or association effect established; infeasible endpoint costs blank.')
 (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
