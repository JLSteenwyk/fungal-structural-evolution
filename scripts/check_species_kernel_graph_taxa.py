#!/usr/bin/env python3
"""Verify complete matched-graph taxon coverage in the audited species-kernel grid."""
import json,hashlib
from pathlib import Path
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
root=Path('results/phylogeny/species-distance-kernels-20260927-v1');r=json.loads((root/'receipt.json').read_text());ap=Path('metadata/species_distance_kernel_readback_20260927.json');a=json.loads(ap.read_text());assert a['source_receipt_sha256']==sha(root/'receipt.json') and a['status']=='passed_full_species_distance_kernel_readback'
assert sha(root/'taxa.json')==r['artifacts']['taxa.json'];tips=set(json.loads((root/'taxa.json').read_text()))
graph=Path('results/orthology/background-match-graph-20260927-v1');gr=json.loads((graph/'receipt.json').read_text());ga=Path('metadata/background_match_graph_completed_readback_20260927.json');assert json.loads(ga.read_text())['producer_receipt_sha256']==sha(graph/'receipt.json')
used=set();counts={};bindings={}
for role in ['target','background']:
 path=graph/(role+'_nodes.jsonl');assert sha(path)==gr['artifacts'][path.name];bindings[str(path)]=sha(path);role_taxa=set();n=0
 for line in path.open():
  row=json.loads(line);role_taxa.update([row['taxon_id']] if role=='target' else [row['taxon_a'],row['taxon_b']]);n+=1
 assert role_taxa<=tips;used.update(role_taxa);counts[role]={'nodes':n,'taxa':len(role_taxa)}
r=dict(status='passed_full_matched_graph_species_kernel_taxon_coverage',tree_taxa=len(tips),graph_taxa=len(used),graph_taxon_ids=sorted(used),role_counts=counts,source_node_bindings=bindings,kernel_readback_sha256=sha(ap),graph_readback_sha256=sha(ga),script_sha256=sha(__file__),scope='Every target/background graph-node taxon is present in the complete common tip grid of all five audited kernels. Covers the full matching graph, not only selected or structurally qualified records. No chosen covariance model or phylogenetic fit.')
with Path('metadata/species_kernel_graph_taxon_coverage_20260927.json').open('x') as f:json.dump(r,f,indent=2);f.write('\n')
print({k:v for k,v in r.items() if k not in ['graph_taxon_ids','source_node_bindings']})
