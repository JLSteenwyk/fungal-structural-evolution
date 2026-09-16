#!/usr/bin/env python3
"""Audit fungal/outgroup separation without choosing a root position or final tree."""
import argparse,csv,json,hashlib
from pathlib import Path
from collections import Counter,defaultdict
from Bio import Phylo
from audit_joint_path_uncertainty import checked,rows,sha
from run_paired_marker_fits import tree_edges
from prepare_paired_phylogenetic_inputs import write_table


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['manifest','snapshot','readback','trees','output']:ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--guides',nargs=2,type=Path,required=True);a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    manifest=list(rows(a.manifest));roles={x['taxon_id']:x['study_role'] for x in manifest}
    if len(roles)!=len(manifest) or set(roles.values())!={'ingroup','outgroup'}:raise ValueError('Invalid taxon-role manifest')
    sr=checked(a.snapshot);rr=checked(a.readback)
    if rr['status']!='passed_all_snapshot_input_and_graph_split_readbacks' or rr['snapshot_receipt_sha256']!=sha(a.snapshot/'receipt.json'):raise ValueError('Snapshot readback differs')
    pins={str(a.manifest):sha(a.manifest),str(a.snapshot/'receipt.json'):sha(a.snapshot/'receipt.json'),str(a.readback/'receipt.json'):sha(a.readback/'receipt.json')}
    sources=[]
    for g in a.guides:
        if checked(g)['status']!='passed_full_species_guide_readback':raise ValueError('Guide audit required')
        pins[str(g/'receipt.json')]=sha(g/'receipt.json');sources.append(('guide',g.name,g/'guide.treefile'))
    for row in sr['inputs']:
        p=a.trees/row['marker']/'tree.treefile'
        if sha(p)!=row['tree_sha256']:raise ValueError('Changed marker tree')
        sources.append(('marker',row['marker'],p))
    records=[];graph_checks=0
    for kind,name,path in sources:
        pins[str(path)]=sha(path);tree=Phylo.read(path,'newick');taxa={t.name for t in tree.get_terminals()}
        if not taxa<=roles.keys() or (kind=='guide' and taxa!=roles.keys()):raise ValueError('Unexpected tree taxa')
        ingroup={t for t in taxa if roles[t]=='ingroup'};outgroup=taxa-ingroup
        edges=tree_edges(path,taxa);side=min([tuple(sorted(ingroup)),tuple(sorted(outgroup))],key=lambda x:(len(x),x));present=bool(ingroup and outgroup and side in edges)
        # Independent edge deletion and component traversal, without clade descendant sets.
        nodes=list(tree.find_clades());index={n:i for i,n in enumerate(nodes)};graph=defaultdict(list);branches=[];leaves={index[t]:t.name for t in tree.get_terminals()}
        for parent in nodes:
            for child in parent.clades:
                u,v=index[parent],index[child];graph[u].append(v);graph[v].append(u);branches.append((u,v,child))
        matching=[]
        for u,v,child in branches:
            visited={u};pending=[v];component=set()
            while pending:
                n=pending.pop()
                if n in visited:continue
                visited.add(n)
                if n in leaves:component.add(leaves[n])
                pending.extend(graph[n])
            if component and component!=taxa and (component==ingroup or component==outgroup):matching.append(child)
        if bool(matching)!=present:raise ValueError('Independent separation check differs')
        if present and abs(sum(n.branch_length for n in matching)-edges[side])>1e-10:raise ValueError('Separating edge length differs')
        graph_checks+=1
        informative=len(ingroup)>=2 and len(outgroup)>=2
        status='insufficient_role_coverage' if not informative else 'exact_fungal_outgroup_split' if present else 'fungal_outgroup_split_absent'
        supports={n.confidence for n in matching if n.confidence is not None}
        if len(supports)>1:raise ValueError('Conflicting split supports')
        support=next(iter(supports)) if supports else None
        records.append({'tree_type':kind,'tree_id':name,'fungal_entries':len(ingroup),'outgroup_entries':len(outgroup),'status':status,'separating_edge_length':edges[side] if present else '', 'marker_sh_alrt_percent':support if kind=='marker' and informative and present and support is not None else '', 'outgroup_taxa_json':json.dumps(sorted(outgroup)),'fungal_taxa_sha256':hashlib.sha256(','.join(sorted(ingroup)).encode()).hexdigest(),'tree_sha256':sha(path)})
    a.output.mkdir(parents=True);write_table(a.output/'root_split_diagnostics.tsv',records)
    counts=Counter(x['status'] for x in records if x['tree_type']=='marker')
    result={'status':'complete_fungal_outgroup_split_diagnostic','guides':2,'completed_marker_trees':sr['completed_markers'],'planned_marker_trees':sr['planned_markers'],'marker_status_counts':dict(counts),'independent_graph_checks':graph_checks,'manifest_role_counts':dict(Counter(roles.values())),'source_pins':pins,'script_sha256':sha(Path(__file__)),'interpretation':'Exact role-separating unrooted bipartition, independently checked by graph traversal. Marker assessment requires >=2 taxa in each role; singleton outgroups give trivial terminal splits. SH-aLRT is marker support, not species-tree support. No root coordinate along the edge, polarity of events, absolute date, final species topology or reconciliation established. Incomplete marker set and entry-level taxonomy uncertainties remain.','artifacts':{'root_split_diagnostics.tsv':sha(a.output/'root_split_diagnostics.tsv')}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_pins'},indent=2))

if __name__=='__main__':main()
