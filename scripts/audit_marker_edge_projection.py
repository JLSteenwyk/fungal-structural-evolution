#!/usr/bin/env python3
"""Independently prune guides and verify marker projection compatibility/path lengths."""
import argparse
import copy
import json
import math
from collections import defaultdict
from pathlib import Path
from Bio import Phylo
from audit_joint_path_uncertainty import checked,rows,sha


def split_lengths(tree,taxa):
    values=defaultdict(float)
    if {t.name for t in tree.get_terminals()}!=taxa:raise ValueError('Pruned taxa differ')
    for node in tree.find_clades():
        if node is tree.root:continue
        side=frozenset(t.name for t in node.get_terminals());other=frozenset(taxa-side)
        key=min((tuple(sorted(side)),tuple(sorted(other))),key=lambda s:(len(s),s))
        values[key]+=node.branch_length or 0.
    return dict(values)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--projection',type=Path,required=True)
    p.add_argument('--guides',type=Path,nargs=2,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=checked(a.projection)
    for file,h in r['source_pins'].items():
        if sha(Path(file))!=h:raise ValueError('Changed source pin')
    trees={};full_taxa={}
    for folder in a.guides:
        checked(folder)
        if str(folder/'receipt.json') not in r['source_pins']:raise ValueError('Unpinned guide')
        trees[folder.name]=Phylo.read(folder/'guide.treefile','newick');full_taxa[folder.name]={t.name for t in trees[folder.name].get_terminals()}
    marker_trees={}
    for file in r['source_pins']:
        f=Path(file)
        if f.name=='aa.treefile':
            key=(f.parent.parent.name.replace('paired-marker-fits-','paired-fit-audit-'),f.parent.name)
            if key in marker_trees:raise ValueError('Ambiguous marker tree')
            marker_trees[key]=f
    full_edges={(x['guide'],x['full_edge_id']):x for x in rows(a.projection/'guide_edges.tsv')}
    groups=defaultdict(list)
    for row in rows(a.projection/'edge_projection.tsv'):groups[row['cohort'],row['marker'],row['guide']].append(row)
    checks=0;compatible=0;pruned=0
    for (cohort,marker,guide),data in sorted(groups.items()):
        mt=Phylo.read(marker_trees[cohort,marker],'newick');taxa={t.name for t in mt.get_terminals()};expected_marker=split_lengths(mt,taxa)
        if {tuple(x['marker_split_taxa'].split(',')) for x in data}!=set(expected_marker) or len(data)!=len(expected_marker):raise ValueError('Marker grid differs')
        tree=copy.deepcopy(trees[guide])
        for tip in sorted(full_taxa[guide]-taxa):tree.prune(tip)
        edges=split_lengths(tree,taxa);pruned+=1
        for row in data:
            split=tuple(row['marker_split_taxa'].split(','));ids=json.loads(row['full_edge_ids_json'])
            if len(ids)!=len(set(ids)) or len(ids)!=int(row['full_guide_edges']):raise ValueError('Path edge identity/count differs')
            if split not in edges:
                if ids or row['status']!='discordant_with_pruned_guide':raise ValueError('Discordance differs')
            else:
                if not ids or row['status']!=('unique_full_guide_edge' if len(ids)==1 else 'collapsed_full_guide_path'):raise ValueError('Compatibility differs')
                path_length=math.fsum(float(full_edges[guide,e]['guide_branch_length']) for e in ids)
                if not math.isclose(path_length,edges[split],rel_tol=1e-10,abs_tol=1e-10):raise ValueError('Pruned length differs from projected path sum')
                compatible+=1
            checks+=1
        if pruned%25==0:print('Pruned guide/marker combinations checked',pruned,flush=True)
    if checks!=r['guide_edge_projection_rows']:raise ValueError('Incomplete readback')
    out=dict(status='passed_full_pruning_projection_readback',pruned_guide_marker_combinations=pruned,projection_rows=checks,compatible_path_sums_checked=compatible,
        projection_receipt_sha256=sha(a.projection/'receipt.json'),script_sha256=sha(Path(__file__)),
        scope='Independent repeated Bio.Phylo taxon pruning, retained marker split universes, guide compatibility/dispositions and all compatible branch-length sums checked. Does not validate guide biology, support, branch-time allocation or causal/evolutionary inference.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
