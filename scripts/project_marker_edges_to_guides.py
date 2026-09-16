#!/usr/bin/env python3
"""Map marker splits to full-guide edges without distributing collapsed-path change."""
import argparse
import hashlib
import json
from collections import defaultdict,Counter
from pathlib import Path
from Bio import Phylo
from audit_joint_path_uncertainty import checked,rows,sha
from run_paired_marker_fits import tree_edges
from prepare_paired_phylogenetic_inputs import write_table


def canonical(side,universe):
    return min((tuple(sorted(side)),tuple(sorted(universe-side))),key=lambda s:(len(s),s))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--guides',nargs=2,type=Path,required=True)
    p.add_argument('--fit-audits',nargs='+',type=Path,required=True)
    p.add_argument('--fits',nargs='+',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or len(a.fit_audits)!=len(a.fits):raise ValueError('Fresh output and paired sources required')
    guide_edges={};full_taxa=None;edge_table=[];pins={}
    for g in a.guides:
        r=checked(g)
        if r['status']!='passed_full_species_guide_readback':raise ValueError('Audited guide required')
        tree=g/'guide.treefile';taxa={t.name for t in Phylo.read(tree,'newick').get_terminals()}
        if full_taxa is not None and taxa!=full_taxa:raise ValueError('Guide tip universes differ')
        full_taxa=taxa;edges=tree_edges(tree,taxa);guide_edges[g.name]=edges;pins[str(g/'receipt.json')]=sha(g/'receipt.json')
        for split,length in edges.items():
            eid=hashlib.sha256(','.join(split).encode()).hexdigest()
            edge_table.append(dict(guide=g.name,full_edge_id=eid,split_taxa=','.join(split),split_size=len(split),guide_branch_length=length))
    output=[];summaries=[];paired=[]
    for audit,fit in zip(a.fit_audits,a.fits):
        r=checked(audit)
        if r['status']!='complete_paired_fit_audit' or r['fit_receipt_sha256']!=sha(fit/'receipt.json'):raise ValueError('Fit audit lineage differs')
        pins[str(audit/'receipt.json')]=sha(audit/'receipt.json');pins[str(fit/'receipt.json')]=sha(fit/'receipt.json')
        groups=defaultdict(list)
        for row in rows(audit/'paired_branches.tsv'):groups[row['marker']].append(row)
        if len(groups)!=r['markers']:raise ValueError('Marker count differs')
        for marker,branches in sorted(groups.items()):
            tree=fit/marker/'aa.treefile';fr=json.loads((fit/marker/'aa.receipt.json').read_text())
            if sha(tree)!=fr['artifacts']['aa.treefile']:raise ValueError('Changed marker topology')
            pins[str(tree)]=sha(tree);taxa={t.name for t in Phylo.read(tree,'newick').get_terminals()}
            if not taxa<=full_taxa:raise ValueError('Marker taxa outside guide')
            marker_edges=tree_edges(tree,taxa)
            if set(marker_edges)!={tuple(x['split_taxa'].split(',')) for x in branches}:raise ValueError('Audited split universe differs')
            projected={}
            for guide,edges in guide_edges.items():
                lookup=defaultdict(list)
                for side,length in edges.items():
                    restricted=set(side)&taxa
                    if not restricted or restricted==taxa:continue
                    lookup[canonical(restricted,taxa)].append(hashlib.sha256(','.join(side).encode()).hexdigest())
                if len(lookup)!=2*len(taxa)-3:raise ValueError('Binary pruned-guide edge count differs')
                projected[guide]=lookup;count=Counter()
                for split in sorted(marker_edges):
                    ids=sorted(lookup.get(split,[]));status='discordant_with_pruned_guide' if not ids else 'unique_full_guide_edge' if len(ids)==1 else 'collapsed_full_guide_path'
                    count[status]+=1
                    output.append(dict(cohort=audit.name,marker=marker,marker_taxa=len(taxa),marker_split_taxa=','.join(split),marker_split_size=len(split),marker_terminal=len(split)==1,guide=guide,status=status,full_guide_edges=len(ids),full_edge_ids_json=json.dumps(ids,separators=(',',':'))))
                summaries.append(dict(cohort=audit.name,marker=marker,guide=guide,marker_taxa=len(taxa),marker_edges=len(marker_edges),**{k:count[k] for k in ['unique_full_guide_edge','collapsed_full_guide_path','discordant_with_pruned_guide']}))
            for split in sorted(marker_edges):
                x,y=[set(projected[g].get(split,[])) for g in guide_edges]
                status=('unique_same_full_split_in_both_guides' if len(x)==len(y)==1 and x==y else
                        'same_collapsed_path_splits_in_both_guides' if len(x)>1 and x==y else
                        'compatible_with_both_but_mapping_differs' if x and y else
                        'compatible_with_one_guide_only' if x or y else 'discordant_with_both_guides')
                paired.append(dict(cohort=audit.name,marker=marker,marker_split_taxa=','.join(split),marker_terminal=len(split)==1,status=status))
    a.output.mkdir(parents=True)
    for name,data in [('guide_edges.tsv',edge_table),('edge_projection.tsv',output),('marker_summary.tsv',summaries),('guide_sensitivity.tsv',paired)]:write_table(a.output/name,data)
    result=dict(status='complete_marker_edge_projection_to_provisional_guides',guide_taxa=len(full_taxa),cohorts=len(a.fit_audits),marker_cohorts=len(summaries)//2,marker_edges=len(paired),guide_edge_projection_rows=len(output),guide_sensitivity_counts=dict(Counter(r['status'] for r in paired)),source_pins=pins,script_sha256=sha(Path(__file__)),
        interpretation='Topology correspondence only. Multiple restricted full-guide splits indicate a collapsed path; no marker branch length is assigned to any constituent edge. Discordance and guide sensitivity retained. Guides are unsupported homogeneous-model trees, not final species phylogeny. No rooted clade, branch-time, acceleration, reconciliation or coupling inference. Independent pruning readback remains required.',
        artifacts={f.name:sha(f) for f in a.output.iterdir()})
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_pins','artifacts']},indent=2))


if __name__=='__main__':main()
