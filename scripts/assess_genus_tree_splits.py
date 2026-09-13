#!/usr/bin/env python3
"""Check genus-label separability in a frozen set of audited unrooted gene trees."""
import argparse
import json
from collections import defaultdict
from pathlib import Path
from Bio import Phylo
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT, sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def incompatible(a, g, universe):
    """Two unrooted splits conflict only when all four intersections occur."""
    b=universe-a;h=universe-g
    return bool(a&g and a&h and b&g and b&h)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['audit','trees','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable split assessment')
    audit=checked_receipt(a.audit)
    if audit['status'] not in ['incomplete_snapshot','complete_support_audit']:raise ValueError('Audited gene trees required')
    membership=ROOT/'metadata/codon_group_membership.tsv'
    groups=defaultdict(set)
    for r in read_table(membership):groups[r['genus_label']].add(r['taxon_id'])
    branch_rows=defaultdict(list)
    for r in read_table(a.audit/'branch_support.tsv'):branch_rows[r['marker']].append(r)
    rows=[]
    for entry in audit['inputs']:
        marker=entry['marker'];receipt_path=a.trees/marker/'receipt.json'
        if sha(receipt_path)!=entry['receipt_sha256']:raise ValueError('Changed gene-tree receipt')
        receipt=json.loads(receipt_path.read_text());path=ROOT/receipt['tree_path']
        if sha(path)!=entry['tree_sha256']:raise ValueError('Changed gene tree')
        universe={n.name for n in Phylo.read(path,'newick').get_terminals()}
        for label,full in sorted(groups.items()):
            group=full&universe;status='assessed' if len(group)>=4 and len(universe-group)>=2 else 'insufficient_taxon_coverage'
            matching=[];conflicting=[]
            if status=='assessed':
                for b in branch_rows[marker]:
                    side=set(b['smaller_side_taxa'].split(';'))
                    if not side<=universe:raise ValueError('Split outside tree universe')
                    if side==group or universe-side==group:matching.append(b)
                    elif incompatible(side,group,universe):conflicting.append(b)
                if len(matching)>1 or (matching and conflicting):raise ValueError('Invalid split partition')
                if not matching and not conflicting:raise ValueError('Resolved tree missing a split without incompatibility')
            reported=[float(b['sh_alrt_percent']) for b in conflicting if b['sh_alrt_percent']!='']
            rows.append({'marker':marker,'genus_label':label,'tree_taxa':len(universe),'manifest_group_taxa':len(full),'group_taxa_in_tree':len(group),
                'missing_group_taxa':';'.join(sorted(full-universe)),'status':status,
                'edge_separable':bool(matching) if status=='assessed' else '',
                'separating_edge_sh_alrt':matching[0]['sh_alrt_percent'] if matching else '',
                'incompatible_edges':len(conflicting) if status=='assessed' else '',
                'incompatible_edges_sh_alrt_ge80':sum(s>=80 for s in reported) if status=='assessed' else '',
                'max_incompatible_edge_sh_alrt':max(reported) if reported else ''})
    summaries=[]
    for label in sorted(groups):
        subset=[r for r in rows if r['genus_label']==label and r['status']=='assessed']
        summaries.append({'genus_label':label,'snapshot_markers':len(audit['inputs']),'assessed_markers':len(subset),
            'edge_separable_markers':sum(r['edge_separable'] for r in subset),
            'separating_edge_sh_alrt_ge80':sum(r['edge_separable'] and r['separating_edge_sh_alrt']!='' and float(r['separating_edge_sh_alrt'])>=80 for r in subset),
            'markers_with_incompatible_edge_sh_alrt_ge80':sum(r['incompatible_edges_sh_alrt_ge80']>0 for r in subset)})
    a.output.mkdir(parents=True);write_table(a.output/'genus_marker_splits.tsv',rows);write_table(a.output/'genus_summary.tsv',summaries)
    result={'status':'complete_genus_split_assessment_of_frozen_tree_snapshot','tree_audit_sha256':sha(a.audit/'receipt.json'),
        'membership_sha256':sha(membership),'script_sha256':sha(Path(__file__)),'completed_gene_trees':len(audit['inputs']),
        'planned_gene_trees':audit['planned_markers'],'pending_gene_trees':audit['pending_markers'],
        'genus_labels':len(groups),'marker_genus_rows':len(rows),
        'interpretation':'Edge separability of observed genus-label members in unrooted gene trees, not rooted species monophyly, orthology validation, independent transitions, or a selection test. Missing members explicit. Early finishing gene trees may be unrepresentative; SH-aLRT is not a bootstrap probability. Topological conflict can reflect biological history or inference/data artifacts.',
        'artifacts':{x.name:sha(x) for x in a.output.iterdir() if x.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summaries,indent=2))

if __name__=='__main__':main()
