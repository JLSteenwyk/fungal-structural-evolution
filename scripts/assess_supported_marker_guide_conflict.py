#!/usr/bin/env python3
"""Descriptive supported split agreement/conflict on provisional full-guide edges."""
import argparse
import json
from pathlib import Path
from collections import defaultdict,Counter
from Bio import Phylo
from audit_joint_path_uncertainty import checked,rows,sha
from run_paired_marker_fits import tree_edges
from prepare_paired_phylogenetic_inputs import write_table


def canonical(side,universe):
    return min((tuple(sorted(side)),tuple(sorted(universe-side))),key=lambda x:(len(x),x))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['snapshot','readback','trees','output']:ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--guides',nargs=2,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    sr=checked(a.snapshot);ar=checked(a.readback)
    if ar['status']!='passed_all_snapshot_input_and_graph_split_readbacks' or ar['snapshot_receipt_sha256']!=sha(a.snapshot/'receipt.json'):raise ValueError('Snapshot readback mismatch')
    pins={str(f/'receipt.json'):sha(f/'receipt.json') for f in [a.snapshot,a.readback]+a.guides}
    guides={};full=None
    for g in a.guides:
        gr=checked(g)
        if gr['status']!='passed_full_species_guide_readback':raise ValueError('Audited guide required')
        taxa={t.name for t in Phylo.read(g/'guide.treefile','newick').get_terminals()}
        if full is not None and full!=taxa:raise ValueError('Guide universes differ')
        full=taxa;guides[g.name]={s:l for s,l in tree_edges(g/'guide.treefile',taxa).items() if len(s)>1}
    support=defaultdict(dict)
    for r in rows(a.snapshot/'branch_support.tsv'):
        support[r['marker']][tuple(r['smaller_side_taxa'].split(';'))]=None if r['sh_alrt_percent']=='' else float(r['sh_alrt_percent'])
    output=[]
    for source in sr['inputs']:
        marker=source['marker'];tree=a.trees/marker/'tree.treefile'
        if sha(tree)!=source['tree_sha256']:raise ValueError('Changed marker tree')
        pins[str(tree)]=sha(tree);taxa={t.name for t in Phylo.read(tree,'newick').get_terminals()}
        if not taxa<=full:raise ValueError('Unknown taxa')
        index={t:i for i,t in enumerate(sorted(taxa))};universe=(1<<len(taxa))-1
        def bits(side):return sum(1<<index[t] for t in side)
        supported=[(bits(side),score,side) for side,score in support[marker].items() if score is not None]
        for guide,edges in guides.items():
            for fullside in edges:
                restricted=set(fullside)&taxa;side=canonical(restricted,taxa)
                eligible=len(side)>=2
                exact=support[marker].get(side)
                strongest=None;witness=''
                if eligible:
                    x=bits(side);nx=universe^x
                    for y,score,other in supported:
                        ny=universe^y
                        if x&y and x&ny and nx&y and nx&ny:
                            if strongest is None or score>strongest:
                                strongest=score;witness=';'.join(other)
                for cutoff in [80,95]:
                    status=('uninformative_taxon_coverage' if not eligible else 'supported_concordance' if exact is not None and exact>=cutoff else 'supported_conflict' if strongest is not None and strongest>=cutoff else 'unresolved_at_support_cutoff')
                    output.append(dict(marker=marker,guide=guide,full_split_taxa=';'.join(fullside),marker_taxa=len(taxa),restricted_split_taxa=';'.join(side),sh_alrt_cutoff=cutoff,status=status,exact_split_sh_alrt=exact if exact is not None else '',maximum_conflicting_sh_alrt=strongest if strongest is not None else '',maximum_conflict_witness_split=witness))
        print(marker,'completed',flush=True)
    grouped=defaultdict(Counter)
    for r in output:grouped[r['guide'],r['full_split_taxa'],r['sh_alrt_cutoff']][r['status']]+=1
    summaries=[]
    states=['uninformative_taxon_coverage','supported_concordance','supported_conflict','unresolved_at_support_cutoff']
    for (guide,split,cutoff),counts in sorted(grouped.items()):
        summaries.append(dict(guide=guide,full_split_taxa=split,sh_alrt_cutoff=cutoff,completed_markers=sr['completed_markers'],**{s:counts[s] for s in states}))
    a.output.mkdir(parents=True)
    write_table(a.output/'marker_guide_conflict.tsv',output);write_table(a.output/'edge_summary.tsv',summaries)
    result=dict(status='complete_provisional_supported_split_diagnostic',completed_markers=sr['completed_markers'],planned_markers=sr['planned_markers'],pending_markers=sr.get('pending_markers',[]),rows=len(output),edge_summary_rows=len(summaries),source_pins=pins,script_sha256=sha(Path(__file__)),interpretation='SH-aLRT cutoffs 80 and 95 are descriptive sensitivity settings, not bootstrap or posterior probabilities. Conflict requires all four bipartition intersections nonempty. Missing or low support is unresolved. Restricted splits can represent collapsed guide paths, so rows are not independent branch replicates. ' + ('Full planned marker batch. ' if sr['completed_markers'] == sr['planned_markers'] and not sr.get('pending_markers') else 'Incomplete completion-order-biased marker snapshot. ') + 'Not a gene concordance factor, quartet statistic, supported species topology, reconciliation or explanation of biological discordance.',artifacts={p.name:sha(p) for p in a.output.iterdir()})
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_pins','pending_markers']},indent=2))

if __name__=='__main__':main()
