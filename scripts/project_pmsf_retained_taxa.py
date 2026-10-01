#!/usr/bin/env python3
"""Recompute full-baseline bootstrap support on four exact retained taxon sets.

Projection is a matched-taxon reference diagnostic, not a native subset refit.
Collapsed edges contribute their summed original lengths; SH-aLRT is unavailable.
Each projected bootstrap split is counted once per raw replicate, including when
multiple original splits collapse to it. No original support label is inherited.
"""
import argparse
from collections import Counter,defaultdict
import csv
import json
import math
from pathlib import Path
from Bio import Phylo
from retained_taxon_projection_sources import load,EDGE_FIELDS,BOOT_FIELDS,BOUNDARY_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def raw_edges(tree,universe):
    names=[n.name for n in tree.get_terminals()]
    assert len(names)==len(universe) and set(names)==universe
    descendants={};result=[]
    for node in tree.find_clades(order='postorder'):
        side={node.name} if node.is_terminal() else set().union(*(descendants[c] for c in node.clades))
        descendants[node]=side
        if node is tree.root:
            assert node.branch_length in [None,0];continue
        length=node.branch_length;assert length is not None and math.isfinite(length) and length>=0
        result.append((side,length))
    assert len(result)==2*len(universe)-3
    return result


def project(edges,retained):
    lengths=defaultdict(list)
    for side,length in edges:
        left=tuple(sorted(side & retained));right=tuple(sorted(retained.difference(side)))
        if not left or not right:continue
        key=min(left,right,key=lambda x:(len(x),x));lengths[key].append(length)
    assert len(lengths)==2*len(retained)-3
    return {key:(math.fsum(values),len(values)) for key,values in lengths.items()}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());sources,universe,policies,bindings=load(plan,args.plan)
    output=Path(plan['output']);output.mkdir(parents=True,exist_ok=False)
    edge_rows=[];boot_rows=[];boundary_rows=[];summary=[]
    for label,source in sources.items():
        root=Path(source['spec']['run']);counts={p:Counter() for p in policies};replicates=0;collapsed_counts=Counter()
        for tree in Phylo.parse(root/'pmsf.ufboot','newick'):
            raw=raw_edges(tree,universe)
            for policy,spec in policies.items():
                projected=project(raw,spec['taxa']);counts[policy].update(projected.keys())
                collapsed_counts[policy]+=sum(n-1 for length,n in projected.values())
            replicates+=1
        assert replicates==1000
        for policy,spec in policies.items():
            n=len(spec['taxa']);outgroups=spec['outgroups'];key=min(tuple(sorted(outgroups)),tuple(sorted(spec['taxa']-outgroups)),key=lambda x:(len(x),x))
            for split,count in sorted(counts[policy].items()):
                assert 0<count<=1000
                boot_rows.append(dict(baseline_run=label,policy=policy,split_taxa_json=json.dumps(split),
                    replicates_with_projected_split=count,empirical_projected_ufboot_percent=count/10))
            merged=0
            for kind,file in [('ml','pmsf.treefile'),('consensus','pmsf.contree')]:
                projected=project(raw_edges(Phylo.read(root/file,'newick'),universe),spec['taxa'])
                for split,(length,components) in sorted(projected.items()):
                    edge_rows.append(dict(baseline_run=label,policy=policy,tree_type=kind,split_taxa_json=json.dumps(split),
                        projected_branch_length_sum=length,original_edge_components=components,
                        empirical_projected_ufboot_percent=counts[policy][split]/10,sh_alrt_percent=None))
                    merged+=components-1
                boundary_rows.append(dict(baseline_run=label,policy=policy,tree_type=kind,retained_taxa=n,
                    retained_ingroup=spec['roles']['ingroup'],retained_outgroup=spec['roles']['outgroup'],
                    boundary_taxa_json=json.dumps(key),boundary_present=key in projected,
                    empirical_projected_ufboot_percent=counts[policy][key]/10,sh_alrt_percent=None))
            summary.append(dict(baseline_run=label,policy=policy,retained_taxa=n,roles=spec['roles'],
                projected_bootstrap_trees=1000,projected_empirical_splits=len(counts[policy]),
                collapsed_bootstrap_edge_components=collapsed_counts[policy],collapsed_ML_consensus_edge_components=merged))
        print('projected_all_original_bootstraps_and_views',label,len(summary),'/16',flush=True)
    assert len(summary)==16 and len(boundary_rows)==32
    for name,rows,fields in [('projected_tree_edges.tsv',edge_rows,EDGE_FIELDS),('projected_bootstrap_splits.tsv',boot_rows,BOOT_FIELDS),('role_boundary.tsv',boundary_rows,BOUNDARY_FIELDS)]:
        with (output/name).open('x') as handle:
            writer=csv.DictWriter(handle,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    verify(bindings)
    receipt=dict(status='complete_full_retained_taxon_baseline_projections_pending_independent_readback',plan_sha256=sha(args.plan),
        original_baseline_runs=4,original_raw_bootstrap_trees=4000,policy_baseline_combinations=16,
        projected_bootstrap_tree_states=16000,tree_views=32,edge_rows=len(edge_rows),bootstrap_split_rows=len(boot_rows),
        role_boundary_rows=len(boundary_rows),boundary_views_with_role_split=sum(r['boundary_present'] for r in boundary_rows),
        summaries=summary,source_hashes=bindings,artifacts={name:sha(output/name) for name in
        ['projected_tree_edges.tsv','projected_bootstrap_splits.tsv','role_boundary.tsv']},scientific_eligibility=False,scope=plan['scope'])
    with (output/'receipt.json').open('x') as handle:handle.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_hashes','artifacts','summaries','scope']},indent=2),flush=True)


if __name__=='__main__':main()
