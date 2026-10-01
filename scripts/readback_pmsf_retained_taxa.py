#!/usr/bin/env python3
"""Reconstruct complete baseline projections through independent tree pruning.

DendroPy deletes taxa and suppresses paths; a second unit-edge copy checks path
component counts. All projected bootstrap frequencies are recomputed from raw
trees. No producer projection or support algorithm is imported.
"""
import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import dendropy
from retained_taxon_projection_sources import load,EDGE_FIELDS,BOOT_FIELDS,BOUNDARY_FIELDS
from reference_measurement_union_sources import verify,bind
from run_ortholog_pair_guide_comparison import sha


def raw_tree(tree,universe):
    names=[n.taxon.label for n in tree.leaf_node_iter()]
    assert len(names)==len(universe) and set(names)==universe
    assert tree.seed_node.edge.length in [None,0]
    edges=[n.edge for n in tree.preorder_node_iter() if n is not tree.seed_node]
    assert len(edges)==2*len(universe)-3
    assert all(e.length is not None and math.isfinite(e.length) and e.length>=0 for e in edges)


def prune(tree,retained,unit_edges=False):
    projected=tree.clone(depth=1)
    if unit_edges:
        for node in projected.preorder_node_iter():
            node.edge.length=None if node is projected.seed_node else 1.0
    projected.retain_taxa_with_labels(sorted(retained))
    projected.deroot()
    leaves=list(projected.leaf_node_iter())
    assert len(leaves)==len(retained) and {n.taxon.label for n in leaves}==retained
    result={}
    for node in projected.preorder_node_iter():
        if node is projected.seed_node:continue
        side={n.taxon.label for n in node.leaf_iter()}
        a,b=tuple(sorted(side)),tuple(sorted(retained.difference(side)))
        key=a if (len(a),a)<=(len(b),b) else b
        assert key and key not in result and node.edge.length is not None and node.edge.length>=0
        result[key]=node.edge.length
    assert len(result)==2*len(retained)-3
    return result


def table(path,fields,keys):
    with Path(path).open() as handle:
        reader=csv.DictReader(handle,delimiter='\t');assert reader.fieldnames==fields
        rows={}
        for row in reader:
            identity=tuple(row[k] for k in keys)
            if 'split_taxa_json' in fields:
                side=json.loads(row['split_taxa_json']);assert side and len(side)==len(set(side))
                identity+=tuple([tuple(side)])
            assert identity not in rows;rows[identity]=row
    return rows


def inspect_tables(sources,universe,policies,output):
    output=Path(output)
    edges=table(output/'projected_tree_edges.tsv',EDGE_FIELDS,['baseline_run','policy','tree_type'])
    boots=table(output/'projected_bootstrap_splits.tsv',BOOT_FIELDS,['baseline_run','policy'])
    boundary=table(output/'role_boundary.tsv',BOUNDARY_FIELDS,['baseline_run','policy','tree_type'])
    edge_seen=set();boot_seen=set();boundary_seen=set();summaries=[]
    for label,source in sources.items():
        root=Path(source['spec']['run']);counts={p:Counter() for p in policies};collapsed=Counter();replicates=0
        for tree in dendropy.Tree.yield_from_files([str(root/'pmsf.ufboot')],schema='newick',rooting='force-unrooted',preserve_underscores=True):
            raw_tree(tree,universe)
            for policy,spec in policies.items():
                projected=prune(tree,spec['taxa'],unit_edges=True)
                assert all(v==int(v) and v>=1 for v in projected.values())
                counts[policy].update(projected.keys());collapsed[policy]+=int(sum(projected.values()))-len(projected)
            replicates+=1
        assert replicates==1000
        raw={kind:dendropy.Tree.get(path=str(root/file),schema='newick',rooting='force-unrooted',preserve_underscores=True)
             for kind,file in [('ml','pmsf.treefile'),('consensus','pmsf.contree')]}
        for tree in raw.values():raw_tree(tree,universe)
        for policy,spec in policies.items():
            retained=spec['taxa'];merged=0
            for split,count in counts[policy].items():
                identity=label,policy,split;row=boots[identity];boot_seen.add(identity)
                assert row['replicates_with_projected_split']==str(count) and 0<count<=1000
                assert float(row['empirical_projected_ufboot_percent'])==count/10
            role_side=spec['outgroups'];other=retained.difference(role_side)
            role_key=min(tuple(sorted(role_side)),tuple(sorted(other)),key=lambda s:(len(s),s))
            for kind,tree in raw.items():
                projected=prune(tree,retained);components=prune(tree,retained,unit_edges=True)
                assert set(projected)==set(components)
                merged+=int(sum(components.values()))-len(components)
                for split,length in projected.items():
                    identity=label,policy,kind,split;row=edges[identity];edge_seen.add(identity)
                    assert math.isclose(float(row['projected_branch_length_sum']),length,rel_tol=1e-12,abs_tol=1e-12)
                    assert components[split]==int(components[split]) and row['original_edge_components']==str(int(components[split]))
                    assert row['sh_alrt_percent']=='' and float(row['empirical_projected_ufboot_percent'])==counts[policy][split]/10
                identity=label,policy,kind;row=boundary[identity];boundary_seen.add(identity)
                assert row['retained_taxa']==str(len(retained)) and row['retained_ingroup']==str(spec['roles']['ingroup']) and row['retained_outgroup']==str(spec['roles']['outgroup'])
                assert tuple(json.loads(row['boundary_taxa_json']))==role_key and row['boundary_present']==str(role_key in projected)
                assert row['sh_alrt_percent']=='' and float(row['empirical_projected_ufboot_percent'])==counts[policy][role_key]/10
            summaries.append(dict(baseline_run=label,policy=policy,retained_taxa=len(retained),roles=spec['roles'],
                projected_bootstrap_trees=1000,projected_empirical_splits=len(counts[policy]),
                collapsed_bootstrap_edge_components=collapsed[policy],collapsed_ML_consensus_edge_components=merged))
        print('independent_pruning_readback_complete',label,len(summaries),flush=True)
    assert set(edges)==edge_seen and set(boots)==boot_seen and set(boundary)==boundary_seen
    return dict(original_baseline_runs=len(sources),original_raw_bootstrap_trees=1000*len(sources),
        policy_baseline_combinations=len(summaries),projected_bootstrap_tree_states=1000*len(summaries),
        tree_views=2*len(summaries),edge_rows=len(edges),bootstrap_split_rows=len(boots),role_boundary_rows=len(boundary),
        boundary_views_with_role_split=sum(row['boundary_present']=='True' for row in boundary.values()),summaries=summaries)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());sources,universe,policies,bindings=load(plan,args.plan)
    root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    assert r['status']=='complete_full_retained_taxon_baseline_projections_pending_independent_readback' and r['plan_sha256']==sha(args.plan) and r['source_hashes']==bindings
    for name,digest in r['artifacts'].items():assert sha(root/name)==digest
    summary=inspect_tables(sources,universe,policies,root)
    assert all(r[k]==v for k,v in summary.items())
    assert summary['original_baseline_runs']==4 and summary['policy_baseline_combinations']==16 and summary['tree_views']==32
    verify(bindings);assert sha(rp)==rh
    for name,digest in r['artifacts'].items():assert sha(root/name)==digest
    bind(bindings,rp)
    proof=dict(status='passed_full_retained_taxon_baseline_projection_independent_readback',plan_sha256=sha(args.plan),
        producer_receipt_sha256=rh,**summary,source_hashes=bindings,dendropy_version=dendropy.__version__,scientific_eligibility=False,
        scope='DendroPy taxon deletion/path suppression, separate unit-edge path component reconstruction, all4000raw baseline replicates projected onto allfour exactpolicysets, all16000bootstrap states/32views andeveryedge/frequency/boundary field reconstructed. No Bio.Phylo projection or support algorithm imported. Recomputed projected UFB is conditional on original526tip inference; no inherited SH-aLRT/native subset refit/root/modeladequacy/accepted framework claim.')
    with args.output.open('x') as handle:handle.write(json.dumps(proof,indent=2)+'\n')
    print(json.dumps({k:v for k,v in proof.items() if k not in ['source_hashes','summaries','scope']},indent=2),flush=True)


if __name__=='__main__':main()
