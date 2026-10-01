#!/usr/bin/env python3
"""Independently reconstruct all retained-reference RF/conflict/presence exports."""
import argparse
import csv
import gzip
import itertools
import json
import math
from pathlib import Path
import dendropy
from dendropy.calculate import treecompare
from retained_tree_comparison_sources import load,PAIR_FIELDS,PRESENCE_FIELDS,CONFLICT_FIELDS,RULE
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def read(path,fields,compressed=False):
    opener=gzip.open if compressed else open
    with opener(path,'rt') as handle:
        reader=csv.DictReader(handle,delimiter='\t');assert reader.fieldnames==fields
        return list(reader)


def inspect(groups,universe,root):
    root=Path(root);pairs=read(root/'comparisons.tsv',PAIR_FIELDS);presence=read(root/'split_presence.tsv',PRESENCE_FIELDS)
    conflicts=read(root/'conflicts.tsv.gz',CONFLICT_FIELDS,True)
    pi={(r['policy'],r['view_a'],r['view_b']):r for r in pairs};assert len(pi)==len(pairs)
    si={(r['policy'],r['view'],tuple(json.loads(r['split_taxa_json']))):r for r in presence};assert len(si)==len(presence)
    ci={(r['policy'],r['view_a'],r['view_b'],tuple(json.loads(r['split_a_taxa_json'])),tuple(json.loads(r['split_b_taxa_json']))):r for r in conflicts};assert len(ci)==len(conflicts)
    pair_seen=set();presence_seen=set();conflict_seen=set();policy_summaries=[]
    for policy,group in groups.items():
        retained=group['taxa'];ns=dendropy.TaxonNamespace(sorted(universe));trees={};nodes={};views=group['views']
        def canonical(side):
            left,right=tuple(sorted(side)),tuple(sorted(retained.difference(side)))
            return left if (len(left),left)<=(len(right),right) else right
        for name,view in views.items():
            file='pmsf.treefile' if view['kind']=='ml' else 'pmsf.contree'
            tree=dendropy.Tree.get(path=str(Path(view['source']['run'])/file),schema='newick',rooting='force-unrooted',preserve_underscores=True,taxon_namespace=ns)
            assert len(list(tree.leaf_node_iter()))==len(universe) and {n.taxon.label for n in tree.leaf_node_iter()}==universe
            tree.retain_taxa_with_labels(sorted(retained));tree.deroot();tree.encode_bipartitions()
            assert len(list(tree.leaf_node_iter()))==len(retained) and {n.taxon.label for n in tree.leaf_node_iter()}==retained
            lookup={}
            for node in tree.preorder_node_iter():
                if node is tree.seed_node or node.is_leaf():continue
                key=canonical({n.taxon.label for n in node.leaf_iter()});assert key not in lookup;lookup[key]=node
                row=view['splits'][key]
                assert math.isclose(float(row['projected_branch_length_sum']),node.edge.length,rel_tol=1e-12,abs_tol=1e-12)
                assert row['sh_alrt_percent']==''
            assert set(lookup)==set(view['splits']) and len(lookup)==len(retained)-3
            trees[name]=tree;nodes[name]=lookup
        union=set().union(*(set(v) for v in nodes.values()))
        for name,lookup in nodes.items():
            for split in union:
                identity=policy,name,split;row=si[identity];presence_seen.add(identity)
                source=views[name]['splits'].get(split);assert row['present']==str(source is not None) and row['sh_alrt_percent']==''
                for field in ['projected_branch_length_sum','empirical_projected_ufboot_percent']:
                    assert row[field]==(source[field] if source else '')
        pair_grid=[(a,b) for a,b in itertools.combinations(views,2) if views[a]['kind']==views[b]['kind'] or views[a]['run']==views[b]['run']]
        assert len(pair_grid)==16
        for a,b in pair_grid:
            identity=policy,a,b;row=pi[identity];pair_seen.add(identity);left,right=nodes[a],nodes[b]
            common=set(left)&set(right);nconf=high=0
            for x,y in itertools.product(set(left)-common,set(right)-common):
                if left[x].edge.bipartition.is_compatible_with(right[y].edge.bipartition):continue
                cid=policy,a,b,x,y;c=ci[cid];conflict_seen.add(cid);nconf+=1
                sa=float(views[a]['splits'][x]['empirical_projected_ufboot_percent']);sb=float(views[b]['splits'][y]['empirical_projected_ufboot_percent'])
                strong=sa>=95 and sb>=95;high+=strong
                assert float(c['empirical_projected_ufboot_percent_a'])==sa and float(c['empirical_projected_ufboot_percent_b'])==sb and c['both_UFB95']==str(strong)
                quartet=c['witness_quartet'].split(';');assert len(set(quartet))==4 and set(quartet)<=retained
                assert {(t in x,t in y) for t in quartet}=={(False,False),(False,True),(True,False),(True,True)}
                xs,ys=set(x),set(y)
                assert quartet==[min(xs&ys),min(xs-ys),min(ys-xs),min(retained-(xs|ys))]
            distance=treecompare.symmetric_difference(trees[a],trees[b],is_bipartitions_updated=True)
            assert distance==int(row['rf_distance'])==len(left)+len(right)-2*len(common)
            assert row['shared_internal_splits']==str(len(common)) and row['unique_internal_splits_a']==str(len(left)-len(common)) and row['unique_internal_splits_b']==str(len(right)-len(common))
            assert math.isclose(float(row['normalized_rf']),distance/(len(left)+len(right)),rel_tol=1e-12,abs_tol=1e-12)
            assert row['incompatible_split_pairs']==str(nconf) and row['both_UFB95_incompatible_pairs']==str(high) and row['support_rule']==RULE
            assert row['comparison_scope']==('within_baseline_ML_vs_consensus' if views[a]['run']==views[b]['run'] else 'cross_baseline_'+views[a]['kind'])
        policy_summaries.append(dict(policy=policy,retained_taxa=len(retained),roles=group['roles'],views=8,comparisons=16,
            shared_all_ml=len(set.intersection(*(set(nodes[n]) for n,v in views.items() if v['kind']=='ml'))),
            shared_all_consensus=len(set.intersection(*(set(nodes[n]) for n,v in views.items() if v['kind']=='consensus'))),
            shared_all_eight=len(set.intersection(*(set(v) for v in nodes.values())))))
        print('independently_checked_every_retained_comparison',policy,flush=True)
    assert pair_seen==set(pi) and presence_seen==set(si) and conflict_seen==set(ci)
    return dict(policies=len(groups),tree_views=sum(len(g['views']) for g in groups.values()),comparison_rows=len(pairs),
        split_presence_rows=len(presence),incompatible_pairs=len(conflicts),policy_summaries=policy_summaries)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());groups,universe,bindings=load(plan,args.plan);root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    assert r['status']=='complete_full_retained_baseline_tree_comparisons_pending_readback' and r['plan_sha256']==sha(args.plan) and r['source_hashes']==bindings
    for name,digest in r['artifacts'].items():assert sha(root/name)==digest
    summary=inspect(groups,universe,root);assert all(r[k]==v for k,v in summary.items())
    assert summary['policies']==4 and summary['tree_views']==32 and summary['comparison_rows']==64
    verify(bindings);assert sha(rp)==rh
    for name,digest in r['artifacts'].items():assert sha(root/name)==digest
    proof=dict(status='passed_full_retained_baseline_tree_comparison_independent_readback',plan_sha256=sha(args.plan),
        producer_receipt_sha256=rh,**summary,source_hashes=bindings,dendropy_version=dendropy.__version__,scientific_eligibility=False,
        scope='RawDendroPy baseline tree parsing andexactretainedtaxon deletion reconstruct everyinternal split/pathlength inall32views. DendroPy RF/bipartitioncompatibility reconstructall64pairs/allincompatible unique splitpairs, quartetwitnesses,allpresencecells andeachprojectedUFB95screen. Complete closedprojection raw-bootstrap sourceproof bound. SH-aLRTunavailable/noformal hypothesiscalibration/native-refit/root/framework/biologicalaimcompletion.')
    with args.output.open('x') as handle:handle.write(json.dumps(proof,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='policy_summaries'}),flush=True)


if __name__=='__main__':main()
