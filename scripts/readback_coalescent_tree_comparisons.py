#!/usr/bin/env python3
"""Independent raw-tree split, RF, support-disposition and quartet-witness reader."""
import argparse
import csv
import gzip
import itertools
import json
import math
from pathlib import Path

import dendropy
from dendropy.calculate import treecompare

from coalescent_tree_comparison_sources import (load,METRICS,PAIR_FIELDS,PRESENCE_FIELDS,CONFLICT_FIELDS,BOUNDARY_FIELDS,SUMMARY_FIELDS)
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def read(path,fields):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as handle:
        reader=csv.DictReader(handle,delimiter='\t');assert reader.fieldnames==fields
        return list(reader)


def numeric(value,expected):
    if expected is None:assert value==''
    elif isinstance(expected,str):assert value==expected
    else:assert math.isclose(float(value),expected,rel_tol=1e-12,abs_tol=1e-12)


def supported(row,rule):
    # Deliberately separate the producer's rule function.
    if rule=='local_PP95':return row['local_posterior']>=.95
    if rule=='original_SH80_and_empirical_UFB95':return row['sh_alrt_percent']>=80 and row['empirical_ufb_percent']>=95
    assert rule in ['original_consensus_empirical_UFB95_SH_unavailable','projected_empirical_UFB95_SH_unavailable']
    assert row['sh_alrt_percent'] is None
    return row['empirical_ufb_percent']>=95


def inspect(groups,universe,root):
    root=Path(root)
    pairs=read(root/'comparisons.tsv',PAIR_FIELDS);presence=read(root/'split_presence.tsv.gz',PRESENCE_FIELDS)
    conflicts=read(root/'conflicts.tsv.gz',CONFLICT_FIELDS);boundaries=read(root/'role_boundaries.tsv',BOUNDARY_FIELDS)
    pi={(r['cohort'],r['view_a'],r['view_b']):r for r in pairs}
    si={(r['cohort'],r['view'],tuple(json.loads(r['split_taxa_json']))):r for r in presence}
    ci={(r['cohort'],r['view_a'],r['view_b'],tuple(json.loads(r['split_a_taxa_json'])),tuple(json.loads(r['split_b_taxa_json']))):r for r in conflicts}
    bi={(r['cohort'],r['view']):r for r in boundaries}
    assert all(len(index)==len(rows) for index,rows in [(pi,pairs),(si,presence),(ci,conflicts),(bi,boundaries)])
    pair_seen=set();presence_seen=set();conflict_seen=set();boundary_seen=set();summaries=[]
    for cohort,group in sorted(groups.items()):
        taxa=group['taxa'];views=group['views'];ns=dendropy.TaxonNamespace(sorted(universe));trees={};nodes={}
        def canon(side):
            a,b=tuple(sorted(side)),tuple(sorted(taxa-set(side)))
            return a if (len(a),a)<=(len(b),b) else b
        for name,view in sorted(views.items()):
            tree=dendropy.Tree.get(path=view['source_tree'],schema='newick',rooting='force-unrooted',preserve_underscores=True,taxon_namespace=ns)
            original={n.taxon.label for n in tree.leaf_node_iter()}
            assert original==(universe if view['prune'] else taxa)
            if view['prune']:tree.retain_taxa_with_labels(sorted(taxa))
            tree.deroot();tree.encode_bipartitions()
            assert {n.taxon.label for n in tree.leaf_node_iter()}==taxa
            lookup={}
            for node in tree.preorder_node_iter():
                if node is tree.seed_node or node.is_leaf():continue
                key=canon({n.taxon.label for n in node.leaf_iter()});assert key not in lookup
                lookup[key]=node
                metric=view['splits'][key]
                assert node.edge.length is not None and math.isclose(node.edge.length,metric['branch_length'],rel_tol=1e-12,abs_tol=1e-12)
                if view['family']=='coalescent':
                    assert node.label.startswith('[') and node.label.endswith(']')
                    native=dict(field.split('=') for field in node.label[1:-1].split(';'))
                    assert float(native['pp1'])==metric['local_posterior'] and float(native['EN'])==metric['effective_genes']
            assert set(lookup)==set(view['splits']) and len(lookup)==len(taxa)-3
            trees[name]=tree;nodes[name]=lookup
        union=set().union(*(set(v) for v in nodes.values()))
        for name,view in views.items():
            for split in union:
                identity=cohort,name,split;row=si[identity];presence_seen.add(identity);source=view['splits'].get(split)
                assert row['present']==str(source is not None) and row['family']==view['family'] and row['support_rule']==view['support_rule']
                for field in METRICS:numeric(row[field],source[field] if source else None)
                assert row['support_rule_met']==(str(supported(source,view['support_rule'])) if source else '')
            boundary=canon(group['outgroups']);r=bi[cohort,name];boundary_seen.add((cohort,name));source=view['splits'].get(boundary)
            assert tuple(json.loads(r['boundary_taxa_json']))==boundary and r['boundary_present']==str(source is not None)
            assert int(r['ingroup'])==group['roles']['ingroup'] and int(r['outgroup'])==group['roles']['outgroup']
            assert r['support_rule']==view['support_rule']
            for field in METRICS:numeric(r[field],source[field] if source else None)
            assert r['support_rule_met']==(str(supported(source,view['support_rule'])) if source else '')
        grid=[(a,b) for a,b in itertools.combinations(sorted(views),2) if
              views[a]['family']=='coalescent' or views[b]['family']=='coalescent']
        for a,b in grid:
            identity=cohort,a,b;row=pi[identity];pair_seen.add(identity);left,right=nodes[a],nodes[b];common=set(left)&set(right)
            count=high=0
            for x,y in itertools.product(set(left)-common,set(right)-common):
                if left[x].edge.bipartition.is_compatible_with(right[y].edge.bipartition):continue
                cid=cohort,a,b,x,y;c=ci[cid];conflict_seen.add(cid);count+=1
                sa=supported(views[a]['splits'][x],views[a]['support_rule']);sb=supported(views[b]['splits'][y],views[b]['support_rule']);high+=sa and sb
                assert c['support_rule_met_a']==str(sa) and c['support_rule_met_b']==str(sb) and c['both_support_rules_met']==str(sa and sb)
                assert c['support_rule_a']==views[a]['support_rule'] and c['support_rule_b']==views[b]['support_rule']
                q=c['witness_quartet'].split(';');assert len(q)==len(set(q))==4 and set(q)<=taxa
                assert {(t in x,t in y) for t in q}=={(False,False),(False,True),(True,False),(True,True)}
                xs,ys=set(x),set(y);assert q==[min(xs&ys),min(xs-ys),min(ys-xs),min(taxa-(xs|ys))]
            distance=treecompare.symmetric_difference(trees[a],trees[b],is_bipartitions_updated=True)
            assert distance==int(row['rf_distance'])==len(left)+len(right)-2*len(common)
            assert int(row['shared_splits'])==len(common) and int(row['unique_splits_a'])==len(left)-len(common) and int(row['unique_splits_b'])==len(right)-len(common)
            assert math.isclose(float(row['normalized_rf']),distance/(len(left)+len(right)),rel_tol=1e-12,abs_tol=1e-12)
            assert int(row['incompatible_pairs'])==count and int(row['both_support_rules_met_pairs'])==high
            assert row['support_rule_a']==views[a]['support_rule'] and row['support_rule_b']==views[b]['support_rule']
            assert row['comparison_scope']==('coalescent_vs_concatenated' if views[a]['family']!=views[b]['family'] else 'coalescent_alignment_support_sensitivity')
        coal=[set(nodes[n]) for n,v in views.items() if v['family']=='coalescent'];refs=[set(nodes[n]) for n,v in views.items() if v['family']=='concatenated']
        summaries.append(dict(cohort=cohort,taxa=len(taxa),roles=group['roles'],coalescent_views=len(coal),concatenated_views=len(refs),
            shared_all_coalescent=len(set.intersection(*coal)),shared_all_references=len(set.intersection(*refs)),shared_all_views=len(set.intersection(*(coal+refs)))))
        print('independent_raw_coalescent_reference_cohort_verified',cohort,len(grid),flush=True)
    assert pair_seen==set(pi) and presence_seen==set(si) and conflict_seen==set(ci) and boundary_seen==set(bi)
    return dict(cohorts=len(groups),tree_views=sum(len(g['views']) for g in groups.values()),comparison_rows=len(pairs),
        split_presence_rows=len(presence),incompatible_pairs=len(conflicts),role_boundary_rows=len(boundaries),
        boundary_views_with_role_split=sum(r['boundary_present']=='True' for r in boundaries),cohort_summaries=summaries)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());groups,universe,bindings=load(plan,args.plan)
    root=Path(plan['output']);rp=root/'receipt.json';producer=json.loads(rp.read_text());rh=sha(rp)
    assert producer['status']=='complete_coalescent_reference_tree_comparisons_pending_independent_readback'
    assert producer['plan_sha256']==sha(args.plan) and producer['source_hashes']==bindings
    for name,digest in producer['artifacts'].items():assert sha(root/name)==digest
    summary=inspect(groups,universe,root);assert all(producer[k]==summary[k] for k in SUMMARY_FIELDS)
    verify(bindings);assert sha(rp)==rh
    for name,digest in producer['artifacts'].items():assert sha(root/name)==digest
    result=dict(status='passed_coalescent_reference_full_raw_tree_comparison_readback',plan_sha256=sha(args.plan),
        producer_receipt_sha256=rh,**summary,source_hashes=bindings,dendropy_version=dendropy.__version__,scientific_eligibility=False,scope=plan['scope'])
    with args.output.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if not isinstance(v,list)}),flush=True)


if __name__=='__main__':main()
