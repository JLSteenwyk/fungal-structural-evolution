#!/usr/bin/env python3
"""Reconstruct PMSF split grids and incompatibility with DendroPy."""
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import dendropy


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    root=Path('results/phylogeny/pmsf-three-run-topology-sensitivity-20260927-v1')
    r=json.loads((root/'receipt.json').read_text())
    for p,d in r['sources'].items():assert sha(p)==d
    for n,d in r['artifacts'].items():assert sha(root/n)==d
    ns=dendropy.TaxonNamespace();trees={};edges={};supports={}
    for name,(run,audit) in r['runs'].items():
        tree=dendropy.Tree.get(path=str(Path('results/phylogeny')/run/'pmsf.treefile'),schema='newick',rooting='force-unrooted',preserve_underscores=True,taxon_namespace=ns)
        taxa={n.taxon.label for n in tree.leaf_node_iter()};assert len(taxa)==526
        tree.encode_bipartitions();lookup={}
        for n in tree.preorder_node_iter():
            if n is tree.seed_node or n.is_leaf():continue
            side={x.taxon.label for x in n.leaf_iter()}
            key=min(tuple(sorted(side)),tuple(sorted(taxa-side)),key=lambda x:(len(x),x))
            assert key not in lookup;lookup[key]=n
        assert len(lookup)==523
        edges[name]=lookup;trees[name]=tree
        supports[name]={}
        for row in read(Path('results/phylogeny')/audit/'branch_support.tsv'):
            if row['tree']!='ml':continue
            side=set(json.loads(row['split_taxa_json']));key=min(tuple(sorted(side)),tuple(sorted(taxa-side)),key=lambda x:(len(x),x))
            node=lookup[key];alrt,ufb=map(float,node.label.split('/'))
            assert alrt==float(row['sh_alrt_percent']) and abs(ufb-float(row['empirical_ufboot_percent']))<=.500001
            assert node.edge.length==float(row['branch_length'])
            supports[name][key]=row
    union=set().union(*(set(s) for s in edges.values()))
    presence=read(root/'split_presence.tsv');assert len(presence)==len(union)*3
    seen=set()
    for row in presence:
        name=row['run'];key=tuple(json.loads(row['split_taxa_json']));assert (name,key) not in seen;seen.add((name,key))
        assert key in union and (row['present']=='True')==(key in edges[name])
        if key in edges[name]:
            source=supports[name][key]
            for col,src in [('branch_length','branch_length'),('sh_alrt','sh_alrt_percent'),('empirical_ufb','empirical_ufboot_percent')]:assert float(row[col])==float(source[src])
        else:assert row['branch_length']==row['sh_alrt']==row['empirical_ufb']==''
    actual=read(root/'conflicts.tsv');indexed={(x['run_a'],x['run_b'],tuple(json.loads(x['split_a_taxa_json'])),tuple(json.loads(x['split_b_taxa_json']))):x for x in actual};assert len(indexed)==len(actual)
    expected=set();summaries=read(root/'comparisons.tsv');assert len(summaries)==3
    for row in summaries:
        a,b=row['run_a'],row['run_b'];left,right=edges[a],edges[b];common=set(left)&set(right);n=high=0
        for x,y in itertools.product(set(left)-common,set(right)-common):
            if left[x].bipartition.is_compatible_with(right[y].bipartition):continue
            key=(a,b,x,y);expected.add(key);entry=indexed[key];n+=1
            sx,sy=supports[a][x],supports[b][y]
            for label,source in [('a',sx),('b',sy)]:
                assert float(entry['sh_alrt_'+label])==float(source['sh_alrt_percent'])
                assert float(entry['empirical_ufb_'+label])==float(source['empirical_ufboot_percent'])
            both=all(float(s['sh_alrt_percent'])>=80 and float(s['empirical_ufboot_percent'])>=95 for s in [sx,sy]);high+=both
            assert (entry['both_meet_existing_support_label']=='True')==both
            q=entry['witness_quartet'].split(';');assert len(set(q))==4 and set(q)<=taxa
            assert {(t in x,t in y) for t in q}=={(False,False),(False,True),(True,False),(True,True)}
        assert int(row['shared_internal_splits'])==len(common)
        assert int(row['unique_internal_splits_a'])==int(row['unique_internal_splits_b'])==523-len(common)
        assert int(row['rf_distance'])==1046-2*len(common)
        assert math.isclose(float(row['normalized_rf']),(1046-2*len(common))/1046)
        assert int(row['incompatible_split_pairs'])==n and int(row['both_supported_incompatible_pairs'])==high
    assert expected==set(indexed)
    assert r['shared_all_three']==len(set.intersection(*(set(s) for s in edges.values())))
    proof=dict(status='passed_full_three_run_pmsf_split_comparison_readback',runs=3,taxa=526,internal_splits_per_run=523,split_presence_rows=len(presence),incompatible_pairs=len(expected),source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),dendropy_version=dendropy.__version__,scope='Independent DendroPy parsing and bipartition compatibility reconstruct all source splits, lengths, reported supports, comparison counts and conflict pairs; empirical support crosschecked against full bootstrap audits. Quartet witnesses and missing cells verified. No inference or SH-aLRT rerun.')
    Path('metadata/pmsf_three_run_topology_readback_20260927.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof,indent=2))


if __name__=='__main__':main()
