#!/usr/bin/env python3
"""Check extracted tree bytes and reconstruct every neighborhood via parent paths."""
import csv,json
from collections import defaultdict,Counter
from pathlib import Path
from Bio import Phylo
from prepare_case_ancestral_neighborhoods import sha,read


def main():
    root=Path('results/ancestral/case-neighborhoods-20260927-v1')
    rp=root/'receipt.json';r=json.loads(rp.read_text())
    for p,h in r['source_hashes'].items():assert sha(p)==h
    for p,h in r['artifacts'].items():assert sha(root/p)==h
    plan=json.loads(Path('metadata/duplication_candidate_tree_review_plan_20260926.json').read_text())
    rows=read(root/'neighborhoods.tsv');desc=read(root/'descendants.tsv');comparison=read(root/'guide_comparison.tsv')
    families={x['family'] for x in rows};assert len(families)==13
    for src in plan['guides']:
        seen=set()
        with Path(src['trees']).open() as f:
            for line in f:
                family,text=line.rstrip().split(': ',1)
                if family in families:
                    assert (root/(src['guide']+'-'+family+'.nwk')).read_text()==text+'\n';seen.add(family)
        assert seen==families
    bytree=defaultdict(list);reported=defaultdict(set)
    for row in rows:bytree[row['tree_file']].append(row)
    for row in desc:
        key=(row['guide'],row['family'],int(row['level']));assert row['gene'] not in reported[key];reported[key].add(row['gene'])
    checked={}
    for filename,group in bytree.items():
        tree=Phylo.read(root/filename,'newick');parents={};nodes={};leaves={};stack=[tree.root]
        while stack:
            node=stack.pop();nodes[node.name]=node
            if not node.clades:leaves[node.name]=node
            for child in node.clades:parents[child]=node;stack.append(child)
        case=group[0]
        def path(gene):
            p=[];n=leaves[gene]
            while True:
                p.append(n)
                if n not in parents:return p
                n=parents[n]
        pa,pb=path(case['gene_a']),set(path(case['gene_b']))
        ancestor=next(n for n in pa if n in pb)
        expected_levels=[]
        for level in range(4):
            expected_levels.append((level,ancestor))
            if ancestor not in parents:break
            ancestor=parents[ancestor]
        assert len(group)==len(expected_levels)
        for level,node in expected_levels:
            row=next(x for x in group if int(x['level'])==level)
            genes={g for g in leaves if node in path(g)}
            key=(row['guide'],row['family'],level)
            assert genes==reported[key] and node.name==row['node']
            assert int(row['proteins'])==len(genes) and int(row['taxa'])==len({g.split('_',1)[0] for g in genes})
            assert int(row['whole_family_proteins'])==len(leaves)
            assert int(row['direct_children'])==len(node.clades)
            assert int(row['is_root'])==(node is tree.root)
            assert int(row['focal_pair_exact'])==(level==0 and genes=={row['gene_a'],row['gene_b']})
            checked[key]=genes
    assert set(checked)==set(reported) and len(checked)==104
    assert sum(map(len,checked.values()))==len(desc)==3722
    counts=Counter()
    for row in comparison:
        a,b=(checked[g,row['family'],int(row['level'])] for g in ['profile','mafft'])
        assert int(row['profile_proteins'])==len(a) and int(row['mafft_proteins'])==len(b)
        assert int(row['shared_proteins'])==len(a&b) and int(row['union_proteins'])==len(a|b)
        assert row['status']==('exact_descendant_set' if a==b else 'different_descendant_set')
        counts[row['status']]+=1
    assert len(comparison)==52
    proof=dict(status='passed_full_case_ancestral_neighborhood_readback',source_receipt_sha256=sha(rp),script_sha256=sha(__file__),trees=len(bytree),neighborhoods=len(checked),descendant_records=len(desc),comparison_counts=dict(counts),scope='Every extracted source tree byte-checked; all MRCA and ancestor memberships rederived from leaf-to-root paths. Same Biopython parser shared. No ancestral inference or independent duplication support.')
    out=Path('metadata/case_ancestral_neighborhoods_completed_20260927.json');assert not out.exists();out.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof,indent=2))


if __name__=='__main__':main()
