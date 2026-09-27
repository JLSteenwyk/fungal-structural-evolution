#!/usr/bin/env python3
"""Check common-tip split projection against independent DendroPy tree pruning."""
import json
from pathlib import Path
import dendropy
from Bio import Phylo
from compare_marker_alignment_topologies import projected_splits,sha


def main():
    source=Path('results/phylogeny/marker-tree-support-complete-v1/receipt.json')
    r=json.loads(source.read_text());cases=0;split_count=0
    for item in r['inputs']:
        path=Path('results/phylogeny/marker-gene-trees-v2')/item['marker']/'tree.treefile'
        assert sha(path)==item['tree_sha256']
        bio=Phylo.read(path,'newick');taxa={x.name for x in bio.get_terminals()}
        for common in [taxa,{x for x in taxa if not x.startswith('O')}]:
            tree=dendropy.Tree.get(path=str(path),schema='newick',rooting='force-unrooted',preserve_underscores=True)
            tree.retain_taxa_with_labels(common);tree.encode_bipartitions()
            other=set()
            for edge in tree.postorder_edge_iter():
                side={x.taxon.label for x in edge.head_node.leaf_iter()}
                complement=common-side
                if min(len(side),len(complement))<2:continue
                other.add(min(tuple(sorted(side)),tuple(sorted(complement)),key=lambda s:(len(s),s)))
            observed=projected_splits(bio,common)
            assert observed==other and len(other)==len(common)-3
            cases+=1;split_count+=len(other)
    out=Path('metadata/marker_topology_projection_readback_20260927.json')
    result=dict(status='passed_full_profile_marker_projection_against_independent_pruning',marker_condition_pairs=cases,internal_splits_checked=split_count,source_receipt_sha256=sha(source),producer_sha256=sha(Path('scripts/compare_marker_alignment_topologies.py')),checker_sha256=sha(Path(__file__)),scope='All 125 completed profile trees, each on all tips and with all outgroups removed; independent DendroPy pruning agrees with Bio.Phylo split restriction. Does not validate unfinished MAFFT results or cross-method comparison.')
    with out.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
