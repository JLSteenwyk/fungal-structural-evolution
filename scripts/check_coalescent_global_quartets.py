#!/usr/bin/env python3
"""Enumerate all 105 six-tip species topologies and actual quartet score checks."""
import argparse
import itertools
import json
from pathlib import Path
import time

from coalescent_global_quartet_audit import matching_quartets,species_tripartitions
from coalescent_quartet_audit_v2 import packed,read_tree,resolved_quartet_denominator
from run_ortholog_pair_guide_comparison import sha


def all_six_tip_trees():
    graphs = [{0:{100},1:{100},2:{100},100:{0,1,2}}]
    for leaf in range(3,6):
        expanded = []
        for graph in graphs:
            for a,b in [(a,b) for a,v in graph.items() for b in v if a < b]:
                g = {k:set(v) for k,v in graph.items()}
                internal = 100+leaf-2
                g[a].remove(b);g[b].remove(a)
                g[a].add(internal);g[b].add(internal)
                g[internal] = {a,b,leaf};g[leaf] = {internal}
                expanded.append(g)
        graphs = expanded
    def subtree(g,node,previous=None):
        if node < 100: return chr(ord('A')+node)
        return '('+','.join(sorted(subtree(g,n,node) for n in g[node] if n != previous))+')'
    result = [subtree(g,100)+';' for g in graphs]
    assert len(result) == len(set(result)) == 105
    return result


def bipartitions(tree):
    descendants = {}
    for node in tree.postorder_node_iter():
        descendants[node] = ({node.taxon.label} if node.is_leaf() else
                             set.union(*(descendants[c] for c in node.child_node_iter())))
    return list(descendants.values()),descendants[tree.seed_node]


def projected_quartet(sides,quad):
    q = set(quad)
    pairs = set()
    for side in sides:
        part = side&q
        if len(part) == 2: pairs.add(tuple(sorted(part if quad[0] in part else q-part)))
    assert len(pairs) <= 1
    return next(iter(pairs)) if pairs else None


def brute_common(species,gene):
    s,_ = bipartitions(species)
    g,tips = bipartitions(gene)
    score = 0
    for quad in itertools.combinations(sorted(tips),4):
        expected = projected_quartet(g,quad)
        if expected is not None and expected == projected_quartet(s,quad): score += 1
    return score


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    genes = [
        '((A,B),(C,D),(E,F));','((A,C),B,D,E,F);','(A,B,C,D,E,F);',
        '(A,B,(C,D,E,F));','((A,B,C),(D,E,F));','(A,((B,C),D,(E,F)));',
        '(A,(((E,F),(C,D)),B));','((A,B),(C,D));','((A,C),B,D,E);',
        '((A,B),C,D,E);','(A,B,C,D);','((A,B),C,D,(E,F));']
    cases=0
    for newick in all_six_tip_trees():
        species=read_tree(data=newick);tips=list('ABCDEF')
        tripartitions=species_tripartitions(species,tips)
        for g in genes:
            gene=read_tree(data=g)
            actual=int(matching_quartets(tripartitions,*packed(gene,tips)))
            expected=brute_common(species,gene)
            assert actual == expected, (newick,g,actual,expected)
            assert actual <= resolved_quartet_denominator(gene)
            cases += 1
    large_species='((((A,B),(C,D)),((E,F),(G,H))),((I,J),(K,L)));'
    large_genes=[large_species,'(A,B,C,D,E,F,G,H,I,J,K,L);',
                 '((A,C,E),(B,D),((F,G,H),(I,J)),K,L);',
                 '(((A,D),(B,C)),(E,(I,J)),(K,L));']
    species=read_tree(data=large_species);tips=list('ABCDEFGHIJKL')
    tripartitions=species_tripartitions(species,tips)
    for g in large_genes:
        gene=read_tree(data=g)
        actual=int(matching_quartets(tripartitions,*packed(gene,tips)))
        assert actual == brute_common(species,gene)
        cases += 1
    result=dict(status='passed_independent_global_matching_quartet_exhaustive_topology_contracts',
                complete_six_tip_species_topologies=105,gene_cases_per_species=12,
                tree_pair_cases=cases,large_missing_multifurcating_gene_cases=4,
                elapsed_seconds=time.monotonic()-started,scientific_eligibility=False,
                source_hashes={p:sha(p) for p in ['scripts/coalescent_global_quartet_audit.py',
                   'scripts/coalescent_quartet_audit_v2.py',str(Path(__file__))]},
                scope='Actual four-tip subset enumeration and original split projection agrees with node-component DP on every 105 six-tip species topology x 12 resolved, partial and multifurcating gene cases, plus four twelve-tip cases. Software contract only; production scoring and biological qualification remain separate.')
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,float,bool))}),flush=True)


if __name__ == '__main__': main()
