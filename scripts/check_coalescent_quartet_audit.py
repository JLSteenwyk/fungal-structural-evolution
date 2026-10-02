#!/usr/bin/env python3
"""Exhaustive independent quartet-DP checks plus actual native numeric contracts."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from coalescent_quartet_audit import (annotated_branches, colored_quartets, close,
                                     numerical_check, packed, read_tree,
                                     resolved_quartet_denominator)
from run_ortholog_pair_guide_comparison import sha


def brute(tree, tips, colors):
    """Enumerate actual four-tip sets, projecting original gene splits."""
    index = {t:i for i,t in enumerate(tips)}
    descendants = {}
    for n in tree.postorder_node_iter():
        descendants[n] = ({n.taxon.label} if n.is_leaf() else
                          set.union(*(descendants[c] for c in n.child_node_iter())))
    found = sorted(descendants[tree.seed_node])
    counts, denominator = [0,0,0], 0
    for quartet in itertools.combinations(found,4):
        taxon_set = set(quartet)
        projected = set()
        for side in descendants.values():
            part = side & taxon_set
            if len(part) == 2:
                pair = frozenset(part if quartet[0] in part else taxon_set-part)
                projected.add(pair)
        assert len(projected) <= 1
        if not projected: continue
        denominator += 1
        if len({int(colors[index[t]]) for t in quartet}) != 4: continue
        pair = set(next(iter(projected)))
        pair_colors = {int(colors[index[t]]) for t in pair}
        if 0 not in pair_colors: pair_colors = set(range(4))-pair_colors
        other = next(c for c in pair_colors if c != 0)
        counts[other-1] += 1
    available = np.prod([sum(colors[index[t]] == c for t in found) for c in range(4)])
    return counts,int(available),denominator


def check_tree(genes, native_tree, expected_tips):
    summaries = []
    tips = sorted(expected_tips)
    packs = [packed(g,tips) for g in genes]
    for branch in annotated_branches(native_tree,tips):
        counts = np.zeros(3)
        for gene, arrays in zip(genes,packs):
            actual,available = colored_quartets(*arrays,branch['colors'])
            expected,expected_available,_ = brute(gene,tips,branch['colors'])
            assert actual.tolist() == expected and available == expected_available
            if available: counts += actual/available
        summary = numerical_check(branch['native'],branch['length'],counts)
        summaries.append(dict(split=branch['split'],groups=branch['groups'],native=branch['native'],**summary))
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--java',default='/usr/bin/java')
    parser.add_argument('--jar',required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    tips = list('ABCDEF')
    trees = [
        '((A,B),(C,D),(E,F));', '((A,C),B,D,E,F);',
        '(A,B,C,D,E,F);', '(A,B,(C,D,E,F));',
        '((A,B,C),(D,E,F));', '(A,((B,C),D,(E,F)));',
        '(A,(((E,F),(C,D)),B));', '((A,B),(C,D));',
        '((A,C),B,D,E);', '((A,B),C,D,E);',
        '(A,B,C,D);', '((A,B),C,D,(E,F));',
    ]
    exhaustive = 0
    started = time.monotonic()
    for newick in trees:
        tree = read_tree(data=newick)
        arrays = packed(tree,tips)
        expected_denominator = None
        for assignment in itertools.product(range(4),repeat=6):
            if len(set(assignment)) != 4: continue
            colors = np.array(assignment,dtype=np.int64)
            actual,available = colored_quartets(*arrays,colors)
            expected,expected_available,denominator = brute(tree,tips,colors)
            assert actual.tolist() == expected, (newick,assignment,actual,expected)
            assert available == expected_available
            assert resolved_quartet_denominator(tree) == denominator
            expected_denominator = denominator
            exhaustive += 1
        print('exhaustive_quartet_case_passed',newick,expected_denominator,flush=True)
    synthetic_cases = {
        'seven_concordant': [trees[0]]*7,
        'fractional_polytomy': [trees[0],trees[1]],
        'missing_and_star': [trees[0],trees[2],trees[7]],
        'zero_evidence': [trees[2]],
        'discordant': [trees[0],'((A,C),(B,D),(E,F));','((A,D),(B,C),(E,F));',trees[1]],
    }
    reference = args.output/'reference.tree'
    reference.write_text(trees[0]+'\n')
    native_results = []
    mutations_rejected = 0
    for name,newicks in synthetic_cases.items():
        folder = args.output/name
        folder.mkdir()
        gp,sp,lp = folder/'genes.tree',folder/'species.tree',folder/'native.log'
        gp.write_text('\n'.join(newicks)+'\n')
        command = [args.java,'-XX:ActiveProcessorCount=1','-Xmx1G','-jar',args.jar,
                   '-i',str(gp.resolve()),'-q',str(reference.resolve()),'-o',str(sp.resolve()),'-t','2','-s','20261002']
        with lp.open('x') as handle:
            subprocess.run(command,stdout=handle,stderr=subprocess.STDOUT,check=True,timeout=120)
        genes = [read_tree(data=n) for n in newicks]
        summaries = check_tree(genes,read_tree(path=sp),tips)
        denominator = sum(resolved_quartet_denominator(g) for g in genes)
        assert 'ASTRAL version 5.7.8' in lp.read_text()
        for summary in summaries:
            native = summary['native']
            counts = summary['independent_counts_canonical_groups']
            for key in ['f1','f2','f3','pp1','pp2','pp3','EN']:
                changed = dict(native)
                changed[key] += .125
                try: numerical_check(changed,summary['independent_map_length'],counts)
                except AssertionError: mutations_rejected += 1
                else: raise AssertionError('Accepted changed native value: '+key)
            try: numerical_check(native,summary['independent_map_length']+.125,counts)
            except AssertionError: mutations_rejected += 1
            else: raise AssertionError('Accepted changed native length')
        native_results.append(dict(case=name,command=command,denominator=denominator,branches=summaries,
                                   artifacts={str(p):sha(p) for p in [gp,sp,lp]}))
        print('native_quartet_numeric_contract_passed',name,flush=True)
    result = dict(status='passed_exhaustive_quartet_dp_and_native_coalescent_numeric_contracts',
                  exhaustive_color_tree_cases=exhaustive,native_cases=len(native_results),
                  altered_native_values_rejected=mutations_rejected,elapsed_seconds=time.monotonic()-started,
                  native_results=native_results,source_hashes={str(p):sha(p) for p in
                    [Path(__file__),Path('scripts/coalescent_quartet_audit.py'),reference,Path(args.jar),Path(args.java)]},
                  scientific_eligibility=False,
                  scope='Software contracts only, including actual native fractional and zero-evidence outputs. Exhaustive color assignments independently enumerate and project all four-tip subsets. No biological pilot or proof of production inference/model adequacy.')
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,float,bool))}),flush=True)


if __name__ == '__main__': main()
