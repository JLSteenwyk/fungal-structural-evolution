#!/usr/bin/env python3
"""Check all node marginals against 70-digit enumeration of complete states."""
import json
from io import StringIO
from itertools import product
from pathlib import Path
import mpmath as mp
import numpy as np
from Bio import Phylo
from scipy.special import gammainc
from scipy.stats import gamma
from stable_indel_posteriors import infer
from prepare_case_ancestral_neighborhoods import sha


def main():
    tree=Phylo.read(StringIO('(a:0.2,(b:0.3,c:0.4)I:0.1)R;'),'newick')
    seq={'a':'01?1','b':'1?0?','c':'001?'}
    nodes=list(tree.find_clades());checks=[]
    with mp.workdps(70):
        for g,l,a in [[.4,.9,.5],[.00001,160000.,100.],[.2,3.,.005]]:
            rates=4*np.diff(gammainc(a+1,a*gamma.ppf([0,.25,.5,.75,1],a=a,scale=1/a)))
            posterior,logs=infer(tree,seq,g,l,rates)
            gain,loss=mp.mpf(g),mp.mpf(l);pi=[loss/(gain+loss),gain/(gain+loss)]
            max_prob,max_ll=0.,0.
            for col in range(4):
                total=mp.mpf(0);ones={n.name:mp.mpf(0) for n in nodes}
                for rate in rates:
                    transitions={}
                    for node in nodes:
                        if node is tree.root:
                            continue
                        decay=mp.exp(-(gain+loss)*mp.mpf(rate)*mp.mpf(node.branch_length))
                        transitions[node]=[[pi[z]+(int(k==z)-pi[z])*decay for z in [0,1]] for k in [0,1]]
                    for states in product([0,1],repeat=len(nodes)):
                        assignment=dict(zip(nodes,states))
                        if any(seq[n.name][col]!='?' and int(seq[n.name][col])!=assignment[n] for n in tree.get_terminals()):
                            continue
                        weight=pi[assignment[tree.root]]/4
                        for node in nodes:
                            for child in node.clades:
                                weight*=transitions[child][assignment[node]][assignment[child]]
                        total+=weight
                        for node in nodes:
                            ones[node.name]+=weight*assignment[node]
                max_ll=max(max_ll,abs(logs[col]-float(mp.log(total))))
                for name,value in ones.items():
                    max_prob=max(max_prob,abs(posterior[name][col]-float(value/total)))
            assert max_ll<1e-11 and max_prob<1e-11,(g,l,a,max_ll,max_prob)
            checks.append(dict(gain=g,loss=l,alpha=a,sites=4,nodes=5,
                maximum_log_likelihood_error=max_ll,maximum_probability_error=max_prob))
    result=dict(status='all_synthetic_node_probabilities_pass_70_digit_enumeration',
        checks=checks,probabilities_checked=60,decimal_precision=70,
        pins={p:sha(p) for p in [__file__,'scripts/stable_indel_posteriors.py']},
        scope='Analytic transition inside/outside marginals versus exhaustive high-precision state enumeration with unknown tips, rare events and gamma extremes. Not production posterior validation, ancestral sequence assembly or model adequacy.')
    Path('metadata/stable_indel_posterior_validation_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],max(r['maximum_probability_error'] for r in checks))


if __name__=='__main__':
    main()
