#!/usr/bin/env python3
"""High-precision replay of a flagged extreme-rate fit; preserve all starts."""
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from Bio import Phylo,SeqIO
from scipy.stats import gamma
from scipy.special import gammainc
from prepare_case_ancestral_neighborhoods import sha


def high_precision(tree,seq,parameters,correction):
    g,l,a=parameters
    rates=4*np.diff(gammainc(a+1,a*gamma.ppf([0,.25,.5,.75,1],a=a,scale=1/a)))
    with mp.workdps(70):
        g,l=mp.mpf(g),mp.mpf(l)
        pi=[l/(g+l),g/(g+l)]
        def prob(pattern,rate):
            def recurse(node):
                if node.is_terminal():
                    return [mp.mpf(int(pattern[node.name] in ['?',str(k)])) for k in [0,1]]
                value=[mp.mpf(1),mp.mpf(1)]
                for child in node.clades:
                    below=recurse(child)
                    changed=-mp.expm1(-(g+l)*mp.mpf(rate)*mp.mpf(child.branch_length))
                    transition=[[1-pi[1]*changed,pi[1]*changed],[pi[0]*changed,1-pi[0]*changed]]
                    for k in [0,1]:
                        value[k]*=sum(transition[k][z]*below[z] for z in [0,1])
                return value
            value=recurse(tree.root)
            return sum(pi[k]*value[k] for k in [0,1])
        total=mp.mpf(0)
        for col in range(len(next(iter(seq.values())))):
            pattern={n:s[col] for n,s in seq.items()}
            excluded={n:('?' if correction=='observed_mask' and x=='?' else '0') for n,x in pattern.items()}
            observed=sum(prob(pattern,r) for r in rates)/4
            exclusion=sum(prob(excluded,r) for r in rates)/4
            total+=mp.log(observed/(1-exclusion))
        return str(total)


def main():
    jid='OG0000152-alignment-famsa-terminal_unknown-observed_mask'
    rp=Path('results/ancestral/conditional-indel-models-20260927-v1')/jid/'receipt.json'
    receipt=json.loads(rp.read_text());job=receipt['job']
    for p,h in job['pins'].items():
        assert sha(p)==h
    tree=Phylo.read(job['tree'],'newick')
    seq={r.id:str(r.seq) for r in SeqIO.parse(job['characters'],'fasta')}
    rows=[]
    for start in receipt['starts']:
        value=high_precision(tree,seq,start['parameters'],job['correction'])
        rows.append(dict(parameters=start['parameters'],high_precision_log_likelihood=value,
            reported=start['log_likelihood'],difference=float(value)-start['log_likelihood']))
    result=dict(job_id=jid,decimal_precision=70,starts=rows,source_receipt_sha256=sha(rp),
        script_sha256=sha(__file__),
        scope='Recursive high-precision binary pruning using the same float gamma category rates. One flagged case, all three fitted starts; no global validation or convergence claim.')
    Path('metadata/conditional_indel_extreme_rate_diagnostic_20260927.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
