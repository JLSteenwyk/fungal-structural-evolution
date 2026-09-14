#!/usr/bin/env python3
"""Plot all planned conditional coupling coefficients and marker-cluster intervals."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from audit_genus_codon_trees import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fits',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--title',default=None)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a new figure path')
    receipt=json.loads((args.fits/'receipt.json').read_text())
    p=args.fits/'focal_coefficients.tsv'
    if sha(p)!=receipt['artifacts'][p.name]:raise ValueError('Changed coefficient table')
    d=pd.read_csv(p,sep='\t');rows=[]
    alphabet_labels={'3di_af':'AF','3di_af_empirical':'AF + empirical F','3di_llm':'LLM'}
    for alphabet in alphabet_labels:
        for rate in ['gamma','freerate']:
            for controls in ['coverage_confidence','composition_adjusted']:
                label=alphabet_labels[alphabet]+' | '+('Gamma' if rate=='gamma' else 'FreeRate')+' | '+('basic' if controls=='coverage_confidence' else '+ composition')
                rows.append((alphabet,rate,controls,label))
    fig,axes=plt.subplots(1,3,figsize=(13,7),sharey=True)
    terms=[('aa_log1p_rate','Sequence-rate association\nat RSA = 0.25'),('rsa_centered','RSA association\nat log(1 + AA rate) = 0'),('aa_by_rsa','Sequence-rate × RSA\ninteraction')]
    for ax,(term,title) in zip(axes,terms):
        ax.axvline(0,color='#777777',lw=.8,ls='--',label='Zero')
        for scale,offset,color,marker in [('tien2013',-.13,'#176B99','o'),('miller1987',.13,'#C46D28','s')]:
            subset=d[(d.term==term)&(d.rsa_scale==scale)]
            for y,(alphabet,rate,controls,label) in enumerate(rows):
                r=subset[(subset.alphabet_model==alphabet)&(subset.rate_model==rate)&(subset.controls==controls)].iloc[0]
                ax.errorbar(r.coefficient,y+offset,xerr=[[r.coefficient-r.ci95_lower],[r.ci95_upper-r.coefficient]],fmt=marker,color=color,markersize=4,capsize=2,lw=1,label=scale if y==0 else None)
        ax.set_title(title,fontsize=11)
        ax.set_xlabel('Coefficient (log-transformed relative rates)')
        ax.grid(axis='y',alpha=.15)
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_yticks(np.arange(len(rows)),[r[3] for r in rows],fontsize=9)
    axes[0].invert_yaxis();fig.legend(*axes[2].get_legend_handles_labels(),fontsize=9,loc='upper center',bbox_to_anchor=(.64,.94),ncol=3)
    fig.suptitle(args.title or f"Conditional associations in {receipt['sites']:,} paired sites from {receipt['markers']} markers",fontsize=13)
    fig.text(.02,.018,'Marker fixed intercepts; 95% marker-cluster intervals conditional on fixed trees and estimated rates.\nBasic controls: coverage and focal CA confidence. Composition model also includes AA proportions and entropy.\nAF and LLM identify alphabet substitution models; every structure here was predicted by ESMFold.',fontsize=9)
    fig.tight_layout(rect=(0,.105,1,.94))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output)
    print('Saved',args.output)


if __name__=='__main__':main()
