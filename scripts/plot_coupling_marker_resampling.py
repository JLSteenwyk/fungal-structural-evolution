#!/usr/bin/env python3
"""Plot observed-mean contrasts and marker influence for adjusted FreeRate fits."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from audit_genus_codon_trees import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--title',default='Marker resampling: composition-adjusted FreeRate models')
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a new figure path')
    receipt=json.loads((args.results/'receipt.json').read_text());path=args.results/'contrast_summary.tsv'
    if sha(path)!=receipt['artifacts'][path.name]:raise ValueError('Changed summary')
    d=pd.read_csv(path,sep='\t');d=d[(d.rate_model=='freerate')&(d.controls=='composition_adjusted')]
    rows=[(a,s) for a in ['3di_af','3di_af_empirical','3di_llm'] for s in ['tien2013','miller1987']]
    labels={'3di_af':'AF','3di_af_empirical':'AF + empirical F','3di_llm':'LLM','tien2013':'Tien','miller1987':'Miller'}
    contrasts=[('aa_log1p_rate_at_observed_mean_RSA','Sequence-rate association\nat observed mean RSA'),('rsa_at_observed_mean_aa_log1p_rate','RSA association\nat observed mean log(1 + AA rate)'),('interaction','Sequence-rate × RSA\ninteraction')]
    fig,axes=plt.subplots(1,3,figsize=(12,5),sharey=True)
    for ax,(contrast,title) in zip(axes,contrasts):
        ax.axvline(0,color='#777777',linestyle='--',linewidth=.8)
        for y,(alphabet,scale) in enumerate(rows):
            r=d[(d.alphabet_model==alphabet)&(d.rsa_scale==scale)&(d.contrast==contrast)].iloc[0]
            color='#176B99' if scale=='tien2013' else '#C46D28'
            ax.plot([r.bootstrap_percentile_lower,r.bootstrap_percentile_upper],[y,y],color=color,lw=1.5)
            ax.plot([r.loo_minimum,r.loo_maximum],[y,y],color='#222222',lw=4,solid_capstyle='butt')
            ax.plot(r.point_estimate,y,'o',color=color,markersize=6)
        ax.set_title(title,fontsize=10);ax.set_xlabel('Conditional coefficient');ax.grid(axis='y',alpha=.15)
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_yticks(range(len(rows)),[labels[a]+' | '+labels[s] for a,s in rows],fontsize=9);axes[0].invert_yaxis()
    fig.suptitle(args.title,fontsize=13)
    handles=[Line2D([0],[0],color='#176B99',lw=1.5,label='Unadjusted 95% marker-bootstrap interval'),Line2D([0],[0],color='#222222',lw=4,label='Leave-one-marker-out range')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.62,.90),ncol=1,fontsize=9)
    fig.text(.02,.02,'2,000 paired marker resamples; 72 omissions per model. Reference means held fixed at the observed site-weighted distribution.\nAll 24 specifications are retained in tables; these panels show the six adjusted FreeRate fits.\nConditional on estimated rates and fixed trees; not causal effects or physical displacement.',fontsize=9)
    fig.tight_layout(rect=(0,.14,1,.77));args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(args.output)
    print('Saved',args.output)


if __name__=='__main__':main()
