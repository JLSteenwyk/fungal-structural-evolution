#!/usr/bin/env python3
"""Show fixed-nuisance likelihood slices for the 12 largest branch dS estimates."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from audit_genus_codon_trees import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--slices',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new figure path')
    receipt=json.loads((a.slices/'receipt.json').read_text())
    for name in ['case_review.tsv','likelihood_slices.tsv']:
        if sha(a.slices/name)!=receipt['artifacts'][name]:raise ValueError('Changed source table')
    cases=pd.read_csv(a.slices/'case_review.tsv',sep='\t').head(12);points=pd.read_csv(a.slices/'likelihood_slices.tsv',sep='\t')
    fig,axes=plt.subplots(4,3,figsize=(12,10))
    for index,((_,case),ax) in enumerate(zip(cases.iterrows(),axes.flat)):
        data=points[points.case_id==case.case_id].sort_values('branch_multiplier')
        ax.axhline(0,color='#999999',lw=.7);ax.axvline(1,color='#999999',ls='--',lw=.7)
        ax.plot(data.branch_multiplier,data.delta_log_likelihood_from_saved_fit,'o-',color='#176B99',markersize=3)
        ax.set_xscale('log');ax.set_xticks([.1,1,10],['0.1','1','10'])
        ax.set_ylim(min(-.1,data.delta_log_likelihood_from_saved_fit.min()*1.08),max(.02,data.delta_log_likelihood_from_saved_fit.max()*1.08))
        genus,marker,_=case.case_id.split('__');ax.set_title(genus+' | '+marker+'\n'+case.target_node+'; fitted dS = '+format(case.target_dS,'.3g'),fontsize=9)
        if index%3==0:ax.set_ylabel('Δ log likelihood')
        if index>=9:ax.set_xlabel('Branch parameter / fitted value')
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Largest-dS branches: conditional likelihood slices',fontsize=14)
    fig.text(.02,.01,'All other parameters held fixed; panels use independent vertical scales. These are not profile likelihoods or confidence intervals.\nThe full 1,655-case grid is retained in tables; large dS and a peaked conditional slice do not establish absence of synonymous saturation.',fontsize=9)
    fig.tight_layout(rect=(0,.065,1,.95));a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output);print('Saved',a.output)


if __name__=='__main__':main()
