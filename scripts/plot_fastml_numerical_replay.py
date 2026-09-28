#!/usr/bin/env python3
"""Plot fixed-parameter replay discrepancies; no fit-acceptance claim."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha,write_json


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--audit',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    receipt=json.loads((args.audit/'receipt.json').read_text());assert receipt['status']=='complete_native_rate_replay_provenance_and_scalar_readback'
    for name,h in receipt['artifacts'].items():assert sha(args.audit/name)==h
    table=pd.read_csv(args.audit/'replay_comparison.tsv',sep='\t');assert len(table)==306
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    for ax,kind,title in zip(axes,['likelihood','probability'],['Log-likelihood discrepancy','Largest marginal probability discrepancy']):
        for variant,color in [('precision-only','#d06b28'),('precision-cache-refresh','#2375a9')]:
            f=table[table.variant.eq(variant)];assert len(f)==153
            x=f['scipy_rate_'+('likelihood_difference' if kind=='likelihood' else 'maximum_probability_difference')].abs().to_numpy()
            y=f['native_rate_'+('likelihood_difference' if kind=='likelihood' else 'maximum_probability_difference')].abs().to_numpy()
            ax.scatter(np.maximum(x,1e-16),np.maximum(y,1e-16),s=20,alpha=.65,label=variant,color=color)
        low,high=(1e-16,1) if kind=='likelihood' else (1e-16,1e-3)
        ax.plot([low,high],[low,high],color='.55',linestyle=':',linewidth=1)
        ax.set(xscale='log',yscale='log',xlim=(low,high),ylim=(low,high),title=title,xlabel='Replay with SciPy rates',ylabel='Replay with native rates')
        ax.spines[['top','right']].set_visible(False)
    axes[0].legend(frameon=False,fontsize=8,loc='upper left')
    fig.suptitle('FastML numerical replay: 153 nonempty fits per variant',fontsize=12)
    args.output.mkdir(parents=True,exist_ok=False)
    for suffix in ['png','pdf']:fig.savefig(args.output/('numerical_replay.'+suffix),dpi=200)
    plt.close(fig)
    write_json(args.output/'receipt.json',dict(status='numerical_replay_figure_generated',source_audit_sha256=sha(args.audit/'receipt.json'),script_sha256=sha(__file__),points_per_panel=len(table),plot_floor=1e-16,artifacts={p.name:sha(p) for p in args.output.iterdir()},caption='Absolute reported-minus-replayed discrepancies at fixed fitted parameters. Each point is one fit; dots can overlap. Values below 1e-16 are displayed at 1e-16 only. Dotted line is equality. Native rates isolate gamma-discretization effects; three uncorrected stale-cache likelihood discrepancies persist. Agreement is numerical reproduction, not optimization or biological qualification.'))


if __name__=='__main__':main()
