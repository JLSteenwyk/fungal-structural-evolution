#!/usr/bin/env python3
"""Visualize matched predictor comparisons with explicit residue selection."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--comparisons',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable figure output')
    receipt=checked_receipt(a.comparisons);rows=read_table(a.comparisons/'comparisons.tsv')
    lookup={(r['sequence_id'],int(r['plddt_cutoff'])):r for r in rows}
    if len(lookup)!=len(rows):raise ValueError('Duplicate control rows')
    chosen=sorted([r for r in rows if r['plddt_cutoff']=='70' and r['status']=='compared'],key=lambda r:r['sequence_id'])
    unfiltered=np.array([float(lookup[(r['sequence_id'],0)]['ca_superposition_rmsd_angstrom']) for r in chosen])
    filtered=np.array([float(r['ca_superposition_rmsd_angstrom']) for r in chosen])
    local=np.array([float(r['pae10_local_mean_absolute_change_angstrom']) for r in chosen])
    fraction=np.array([float(r['retained_fraction']) for r in chosen])
    if not len(chosen) or any(not np.isfinite(v).all() or (v<=0).any() for v in [unfiltered,filtered,local]):raise ValueError('Log plot requires finite positive values; do not silently drop zeros')
    if (fraction<.5).any() or (fraction>1).any():raise ValueError('Invalid retained residue fraction')
    fig,axes=plt.subplots(1,2,figsize=(12,5.7))
    kwargs={'c':fraction,'cmap':'viridis','vmin':.5,'vmax':1.,'s':24,'alpha':.8,'linewidths':.25,'edgecolors':'white'}
    axes[0].scatter(unfiltered,filtered,**kwargs)
    low=min(unfiltered.min(),filtered.min())*.7;high=max(unfiltered.max(),filtered.max())*1.4
    axes[0].plot([low,high],[low,high],color='#8c969f',ls='--',lw=1)
    axes[0].set_xlim(low,high);axes[0].set_ylim(low,high)
    axes[0].set_xlabel('Full-protein CA RMSD (Å)');axes[0].set_ylabel('Joint pLDDT ≥70 CA RMSD (Å)')
    axes[0].set_title(f'Same {len(chosen)} proteins; residue sets change',fontsize=11)
    plot=axes[1].scatter(filtered,local,**kwargs)
    axes[1].set_xlabel('Joint pLDDT ≥70 CA RMSD (Å)');axes[1].set_ylabel('PAE≤10 local distance change (Å)')
    axes[1].set_title('Global versus confident local disagreement',fontsize=11)
    for ax in axes:
        ax.set_xscale('log');ax.set_yscale('log');ax.grid(alpha=.15);ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    fig.suptitle('AlphaFold–ESMFold agreement on exact protein sequences',fontsize=14,y=.97)
    fig.subplots_adjust(left=.08,right=.86,bottom=.2,top=.86,wspace=.36)
    cb=fig.colorbar(plot,cax=fig.add_axes([.9,.26,.018,.5]));cb.set_label('Fraction of residues retained at pLDDT ≥70')
    fig.text(.08,.055,'57 of 266 controls fail the pLDDT ≥70 coverage rule and are absent from both panels.\nPAE applies in both directions in both models. Points are not independent replicates or measurements of experimental accuracy.',fontsize=9,color='#43535f')
    a.output.mkdir(parents=True)
    for ext in ['svg','png','pdf']:fig.savefig(a.output/('predictor_agreement.'+ext),dpi=180)
    plt.close(fig)
    ranked=sorted(chosen,key=lambda r:(-float(r['ca_superposition_rmsd_angstrom']),r['sequence_id']))
    write_table(a.output/'plddt70_geometry_review_order.tsv',[dict(review_rank=i+1,**r) for i,r in enumerate(ranked)])
    result={'status':'complete_predictor_agreement_figure','comparison_receipt_sha256':sha(a.comparisons/'receipt.json'),'script_sha256':sha(Path(__file__)),
        'matplotlib_version':matplotlib.__version__,'plotted_controls':len(chosen),'total_controls':receipt['unique_sequences'],
        'interpretation':'Same pLDDT70-eligible protein cohort in both panels; left contrasts changed residue masks, right contrasts different global/local metrics. Ranking is for artifact/domain inspection, not acceleration, accuracy or significance.',
        'artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
