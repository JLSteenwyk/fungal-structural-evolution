#!/usr/bin/env python3
"""Plot audited protein-balanced predictor agreement and threshold attrition."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from audit_joint_path_uncertainty import checked,rows,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['summary','audit','output']:ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=checked(a.summary);audit=json.loads(a.audit.read_text())
    if audit['status']!='passed_hierarchical_control_summary_readback' or audit['summary_receipt_sha256']!=sha(a.summary/'receipt.json'):raise ValueError('Matching summary readback required')
    data=list(rows(a.summary/'protein_summary.tsv'));cohorts=list(rows(a.summary/'cohort_summary.tsv'))
    common=set.intersection(*[{x['sequence_sha256'] for x in data if x['joint_predicted_plddt_cutoff']==c} for c in ['0','70','90']])
    fig,axes=plt.subplots(1,3,figsize=(12,4.7));plt.rcParams['svg.fonttype']='none'
    plotted=[]
    for ax,cutoff,title in zip(axes[:2],['0','90'],['No confidence filter','Both predictors pLDDT ≥90']):
        subset=[x for x in data if x['joint_predicted_plddt_cutoff']==cutoff]
        x=np.array([float(v['af_experiment_rmsd_angstrom']) for v in subset]);y=np.array([float(v['esm_experiment_rmsd_angstrom']) for v in subset])
        if np.any(x<=0) or np.any(y<=0):raise ValueError('Nonpositive RMSD requires different scale')
        lo=min(x.min(),y.min())*.7;hi=max(x.max(),y.max())*1.4
        ax.scatter(x,y,s=28,color='#246b91',alpha=.8,edgecolors='none');ax.plot([lo,hi],[lo,hi],color='#666666',ls='--',lw=1)
        ax.set(xscale='log',yscale='log',xlim=(lo,hi),ylim=(lo,hi),xlabel='AlphaFold–experiment RMSD (Å)',ylabel='ESMFold–experiment RMSD (Å)',title=f'{title}\nn = {len(subset)} proteins');ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.12)
        plotted.append({'panel':title,'points':len(subset),'cutoff':int(cutoff)})
    ax=axes[2]
    for sequence in sorted(common):
        lookup={x['joint_predicted_plddt_cutoff']:float(x['paired_delta_rmsd_angstrom']) for x in data if x['sequence_sha256']==sequence}
        ax.plot([0,1,2],[lookup[c] for c in ['0','70','90']],color='#b9c7ce',alpha=.65,lw=.8)
    medians=[float(next(x for x in cohorts if x['cohort']=='same_proteins_all_thresholds' and x['joint_predicted_plddt_cutoff']==c and x['metric']=='paired_delta_rmsd_angstrom')['median']) for c in ['0','70','90']]
    ax.plot([0,1,2],medians,'o-',color='#b64c2e',lw=2,label='Equal-protein median');ax.axhline(0,color='#555555',ls='--',lw=1)
    ax.set(xticks=[0,1,2],xticklabels=['0','70','90'],xlabel='Joint predicted pLDDT cutoff',ylabel='Paired RMSD difference (Å)\nESMFold minus AlphaFold',title=f'Same {len(common)} proteins at all cutoffs');ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.12)
    fig.suptitle('Matched experimental controls: descriptive protein-level agreement',fontsize=13,y=.98)
    fig.text(.5,.015,'Medians: chains → deposited models → entries → proteins. Thresholds change residues and may change entries.\nSelected reference set; experimental context and training overlap unresolved. Positive paired difference means lower AlphaFold discrepancy.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.10,1,.92));a.output.mkdir(parents=True)
    for ext in ['svg','png']:fig.savefig(a.output/('experimental_predictor_controls.'+ext),dpi=180)
    plt.close(fig)
    receipt={'status':'complete_experimental_predictor_control_figure','summary_receipt_sha256':sha(a.summary/'receipt.json'),'audit_sha256':sha(a.audit),'script_sha256':sha(Path(__file__)),'scatter_panels':plotted,'matched_protein_lines':len(common),'matched_median_delta_rmsd':medians,'interpretation':'Descriptive nested-median summaries; log scatter axes. Matched-protein panel does not hold residue or entry cohort fixed across thresholds. No significance or unbiased accuracy claim.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
