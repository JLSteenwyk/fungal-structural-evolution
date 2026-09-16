#!/usr/bin/env python3
"""Plot audited combined controls with length tiers and confidence-dependent coverage."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from audit_joint_path_uncertainty import checked,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['summary','audit','output']:ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=checked(a.summary);audit=json.loads(a.audit.read_text())
    if audit['status']!='passed_combined_control_summary_readback' or audit['summary_receipt_sha256']!=sha(a.summary/'receipt.json'):raise ValueError('Matching combined summary readback required')
    data=pd.read_csv(a.summary/'protein_summary.tsv',sep='\t');dis=pd.read_csv(a.summary/'protein_disposition.tsv',sep='\t');cohorts=pd.read_csv(a.summary/'cohort_summary.tsv',sep='\t')
    tiers=[('predictor-controls-summary-v1','≤512 aa','#0072B2'),('predictor-controls-summary-513-768-v1','513–768 aa','#D55E00'),('predictor-controls-summary-769-1152-v1','769–1,152 aa','#009E73'),('predictor-controls-summary-1153-1536-v1','1,153–1,536 aa','#CC79A7')]
    if set(dis.source_tier)!={t[0] for t in tiers}:raise ValueError('Unexpected length tier')
    cut='joint_predicted_plddt_cutoff';metric='paired_delta_rmsd_angstrom';common=set.intersection(*[set(data.loc[data[cut]==c,'sequence_sha256']) for c in [0,70,90]])
    plt.rcParams.update({'svg.fonttype':'none','font.size':9})
    fig,axs=plt.subplots(2,3,figsize=(14,9));coords=[];coverage=[]
    full_values=data[['af_experiment_rmsd_angstrom','esm_experiment_rmsd_angstrom']].to_numpy();assert np.all(full_values>0)
    lo=float(full_values.min()*.7);hi=float(full_values.max()*1.3)
    for index,cutoff in enumerate([0,70,90]):
        ax=axs[0,index];subset=data[data[cut]==cutoff]
        for tier,label,color in tiers:
            part=subset[subset.source_tier==tier]
            ax.scatter(part.af_experiment_rmsd_angstrom,part.esm_experiment_rmsd_angstrom,s=25,color=color,alpha=.75,edgecolors='none')
            for row in part.to_dict('records'):
                coords.append({'sequence_sha256':row['sequence_sha256'],'source_tier':tier,'cutoff':cutoff,'af_experiment_rmsd_angstrom':row['af_experiment_rmsd_angstrom'],'esm_experiment_rmsd_angstrom':row['esm_experiment_rmsd_angstrom'],'paired_delta_rmsd_angstrom':row[metric],'common_all_cutoffs':row['sequence_sha256'] in common})
        ax.plot([lo,hi],[lo,hi],ls='--',color='.5',lw=1);ax.set(xscale='log',yscale='log',xlim=(lo,hi),ylim=(lo,hi),xlabel='AlphaFold–experiment RMSD (Å)',ylabel='ESMFold–experiment RMSD (Å)',title=f'{chr(65+index)}  Joint pLDDT ≥{cutoff} · {len(subset)} proteins');ax.set_aspect('equal');ax.grid(alpha=.12)
    ax=axs[1,0];locations=np.arange(4)
    for j,cutoff in enumerate([0,70,90]):
        fractions=[]
        for tier,_,_ in tiers:
            group=dis[(dis.source_tier==tier)&(dis[cut]==cutoff)];n=int((group.status=='included').sum());fractions.append(100*n/len(group));coverage.append({'source_tier':tier,'cutoff':cutoff,'included':n,'universe':len(group)})
        bars=ax.barh(locations+(j-1)*.24,fractions,height=.22,color=['#b9cbd9','#668ca8','#254d6b'][j],label=str(cutoff))
        for b,item in zip(bars,coverage[-4:]):ax.text(b.get_width()+1,b.get_y()+b.get_height()/2,f"{item['included']}/{item['universe']}",va='center',fontsize=8)
    ax.set(yticks=locations,yticklabels=[t[1] for t in tiers],xlim=(0,125),xticks=[0,25,50,75,100],xlabel='Eligible proteins (% of each tier)\nTop/middle/bottom bars: pLDDT 0/70/90',title='D  Coverage after residue and confidence filters');ax.invert_yaxis();ax.spines[['top','right']].set_visible(False)
    for column,only_common,title in [(1,True,f'E  Same {len(common)} proteins at all cutoffs'),(2,False,'F  All eligible proteins at each cutoff')]:
        ax=axs[1,column]
        subset=data[data.sequence_sha256.isin(common)] if only_common else data
        if only_common:
            for sid,group in subset.groupby('sequence_sha256'):
                ax.plot(range(3),group.set_index(cut).loc[[0,70,90],metric],color='.75',lw=.7,alpha=.7,zorder=1)
        for tier,_,color in tiers:
            group=subset[subset.source_tier==tier]
            x=[([0,70,90].index(int(c))+(int(hashlib.sha256(s.encode()).hexdigest()[:8],16)/2**32-.5)*.26) for c,s in zip(group[cut],group.sequence_sha256)]
            ax.scatter(x,group[metric],s=19,color=color,alpha=.7,zorder=2)
        cohort='same_proteins_all_thresholds' if only_common else 'all_available'
        med=cohorts[(cohorts.source_tier=='combined')&(cohorts.cohort==cohort)&(cohorts.metric==metric)].set_index(cut).loc[[0,70,90],'median']
        ax.plot(range(3),med,'k_-',ms=16,lw=1.4,label='Equal-protein median',zorder=3)
        ax.axhline(0,color='.45',ls='--',lw=.8);ax.set(yscale='symlog',xticks=[0,1,2],xticklabels=['0','70','90'],xlabel='Joint predicted pLDDT cutoff',ylabel='Paired RMSD difference (Å)\nESMFold–experiment minus AlphaFold–experiment',title=title);ax.set_yscale('symlog',linthresh=1);ax.set_ylim(float(data[metric].min()*1.3),float(data[metric].max()*1.3));ax.grid(alpha=.12);ax.legend(frameon=False,fontsize=8)
    handles=[Line2D([0],[0],marker='o',lw=0,color=color,label=label) for _,label,color in tiers]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.951),ncol=4,frameon=False,title='Full-protein length tier')
    fig.suptitle('79 experimental-reference controls: prediction agreement and eligibility',fontsize=15,y=.99)
    fig.text(.5,.012,'Each point is a canonical protein; medians aggregate chains → deposited models → entries → protein. Scatter axes are logarithmic.\nDifference axes use a symmetric log scale, linear within ±1 Å. Confidence filters change residues and may change entries.\nSelected reference set; experimental context and training overlap unresolved. Positive paired difference means lower AlphaFold discrepancy.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.09,1,.88),h_pad=2.6,w_pad=2.2);a.output.mkdir(parents=True)
    pd.DataFrame(coords).to_csv(a.output/'plotted_proteins.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(a.output/'plotted_coverage.tsv',sep='\t',index=False)
    for ext in ['svg','png']:fig.savefig(a.output/('combined_experimental_controls.'+ext),dpi=180)
    plt.close(fig)
    result={'status':'complete_combined_experimental_control_figure','proteins':r['control_sequences'],'plotted_protein_threshold_rows':len(coords),'common_proteins':len(common),'source_summary_sha256':sha(a.summary/'receipt.json'),'source_audit_sha256':sha(a.audit),'script_sha256':sha(Path(__file__)),'interpretation':'Descriptive paired protein summaries; length tiers and exclusions explicit. No confidence interval, significance, causal filter/length effect, unbiased accuracy or training-independence claim.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
