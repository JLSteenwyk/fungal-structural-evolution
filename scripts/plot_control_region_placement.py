#!/usr/bin/env python3
"""Plot a selected experimental-control placement discrepancy with its PAE context."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from audit_joint_path_uncertainty import checked,rows,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['summary','summary-audit','pae-summary','pae-audit','controls','output']:
        p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--alphafold-model',required=True);p.add_argument('--cutoff',type=int,default=90)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    summary_audit=json.loads(a.summary_audit.read_text());pae_audit=json.loads(a.pae_audit.read_text())
    checked(a.controls);pr=checked(a.pae_summary)
    if (summary_audit['status']!='passed_full_segment_grid_count_eligibility_readback'
        or summary_audit['summary_sha256']!=sha(a.summary)
        or pae_audit['status']!='passed_all_control_pae_summary_rows'
        or pae_audit['source_receipt_sha256']!=sha(a.pae_summary/'receipt.json')
        or pr['source_receipts']['controls']!=sha(a.controls/'receipt.json')):
        raise ValueError('Matching reviewed sources required')
    data=sorted([r for r in rows(a.summary) if r['alphafold_model_id']==a.alphafold_model and int(r['joint_predicted_plddt_cutoff'])==a.cutoff],key=lambda r:int(r['alignment_start']))
    if not data or len({r['sequence_sha256'] for r in data})!=1:raise ValueError('Missing or ambiguous selected protein')
    model=next(m for m in json.loads((a.controls/'model_provenance.json').read_text()) if m['sequence_sha256']==data[0]['sequence_sha256'])
    path=Path(model['local_pae_npz_path'])
    if sha(path)!=model['local_pae_npz_sha256']:raise ValueError('Changed PAE source')
    with np.load(path,allow_pickle=False) as d:pae=d['pae'].copy()
    n=model['length'];plt.rcParams['svg.fonttype']='none'
    fig,axes=plt.subplots(1,2,figsize=(11,5.3),gridspec_kw={'width_ratios':[1,1.15]})
    ax=axes[0];im=ax.imshow(pae,origin='upper',extent=(.5,n+.5,n+.5,.5),vmin=0,vmax=float(np.ceil(pae.max())),cmap='viridis',interpolation='nearest',rasterized=True)
    for i,r in enumerate(data,1):
        start,end=int(r['alignment_start']),int(r['alignment_end']);width=end-start+1
        ax.add_patch(Rectangle((start-.5,start-.5),width,width,fill=False,color='white',lw=.8))
        ax.text((start+end)/2,(start+end)/2,str(i),color='white',ha='center',va='center',fontsize=8)
    ax.set(xlabel='Column residue',ylabel='Row residue',title='ESMFold predicted aligned error\nFull sequence; no confidence filter')
    fig.colorbar(im,ax=ax,fraction=.046,pad=.04,label='Predicted aligned error (Å)')
    ax=axes[1];x=np.arange(len(data));plotted=[]
    for predictor,color,label in [('af','#246b91','AlphaFold'),('esm','#bc4b2b','ESMFold')]:
        for suffix,marker,style,description in [('segment_rmsd','o','-','Separate region fit'),('after_whole_mask_fit_rmsd','s','--','Whole-mask fit')]:
            values=[float(r[predictor+'_experiment_'+suffix]) for r in data]
            if any(v<=0 or not np.isfinite(v) for v in values):raise ValueError('Invalid log-axis geometry')
            ax.plot(x,values,marker=marker,ls=style,color=color,label=label+': '+description,markersize=5,lw=1.2,markerfacecolor=color if marker=='o' else 'white')
            plotted.append(dict(predictor=label,fit=description,values=values))
    ax.set(yscale='log',ylabel='Region–experiment CA RMSD (Å)',xticks=x,
        xticklabels=[str(i)+'  '+r['pfam_name']+'\n'+r['alignment_start']+'–'+r['alignment_end'] for i,r in enumerate(data,1)],
        title=f'Annotated-region agreement\nJoint focal pLDDT ≥{a.cutoff}')
    ax.tick_params(axis='x',labelsize=8);ax.grid(axis='y',alpha=.2)
    ax.legend(loc='upper left',bbox_to_anchor=(0,-.18),fontsize=8,frameon=False,ncol=2)
    fig.suptitle(a.alphafold_model+': local agreement and uncertain relative placement',fontsize=13,y=.98)
    fig.text(.5,.015,'Exploratory selected case. Regions include Pfam repeats/families; separate fits mechanically favor lower RMSD.\nGeometry uses observed matched residues and nested chain/model/entry medians. PAE is predicted uncertainty; neither panel proves biological motion.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.09,1,.92));a.output.mkdir(parents=True)
    for ext in ['png','svg']:fig.savefig(a.output/('control_region_placement.'+ext),dpi=180)
    plt.close(fig)
    result=dict(status='complete_control_region_placement_figure',alphafold_model=a.alphafold_model,sequence_sha256=data[0]['sequence_sha256'],regions=len(data),cutoff=a.cutoff,plotted_values=plotted,
        source_hashes={str(p):sha(p) for p in [a.summary,a.summary_audit,a.pae_summary/'receipt.json',a.pae_audit,a.controls/'receipt.json',path]},script_sha256=sha(Path(__file__)),
        artifacts={f.name:sha(f) for f in a.output.iterdir()})
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
