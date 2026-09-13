#!/usr/bin/env python3
"""Plot changes in paired coverage and distinguish presence from marker depth."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['previous','expanded','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new figure output')
    old_receipt=checked_receipt(a.previous);new_receipt=checked_receipt(a.expanded)
    if old_receipt['manifest_sha256']!=new_receipt['manifest_sha256']:raise ValueError('Different taxon manifests')
    old={(r['study_role'],r['lineage_group']):r for r in read_table(a.previous/'lineage_coverage.tsv')}
    new={(r['study_role'],r['lineage_group']):r for r in read_table(a.expanded/'lineage_coverage.tsv')}
    if set(old)!=set(new):raise ValueError('Different lineage universe')
    keys=sorted(new);den=np.array([int(new[k]['taxa']) for k in keys])
    if not np.array_equal(den,[int(old[k]['taxa']) for k in keys]) or (den<=0).any():raise ValueError('Different or empty denominators')
    before=np.array([int(old[k]['taxa_with_any_usable_marker']) for k in keys])
    any_marker=np.array([int(new[k]['taxa_with_any_usable_marker']) for k in keys])
    ten=np.array([int(new[k]['taxa_with_at_least10_usable_markers']) for k in keys])
    fifty=np.array([int(new[k]['taxa_with_at_least50_usable_markers']) for k in keys])
    bins=[den-any_marker,any_marker-ten,ten-fifty,fifty]
    if any((b<0).any() for b in bins) or not np.array_equal(sum(bins),den):raise ValueError('Invalid coverage depth partition')
    fig,axes=plt.subplots(1,2,figsize=(14,11),sharey=True)
    y=np.arange(len(keys));labels=[f'{k[1]}  (n={den[i]})' for i,k in enumerate(keys)]
    ax=axes[0]
    ax.hlines(y,100*before/den,100*any_marker/den,color='#b9c1c6',lw=2,zorder=1)
    ax.scatter(100*before/den,y,s=34,facecolors='white',edgecolors='#cc7049',label='Earlier snapshot',zorder=3)
    ax.scatter(100*any_marker/den,y,s=18,color='#153f56',label='Expanded snapshot',zorder=4)
    ax.set_yticks(y,labels,fontsize=9);ax.set_xlabel('Taxa with ≥1 usable paired marker (%)')
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2,frameon=False,fontsize=9)
    left=np.zeros(len(keys))
    for values,color,label in zip(bins,['#e1e5e8','#a7d8d5','#398f95','#153f56'],['0 markers','1–9 markers','10–49 markers','≥50 markers']):
        width=100*values/den;axes[1].barh(y,width,left=left,height=.7,color=color,label=label);left+=width
    axes[1].set_xlabel('Expanded snapshot: taxa by marker depth (%)')
    axes[1].legend(loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2,frameon=False,fontsize=9)
    boundary=sum(k[0]=='ingroup' for k in keys)-.5
    for ax in axes:
        ax.set_xlim(-3,103);ax.set_xticks([0,25,50,75,100]);ax.set_axisbelow(True);ax.grid(axis='x',color='#eceff1')
        ax.axhline(boundary,color='#687780',ls='--',lw=.8)
        for spine in ax.spines.values():spine.set_visible(False)
        ax.tick_params(axis='both',length=0)
    axes[0].invert_yaxis()
    fig.suptitle('Paired sequence–structure coverage across the full sampling design',fontsize=15,x=.58,y=.975)
    fig.text(.24,.025,'All 526 manifest entries retained. Dashed line separates fungal groups (above) from outgroups (below).\nUsable: confidence/coverage eligible in a marker with ≥4 eligible taxa; not proof of orthology or supported evolutionary change.',fontsize=9,color='#425563')
    fig.subplots_adjust(left=.24,right=.98,bottom=.09,top=.90,wspace=.14)
    a.output.mkdir(parents=True)
    for extension in ['svg','png','pdf']:fig.savefig(a.output/('paired_lineage_coverage.'+extension),dpi=180)
    plt.close(fig)
    result={'status':'complete_paired_lineage_coverage_figure','previous_receipt_sha256':sha(a.previous/'receipt.json'),'expanded_receipt_sha256':sha(a.expanded/'receipt.json'),'script_sha256':sha(Path(__file__)),'matplotlib_version':matplotlib.__version__,'groups':len(keys),'taxa':int(den.sum()),'previous_taxa_with_any_marker':int(before.sum()),'expanded_depth_counts':{label:int(v.sum()) for label,v in zip(['zero','one_to_nine','ten_to_fortynine','fifty_or_more'],bins)},'artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
