#!/usr/bin/env python3
"""Plot every paired FCS difference interval under both copy-review policies."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['full','omission','readback','output']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    audit=json.loads(a.readback.read_text())
    if audit['status']!='passed_all_saved_difference_arrays_and_sorted_percentiles':
        raise ValueError('Numerical readback required')
    sources={};pins={str(a.readback):sha(a.readback)}
    for label,folder,count in [('full',a.full,89),('omission',a.omission,88)]:
        r=json.loads((folder/'receipt.json').read_text())
        if (r['status']!='complete_paired_conditional_fcs_bootstrap_comparison'
                or r['markers']!=count or r['contrasts']!=120
                or audit['input_receipts'].get(str(folder/'receipt.json'))!=sha(folder/'receipt.json')):
            raise ValueError('Wrong or unvalidated figure source')
        path=folder/'paired_difference_summary.tsv'
        if sha(path)!=r['artifacts'][path.name]:raise ValueError('Changed summary')
        pins[str(folder/'receipt.json')]=sha(folder/'receipt.json');pins[str(path)]=sha(path)
        sources[label]=pd.read_csv(path,sep='\t',float_precision='round_trip').set_index(['model_id','contrast'],verify_integrity=True)
    grid=list(itertools.product(['3di_af','3di_af_empirical','3di_llm'],['gamma','freerate'],['tien2013','miller1987'],['coverage_confidence','composition_adjusted']))
    contrasts=[('aa_at_RSA_025','AA-rate association\nat RSA = 0.25'),
               ('aa_at_baseline_mean_RSA','AA-rate association\nat baseline mean RSA'),
               ('rsa_at_aa_log_rate_0','RSA association\nat log(1 + AA rate) = 0'),
               ('rsa_at_baseline_mean_aa_log_rate','RSA association\nat baseline mean log(1 + AA rate)'),
               ('interaction','AA-rate × RSA\ninteraction')]
    colors={'full':'#156B8A','omission':'#BA5726'}
    labels=dict(zip(['3di_af','3di_af_empirical','3di_llm','gamma','freerate','tien2013','miller1987','coverage_confidence','composition_adjusted'],
                    ['3Di-AF','3Di-AF+F','3Di-LLM','G4','R4','Tien','Miller','Base','+Comp']))
    fig,axes=plt.subplots(1,5,figsize=(19,10.5),sharey=True)
    plotted=[]
    for ax,(contrast,title) in zip(axes,contrasts):
        ax.axvline(0,color='#666666',ls='--',lw=.85)
        for start in [0,16]:ax.axhspan(start-.5,start+7.5,color='#F0F3F5',zorder=0)
        for y,spec in enumerate(grid):
            model='__'.join(spec)
            for group,offset,marker in [('full',-.14,'o'),('omission',.14,'s')]:
                row=sources[group].loc[(model,contrast)]
                values=[float(row[k]) for k in ['sensitivity_minus_baseline','percentile_lower','percentile_upper']]
                point,low,high=values
                if not np.isfinite(values).all() or low>high:raise ValueError('Invalid plotting values')
                ax.plot([low,high],[y+offset]*2,color=colors[group],lw=1)
                ax.plot(point,y+offset,marker,color=colors[group],ms=3)
                plotted.append(dict(cohort=group,model_id=model,contrast=contrast,point=point,lower=low,upper=high))
        ax.set_title(title,fontsize=10,pad=12)
        ax.set_xlabel('FCS omission − baseline',fontsize=9)
        ax.tick_params(axis='x',labelsize=8)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
        ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0),useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(8)
        ax.grid(axis='y',alpha=.15)
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_yticks(range(24),[' | '.join(labels[x] for x in spec) for spec in grid],fontsize=8)
    axes[0].set_ylim(23.7,-.7)
    fig.suptitle('FCS omission sensitivity: paired changes in conditional coupling estimates',fontsize=15,y=.978)
    fig.legend(handles=[Line2D([0],[0],color=colors['full'],marker='o',label='Full cohort: 89 markers'),
                        Line2D([0],[0],color=colors['omission'],marker='s',label='Copy-review omission: 88 markers')],
               loc='upper center',bbox_to_anchor=(.59,.94),ncol=2,frameon=False,fontsize=10)
    notes=[
        'Points: fitted differences. Lines: unadjusted 95% percentile intervals from 2,000 paired whole-marker draws. All 24 specifications and five contrasts shown.',
        "Each baseline/FCS pair uses its cohort's fixed baseline references; panel scales differ. Base: coverage/confidence controls; +Comp: amino-acid composition controls.",
        'Copy-review omission excludes marker 4986044at2759. Structural alphabet labels are fitted models; all predictions in these cohorts are ESMFold.',
        'Intervals are conditional on estimated rates and fixed gene trees. No multiple-testing significance, equivalence, causal effect or structural-acceleration claim.'
    ]
    for i,line in enumerate(notes):fig.text(.02,.087-i*.022,line,fontsize=9)
    fig.subplots_adjust(left=.225,right=.988,top=.85,bottom=.155,wspace=.27)
    a.output.mkdir(parents=True)
    table=a.output/'plotted_intervals.tsv';pd.DataFrame(plotted).to_csv(table,sep='\t',index=False)
    for row in pd.read_csv(table,sep='\t',float_precision='round_trip').itertuples():
        source=sources[row.cohort].loc[(row.model_id,row.contrast)]
        np.testing.assert_allclose([row.point,row.lower,row.upper],source[['sensitivity_minus_baseline','percentile_lower','percentile_upper']].to_numpy(dtype=float),rtol=0,atol=1e-15)
    if len(plotted)!=240:raise ValueError('Incomplete figure grid')
    fig.savefig(a.output/'paired_fcs_differences.svg',metadata={'Date':None})
    fig.savefig(a.output/'paired_fcs_differences.png',dpi=160)
    fig.savefig(a.output/'paired_fcs_differences.pdf',metadata={'CreationDate':None,'ModDate':None})
    plt.close(fig)
    for path,digest in pins.items():
        if sha(path)!=digest:raise ValueError('Figure source changed')
    result=dict(status='complete_paired_fcs_difference_figure',plotted_intervals=240,source_pins=pins,script_sha256=sha(__file__),
                artifacts={p.name:sha(p) for p in a.output.iterdir()},scope='Every matched specification and contrast displayed; exported points/endpoints read back against validated source tables. Visual inspection remains separate. Unadjusted conditional intervals, no new significance tests.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
