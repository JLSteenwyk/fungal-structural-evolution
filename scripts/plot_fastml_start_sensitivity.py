#!/usr/bin/env python3
"""Plot the fully checked start-count and internal-posterior range diagnostics."""
import argparse
from collections import Counter
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from screen_duplication_alignment_reuse import sha


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--proof',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    proof=json.loads(args.proof.read_text())
    assert proof['status']=='passed_full_fastml_start_posterior_sensitivity_readback'
    for path,h in proof['source_hashes'].items():assert sha(path)==h,path
    table=next(Path(p) for p in proof['source_hashes'] if Path(p).name=='start_sensitivity.jsonl')
    with table.open() as stream:rows=[json.loads(line) for line in stream]
    assert len(rows)==proof['groups']==156
    rows=[r for r in rows if r['status']!='no_coded_characters']
    assert len(rows)==proof['nonempty_groups']==153
    counts=Counter(r['near_best_starts'] for r in rows)
    assert {str(k):v for k,v in counts.items()}==proof['near_best_start_counts']
    singleton=[r for r in rows if r['near_best_starts']==1]
    multiple=[r for r in rows if r['near_best_starts']>1]
    assert all(r['near_best_maximum_internal_probability_range']==0 for r in singleton)
    args.output.mkdir(parents=True,exist_ok=False)
    fig,axes=plt.subplots(1,2,figsize=(12,5.2),gridspec_kw={'width_ratios':[1,1.35]})
    positions=range(1,6);values=[counts[n] for n in positions]
    axes[0].bar(positions,values,color=['#9a9fa5']+['#247b91']*4,width=.7)
    for x,y in zip(positions,values):axes[0].text(x,y+1.5,str(y),ha='center',fontsize=11)
    axes[0].set_xticks(list(positions));axes[0].set_ylim(0,max(values)*1.2)
    axes[0].set_xlabel('Starts within 10⁻⁵ log-likelihood units of best')
    axes[0].set_ylabel('Input groups');axes[0].set_title('Retained starts per input')
    for selected,marker,color,label in [(multiple,'o','#247b91',f'2–5 near-best starts ({len(multiple)})'),
                                        (singleton,'s','#777d84',f'1 near-best start ({len(singleton)})')]:
        x=[100*r['all_maximum_internal_probability_range'] for r in selected]
        y=[100*r['near_best_maximum_internal_probability_range'] for r in selected]
        axes[1].scatter(x,y,s=28,marker=marker,edgecolors=color,
                        facecolors=color if marker=='o' else 'none',alpha=.75,label=label)
    axes[1].set_xlabel('All-start maximum spread (percentage points)')
    axes[1].set_ylabel('Near-best maximum spread (percentage points)')
    axes[1].set_title('Internal-node presence-probability spread')
    ymax=100*max(r['near_best_maximum_internal_probability_range'] for r in rows)
    axes[1].set_ylim(-.015,max(.3,ymax*1.15));axes[1].set_yticks([0,.1,.2,.3])
    axes[1].legend(loc='upper right',fontsize=9,frameon=False)
    for axis in axes:
        axis.spines[['top','right']].set_visible(False)
        axis.grid(axis='y',alpha=.18);axis.set_axisbelow(True)
    fig.suptitle('Ancestral presence probabilities across optimizer starts',fontsize=14)
    fig.text(.5,.02,'153 inputs; five starts each. One retained start has zero spread automatically. Ranges are not confidence intervals.',
             ha='center',fontsize=9,color='#444444')
    fig.tight_layout(rect=(0,.06,1,.94))
    for suffix in ['png','pdf']:
        fig.savefig(args.output/('fastml_start_sensitivity.'+suffix),dpi=180,bbox_inches='tight')
    plt.close(fig)
    receipt=dict(status='complete_source_checked_fastml_start_sensitivity_figure_pending_visual_review',
        source_hashes={str(args.proof):sha(args.proof),str(table):sha(table),str(Path(__file__)):sha(__file__)},
        groups=len(rows),near_best_start_counts={str(k):v for k,v in counts.items()},
        plotted_rows=[dict(job_id=r['job_id'],near_best_starts=r['near_best_starts'],
                          all_maximum_internal_probability_range=r['all_maximum_internal_probability_range'],
                          near_best_maximum_internal_probability_range=r['near_best_maximum_internal_probability_range']) for r in rows],
        artifacts={p.name:sha(p) for p in args.output.iterdir()},
        scope='All 153 nonempty inputs; raw probability spreads converted to percentage points. Single retained starts distinguished. No independence, confidence interval, globally converged optimizer or biologically qualified ancestral estimate inferred.')
    (args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
