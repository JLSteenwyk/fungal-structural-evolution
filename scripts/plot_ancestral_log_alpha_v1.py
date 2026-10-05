#!/usr/bin/env python3
"""Plot complete diagnostic alpha traces without treating them as posterior draws."""
import argparse
from collections import defaultdict
import csv
from datetime import datetime,timezone
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    args.output.mkdir(parents=True)
    pp=Path('metadata/ancestral_log_alpha_diagnostic_20261005_v2.json')
    rp=Path('metadata/ancestral_log_alpha_readback_20261005_v1.json')
    tp=Path('metadata/ancestral_log_alpha_readback_transport_20261005_v1.json')
    producer,reader,transport=[json.loads(p.read_text()) for p in [pp,rp,tp]]
    assert reader['producer_receipt_sha256']==sha(pp)
    assert reader['status']=='passed_independent_full_current_alpha_trace_and_native_decimal_readback'
    assert transport['validation_sha256']==sha(rp) and transport['original_tool_terminal_exit_code']==0
    pins={}
    for path in [pp,rp,tp,Path(__file__)]:bind(pins,path)
    trace=Path('results/ancestral/full-current-log-alpha-diagnostic-20261005-v2/all_current_alpha_traces.jsonl')
    bind(pins,trace,producer['source_hashes'][str(trace)])
    groups=defaultdict(list);total=0
    with trace.open() as handle:
        for line in handle:
            row=json.loads(line);groups[(row['stage'],row['prior_label'],row['iteration'])].append(row);total+=1
    assert total==34020 and len(groups)==126
    summary=[]
    for (stage,prior,iteration),rows in sorted(groups.items()):
        values=[np.log10(row['alpha']) for row in rows if row['alpha_state']=='finite']
        missing=sum(row['alpha_state']=='positive_infinity' for row in rows)
        assert len(values)+missing==len(rows)
        assert len(rows)==(532 if stage=='v6' else 8)
        q=np.quantile(values,[.1,.5,.9],method='linear')
        summary.append(dict(stage=stage,prior=prior,iteration=iteration,role_rows=len(rows),
            finite_rows=len(values),positive_infinity_rows=missing,
            finite_log10_alpha_q10=float(q[0]),finite_log10_alpha_median=float(q[1]),
            finite_log10_alpha_q90=float(q[2])))
    table=args.output/'alpha_trace_summary.tsv'
    with table.open('x') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(summary[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(summary)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(11.5,6.8),sharex=True,sharey='row',constrained_layout=True)
    for column,prior in enumerate(['package','centered','broad']):
        subset=[row for row in summary if row['stage']=='v6' and row['prior']==prior]
        x=[row['iteration'] for row in subset]
        axes[0,column].fill_between(x,[row['finite_log10_alpha_q10'] for row in subset],
            [row['finite_log10_alpha_q90'] for row in subset],alpha=.22,color='#277DA1')
        axes[0,column].plot(x,[row['finite_log10_alpha_median'] for row in subset],color='#277DA1',lw=2)
        axes[0,column].set_title(prior.capitalize()+' prior')
        for stage,color,label in [('v6','#277DA1','Full V6: 532 roles per prior'),
                                  ('v7_comparison','#D1495B','V7 comparison: 8 roles per prior')]:
            rows=[row for row in summary if row['stage']==stage and row['prior']==prior]
            axes[1,column].plot([row['iteration'] for row in rows],
                [100*row['positive_infinity_rows']/row['role_rows'] for row in rows],
                marker='o',ms=3,color=color,label=label)
        axes[1,column].set_ylim(-.3,14)
        axes[1,column].set_xlabel('Saved scalar iteration')
        axes[1,column].set_xticks([0,5,10,15,20])
        for ax in axes[:,column]:ax.grid(alpha=.2)
    axes[0,0].set_ylabel('log10(alpha), finite values only')
    axes[1,0].set_ylabel('Roles logging positive infinity (%)')
    axes[1,2].legend(loc='upper right',fontsize=8,frameon=False)
    fig.suptitle('Short-run alpha diagnostics: full design, not qualified posterior samples',fontsize=13)
    fig.supxlabel('Top: median and central 80% of finite full-V6 role values. Bottom: nonfinite values remain explicit.\n'
        'Families and comparison roles are not independent biological replicates.',fontsize=9)
    outputs=[]
    for suffix in ['png','pdf']:
        path=args.output/('ancestral_log_alpha_diagnostics_20261005.'+suffix)
        fig.savefig(path,dpi=180);outputs.append(path);bind(pins,path)
    plt.close(fig)
    bind(pins,table)
    verify(pins)
    result=dict(status='completed_full_current_alpha_trace_diagnostic_figure',checked_utc=datetime.now(timezone.utc).isoformat(),
        scalar_rows=total,summary_bins=126,plot_full_v6_native_zero_roles=1596,
        separate_v7_comparison_roles=24,positive_infinity_observations=26,outputs=[str(p) for p in outputs],
        table=str(table),source_hashes=pins,scientific_eligibility=False,new_mcmc_runs=0,gpu=False,
        scope='All34,020readback-verified rows summarized by original phase/prior/iteration.126bins retain '
              'every finite and infinite alpha observation. Finite q10/median/q90 describe532full-V6roles '
              'per prior; separate8-role V7comparison has its own denominator. Failed originals are '
              'not silently replaced. No pooled posterior, confidence interval, independent biological '
              'replicate, convergence diagnostic, review admission or biological effect claim.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
