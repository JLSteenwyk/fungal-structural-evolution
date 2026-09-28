"""Plot empirical residual fourth moments for every synthetic configuration."""
import json
from pathlib import Path
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha,write_json

root=Path('results/model_validation/matched-simulation-residual-reference-20260928-v1')
r=json.loads((root/'receipt.json').read_text())
audit_path=Path('metadata/matched_simulation_residual_reference_readback_20260928.json')
a=json.loads(audit_path.read_text())
assert a['status']=='passed_all_residual_reference_summary_arithmetic'
assert a['source_receipt_sha256']==sha(root/'receipt.json') and a['unresolved']==0
source=root/'diagnostic_reference_summary.tsv'
assert sha(source)==r['artifacts'][source.name]
frame=pd.read_csv(source,sep='\t');table=frame[frame.metric=='fourth_raw_moment']
assert len(table)==32 and (table.available==999).all()
assert (table.q025<=table['median']).all() and (table['median']<=table.q975).all()
plan=json.loads(Path('metadata/matched_kr_simulations_plan_20260928.json').read_text())
cases={c['id']:c for c in plan['cases']}
xmax=float(np.ceil(table.q975.max()*1.06))
out=Path('results/figures/matched-residual-reference-20260928-v1');out.mkdir(parents=True,exist_ok=False)
shutil.copyfile(source,out/source.name)
fig,axes=plt.subplots(1,2,figsize=(11,6.5),sharex=True,sharey=True)
for ax,n in zip(axes,[48,240]):
    labels=[]
    for index in range(8):
        name=f'n{n}-variance-{index}';labels.append(', '.join(f'{v:g}' for v in cases[name]['ratios']))
        for mode,color,offset,marker,label in [('generating_covariance','#287d9c',-.15,'o','Generating covariance'),('fitted_covariance','#b65e39',.15,'s','Fitted covariance')]:
            row=table[(table.case==name)&(table['mode']==mode)].iloc[0]
            ax.errorbar(row['median'],index+offset,xerr=[[row['median']-row.q025],[row.q975-row['median']]],
                fmt=marker,color=color,capsize=2,ms=4,lw=1.2,label=label if index==0 else None)
    ax.axvline(3,color='#777777',ls='--',lw=1,label='Standard normal fourth moment')
    ax.set_title(f'{n} observations');ax.set_yticks(range(8),labels)
    ax.set_ylim(7.6,-.6);ax.set_xlim(0,xmax);ax.set_xlabel('Mean of standardized residuals raised to the fourth power')
    ax.grid(axis='x',color='#ededed');ax.set_axisbelow(True);ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('Generating variance ratios\n(background, family, species)')
fig.suptitle('Residual diagnostic spread depends on covariance fitting',fontsize=14)
h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='lower center',bbox_to_anchor=(.53,.12),ncol=3,frameon=False,fontsize=9)
fig.text(.52,.055,'Points: empirical medians. Bars: empirical 2.5%–97.5% ranges across 999 responses per configuration.\nSame 15,984 Gaussian responses in both modes; fixed effects estimated in both. Residual correlations remain.\nSynthetic design references only: ranges are not confidence intervals or fungal-data rejection thresholds.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.21,1,.94))
fig.savefig(out/'residual_fourth_moment.png',dpi=180);fig.savefig(out/'residual_fourth_moment.pdf');plt.close(fig)
receipt=dict(status='completed_residual_reference_figure_pending_visual_inspection',source_table_sha256=sha(source),
    source_readback_sha256=sha(audit_path),script_sha256=sha(__file__),plotted_points=32,source_summary_rows=len(frame),
    axis_limits=[0,xmax],all_ranges_within_axes=bool(table.q025.min()>=0 and table.q975.max()<=xmax),
    artifacts={p.name:sha(p) for p in out.iterdir()},scope='Descriptive synthetic fourth-moment summaries; complete256-row source table retained. No calibrated hypothesis test.')
assert receipt['all_ranges_within_axes']
write_json(out/'receipt.json',receipt);write_json(Path('metadata/matched_residual_reference_figure_20260928.json'),receipt)
print(receipt['status'])
