"""Plot fully accounted synthetic coverage after the audited numerical overlay."""
import json
from pathlib import Path
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha,write_json
from matched_calibration_intervals import exact_coverage_bounds


def main():
    root=Path('results/model_validation/matched-kr-continuation-audit-20260928-v1')
    receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='passed_all_184_qualification_overlays_and_15984_denominator_readback'
    assert receipt['attempted']==15984 and receipt['unresolved_coefficient_outcomes']==0
    path=root/'all_coefficient_coverage.tsv';assert sha(path)==receipt['artifacts'][path.name]
    table=pd.read_csv(path,sep='\t')
    assert len(table)==112 and (table.unresolved==0).all() and (table.attempted==999).all()
    for row in table.itertuples(index=False):
        expected=exact_coverage_bounds(row.covered,row.unresolved,row.attempted,row.alpha)
        for field in ['lower','upper','observed_fraction_lower','observed_fraction_upper']:
            np.testing.assert_allclose(getattr(row,field),expected[field],rtol=1e-12,atol=1e-14)
    focal=table[table.coefficient==0]
    assert len(focal)==32 and focal.lower.min()>=.8 and focal.upper.max()<=1
    out=Path('results/figures/matched-kr-coverage-resolved-20260928-v1');out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(path,out/path.name)
    fig,axes=plt.subplots(1,2,figsize=(11,6.4),sharex=True,sharey=True)
    for ax,n in zip(axes,[48,240]):
        labels=[]
        for index in range(8):
            subset=focal[focal.case==f'n{n}-variance-{index}'];assert len(subset)==2
            first=subset.iloc[0]
            labels.append(', '.join(f'{first[k]:g}' for k in ['background_ratio','family_ratio','species_ratio']))
            for method,color,offset,marker in [('KR candidate','#167d9a',-.14,'o'),('Conditional t','#a74c39',.14,'s')]:
                r=subset[subset.method==method].iloc[0];point=r.covered/r.attempted
                ax.errorbar(point,index+offset,xerr=[[point-r.lower],[r.upper-point]],fmt=marker,
                            color=color,ms=4,capsize=2,lw=1.2,label=method if index==0 else None)
        ax.axvline(.95,color='#555555',ls='--',lw=1,label='Nominal 95%')
        ax.set_yticks(range(8),labels);ax.set_xlim(.8,1.002);ax.set_ylim(7.6,-.6)
        ax.set_xticks([.8,.85,.9,.95,1.],['80%','85%','90%','95%','100%'])
        ax.set_title(f'{n} observations');ax.set_xlabel('Intercept coverage (axis starts at 80%)')
        ax.grid(axis='x',color='#eeeeee');ax.set_axisbelow(True);ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Generating variance ratios\n(background, family, species)')
    fig.suptitle('Synthetic interval coverage after numerical qualification',fontsize=14)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.52,.12),ncol=3,frameon=False)
    fig.text(.52,.055,'999 fixed replicates per configuration; all 15,984 responses retained. Bars: marginal 95% Monte Carlo intervals.\n184 reviews resolved by verified optimizer continuations; selected estimates unchanged.\nSynthetic Gaussian designs only; no simultaneous or fungal-data coverage guarantee.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.20,1,.94));fig.savefig(out/'intercept_coverage.png',dpi=180);fig.savefig(out/'intercept_coverage.pdf');plt.close(fig)
    result=dict(status='completed_resolved_synthetic_coverage_figure',source_receipt_sha256=sha(root/'receipt.json'),
        source_table_sha256=sha(path),script_sha256=sha(__file__),attempted=15984,table_rows=112,plotted_intercept_rows=32,
        axis_limits=[.8,1.002],all_plotted_intervals_within_axis=True,
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='Updated figure after audited qualification overlay; original unresolved figure preserved. '
              'All coefficients retained in copied table; figure shows intercept only. '
              'Marginal Monte Carlo uncertainty, not multiplicity or biological inference. Visual inspection pending.')
    write_json(out/'receipt.json',result);write_json(Path('metadata/matched_kr_resolved_figure_20260928.json'),result)
    print(result['status'])


if __name__=='__main__':main()
