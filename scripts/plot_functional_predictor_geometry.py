"""Plot every qualified cross-predictor functional-context comparison."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha

root=Path('results/functional_sites/prediction-context-geometry-20260927-v1')
rp=root/'receipt.json';receipt=json.loads(rp.read_text())
auditpath=Path('metadata/functional_context_geometry_summary_readback_20260927.json')
audit=json.loads(auditpath.read_text())
assert audit['status']=='passed_complete_functional_context_geometry_summary_readback'
assert sha(rp)==audit['source_receipt_sha256']
table=root/'functional_context_geometry.tsv'
assert sha(table)==receipt['artifacts'][table.name]==audit['source_table_sha256']
summarypath=Path('metadata/functional_context_geometry_summary_20260927.tsv')
assert sha(summarypath)==audit['output_sha256']
data=pd.read_csv(table,sep='\t');summary=pd.read_csv(summarypath,sep='\t')
keys=['marker','taxon_id','protein_id','protein_residue_1based']
assert len(data)==600 and len(data[keys].drop_duplicates())==150
assert not data.duplicated(keys+['context_definition']).any()
definitions=['focal_triplet','afdb_context','esmfold_context','context_union']
titles=['Focal residue and sequence neighbors','AlphaFold-selected context','ESMFold-selected context','Union of both contexts']
groups=[('same','same'),('same','different'),('different','same'),('different','different')]
labels=['Same state\nsame partner','Same state\ndifferent partner','Different state\nsame partner','Different state\ndifferent partner']
colors=['#4477AA','#66CCEE','#CCBB44','#AA3377']
fig,axes=plt.subplots(2,2,figsize=(11,8),sharey=True)
rng=np.random.default_rng(927);total=0
for ax,definition,title in zip(axes.flat,definitions,titles):
    panel=data[data.context_definition.eq(definition)]
    assert len(panel)==150
    for i,((state,partner),color) in enumerate(zip(groups,colors)):
        part=panel[panel.paired_state_comparison.eq(state)&panel.partner_comparison.eq(partner)].sort_values(keys)
        values=part.ca_rmsd_angstrom.to_numpy();assert len(values)>0
        sr=summary[summary.context_definition.eq(definition)&summary.paired_state_comparison.eq(state)&summary.partner_comparison.eq(partner)].iloc[0]
        assert len(values)==sr.residues
        np.testing.assert_allclose([min(values),np.median(values),max(values)],[sr.rmsd_min,sr.rmsd_median,sr.rmsd_max],rtol=1e-12,atol=1e-14)
        ax.scatter(i+rng.uniform(-.17,.17,len(values)),values,s=19,c=color,alpha=.65,linewidths=0)
        ax.plot([i-.23,i+.23],[np.median(values)]*2,color='black',lw=2)
        ax.text(i,.47,f'n = {len(values)}',ha='center',fontsize=9)
        total+=len(values)
    ax.set(title=title,ylim=(-.015,.505),xticks=range(4),xticklabels=labels)
    ax.tick_params(axis='x',labelsize=9);ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
for ax in axes[:,0]: ax.set_ylabel('Cα RMSD between predictors (Å)')
fig.suptitle('Local geometry across AlphaFold and ESMFold predictions',fontsize=15,y=.99)
fig.text(.5,.025,'150 exact protein positions; each appears in all four panels. Black bars: medians.\nState = structural-alphabet state; partner = descriptor-selected residue. Points are dependent observations.',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.075,1,.955))
assert total==600
prefix=Path('docs/figures/functional_predictor_geometry_20260927');prefix.parent.mkdir(exist_ok=True,parents=True)
artifacts={}
for ext in ['png','pdf','svg']:
    path=prefix.with_suffix('.'+ext);assert not path.exists();fig.savefig(path,dpi=180);artifacts[str(path)]=sha(path)
plt.close(fig)
result=dict(status='complete_full_functional_predictor_geometry_figure',plotted_comparisons=total,unique_focal_positions=150,panels=4,verified_summary_groups=16,
    source_receipt_sha256=sha(rp),source_table_sha256=sha(table),summary_audit_sha256=sha(auditpath),script_sha256=sha(__file__),artifacts=artifacts,
    scope='All qualified observations shown without additional filtering. Four context definitions reuse the same positions. Medians are descriptive; no independent replication, predictor accuracy, causal effect, or evolutionary change is established.')
prefix.with_suffix('.receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
