#!/usr/bin/env python3
"""Plot verified predictor disagreement alongside changing coverage denominators."""
import argparse,csv,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_ortholog_pair_guide_comparison import sha

p=argparse.ArgumentParser(description=__doc__)
for name in ['comparison','readback','output','figure-prefix']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();rp=a.comparison/'receipt.json';r=json.loads(rp.read_text());audit=json.loads(a.readback.read_text())
if r['status']!='complete_overlap_model_feature_comparison' or audit['status']!='passed_full_overlap_model_feature_comparison_readback' or audit['producer_receipt_sha256']!=sha(rp):raise ValueError('Unverified full-model comparison')
table=a.comparison/'threshold_summary.tsv'
if sha(table)!=r['artifacts'][table.name]:raise ValueError('Changed summary')
rows=list(csv.DictReader(table.open(),delimiter='\t'))
if rows!=audit['summary']:raise ValueError('Summary differs from readback')
if len(rows)!=12 or {(x['plddt_cutoff'],x['pae_cutoff']) for x in rows}!={(str(p),str(q)) for p in [0,70,90] for q in ['unfiltered',5,10,15]}:raise ValueError('Incomplete threshold grid')
a.output.mkdir(parents=True,exist_ok=False);a.figure_prefix.parent.mkdir(parents=True,exist_ok=True)
for suffix in ['.png','.pdf']:
 if a.figure_prefix.with_suffix(suffix).exists():raise FileExistsError('Fresh figure prefix required')
fig,axes=plt.subplots(1,3,figsize=(13,4.8));colors={0:'#777777',70:'#167D9A',90:'#B95728'};paes=['unfiltered','5','10','15'];plotted=[]
for cutoff in [0,70,90]:
 subset=[next(x for x in rows if x['plddt_cutoff']==str(cutoff) and x['pae_cutoff']==pae) for pae in paes]
 for ax,key,mult in [(axes[0],'pooled_state_mismatch_fraction',100),(axes[1],'compared',1),(axes[2],'retained_residues',.001)]:
  values=[float(x[key])*mult for x in subset];ax.plot(range(4),values,'o-',label=f'pLDDT ≥{cutoff}',color=colors[cutoff],lw=1.6,ms=4)
  for pae,value in zip(paes,values):plotted.append(dict(plddt_cutoff=cutoff,pae_cutoff=pae,metric=key,plotted_value=value))
for ax,title,ylabel in zip(axes,['Structural-state disagreement','Eligible model pairs','Retained residues'],['Pooled disagreement (%)','Pairs (of 643)','Residues (thousands)']):
 ax.set_title(title,fontsize=12);ax.set_xticks(range(4),['None','5','10','15']);ax.set_xlabel('Maximum context PAE cutoff');ax.set_ylabel(ylabel);ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False);ax.set_ylim(bottom=0)
axes[1].set_ylim(0,680)
fig.suptitle('Same-sequence AlphaFold–ESMFold comparisons',fontsize=15,y=.98)
fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.915),ncol=3,frameon=False)
fig.text(.03,.025,'Full encoded proteins; eligibility requires ≥50 residues and ≥50% coverage in both models.\nCohorts and retained residues change across thresholds. Disagreement is not an experimental error rate or evolutionary change.',fontsize=9)
fig.subplots_adjust(left=.065,right=.98,top=.77,bottom=.25,wspace=.35)
artifacts={}
for suffix in ['.png','.pdf']:
 path=a.figure_prefix.with_suffix(suffix);fig.savefig(path,dpi=180);artifacts[str(path)]=sha(path)
plt.close(fig)
path=a.output/'plotted_values.tsv'
with path.open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(plotted[0]),delimiter='\t');w.writeheader();w.writerows(plotted)
artifacts[str(path)]=sha(path)
result=dict(status='complete_overlap_model_feature_figure',producer_receipt_sha256=sha(rp),readback_sha256=sha(a.readback),source_table_sha256=sha(table),script_sha256=sha(__file__),plotted_values=len(plotted),artifacts=artifacts,scope='All twelve verified confidence alternatives with explicit pair/residue denominators. Lines connect sensitivity settings, not a calibrated response to confidence; no causal or accuracy claim.')
(a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
