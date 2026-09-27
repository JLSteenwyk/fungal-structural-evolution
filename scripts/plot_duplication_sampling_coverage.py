#!/usr/bin/env python3
"""Plot verified frozen structural coverage with explicit event and taxon denominators."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from run_ortholog_pair_guide_comparison import sha

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--membership',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--figure-prefix',type=Path,required=True);a=p.parse_args()
r=json.loads((a.membership/'receipt.json').read_text());table=a.membership/'taxon_coverage_with_membership.tsv'
if r['status']!='complete_duplication_coverage_reconciliation_membership' or sha(table)!=r['artifacts'][table.name]:raise ValueError('Unverified membership table')
for path,h in r['source_sha256'].items():
    if sha(path)!=h:raise ValueError('Changed membership source')
data=pd.read_csv(table,sep='\t',keep_default_na=False);included=data[data.in_reconciliation.eq(1)].copy()
if set(included.guide)!= {'profile','mafft'} or included.groupby('guide').taxon.nunique().to_dict()!= {'mafft':526,'profile':526}:raise ValueError('Unexpected reconciliation taxon scope')
included['lineage_group']=included.lineage.str.split(';').str[0]
included['has_both_models']=included.both_models.gt(0).astype(int)
group=included.groupby(['guide','study_role','lineage_group'],sort=True).agg(taxa=('taxon','size'),taxa_with_both_models=('has_both_models','sum'),terminal_events=('terminal_singleton_events','sum'),neither_model=('neither_model','sum'),one_model=('one_model','sum'),both_models=('both_models','sum')).reset_index()
if not (group.terminal_events==group.neither_model+group.one_model+group.both_models).all():raise ValueError('Event partition differs')
group['event_coverage_pct']=100*group.both_models/group.terminal_events
group['taxon_coverage_pct']=100*group.taxa_with_both_models/group.taxa
for guide,part in included.groupby('guide'):
    observed=group[group.guide.eq(guide)]
    if int(observed.terminal_events.sum())!=int(part.terminal_singleton_events.sum()) or int(observed.both_models.sum())!=int(part.both_models.sum()):raise ValueError('Lineage aggregation differs')
a.output.mkdir(parents=True,exist_ok=False);a.figure_prefix.parent.mkdir(parents=True,exist_ok=True)
for suffix in ['.png','.pdf','.svg']:
    if a.figure_prefix.with_suffix(suffix).exists():raise FileExistsError('Fresh figure prefix required')
group.to_csv(a.output/'lineage_coverage.tsv',sep='\t',index=False)
order=group[group.guide.eq('profile')].sort_values(['study_role','taxa','lineage_group'],ascending=[True,False,True])
keys=list(zip(order.study_role,order.lineage_group));labels=[f'{row.lineage_group}  (n={row.taxa})'+(' *' if row.study_role=='outgroup' else '') for row in order.itertuples()]
fig,axes=plt.subplots(1,2,figsize=(12,11),sharey=True)
colors={'profile':'#176B87','mafft':'#C86020'}
for guide,marker,offset in [('profile','o',-.09),('mafft','x',.09)]:
    part=group[group.guide.eq(guide)].set_index(['study_role','lineage_group']).loc[keys]
    for ax,column in zip(axes,['event_coverage_pct','taxon_coverage_pct']):
        ax.scatter(part[column], [i+offset for i in range(len(keys))],s=27,marker=marker,color=colors[guide],label=guide.upper(),zorder=3)
for ax,title in zip(axes,['Events with both copies modeled','Taxa contributing ≥1 modeled pair']):
    ax.set_title(title,fontsize=12,pad=12);ax.set_xlim(-2,102);ax.set_xticks([0,25,50,75,100]);ax.set_xlabel('Coverage (%)');ax.grid(axis='x',alpha=.25);ax.spines[['top','right']].set_visible(False)
axes[0].set_yticks(range(len(keys)),labels,fontsize=9);axes[0].invert_yaxis();fig.legend(*axes[1].get_legend_handles_labels(),loc='upper right',bbox_to_anchor=(.98,.957),ncol=2,frameon=False,fontsize=9)
fig.suptitle('Structural coverage of reported terminal duplication candidates',fontsize=15,y=.98)
fig.text(.29,.937,'526 reconciled taxa; frozen AlphaFold model inventory',fontsize=10)
fig.text(.02,.02,'n = reconciled taxa in each manifest lineage group; * = outgroup. Candidate-manifest exclusion S. jurei omitted.\nEvent percentages use all reported terminal singleton-side events as denominators. Missing models are not biological absences.\nTree guides overlap and are sensitivity analyses, not independent replicates. These summaries do not correct ascertainment bias.',fontsize=9)
fig.subplots_adjust(left=.29,right=.98,bottom=.12,top=.91,wspace=.12)
paths=[]
for suffix in ['.png','.pdf','.svg']:
    path=a.figure_prefix.with_suffix(suffix);fig.savefig(path,dpi=180);paths.append(path)
plt.close(fig)
result=dict(status='complete_verified_duplication_lineage_coverage_figure',membership_receipt_sha256=sha(a.membership/'receipt.json'),source_table_sha256=sha(table),lineage_rows=len(group),reconciled_taxa_per_guide=526,excluded_candidate_taxa=sorted(set(data.loc[data.in_reconciliation.eq(0),'taxon'])),artifacts={str(path):sha(path) for path in [a.output/'lineage_coverage.tsv',*paths]},scope='Descriptive broad-lineage aggregation of independently checked terminal-event coverage. Separate event-weighted and taxon-level availability denominators; no phylogenetic correction, biological absence, missingness adjustment or independent-replicate claim.')
(a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(group[['guide','study_role','lineage_group','taxa','both_models','event_coverage_pct']].to_string(index=False))
