#!/usr/bin/env python3
"""Plot verified control availability and taxonomic concentration, not effect estimates."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from run_ortholog_pair_guide_comparison import sha
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--proof',type=Path,required=True);p.add_argument('--prefix',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
rp=a.source/'receipt.json';r=json.loads(rp.read_text());proof=json.loads(a.proof.read_text());f=a.source/'taxon_support.tsv'
assert proof['status']=='passed_full_taxon_architecture_support_readback' and proof['producer_receipt_sha256']==sha(rp) and sha(f)==r['artifacts'][f.name]
bindings={str(p):sha(p) for p in [rp,a.proof,f]}
d=pd.read_csv(f,sep='\t');d=d[d.background_set.eq('both_guides_unreported_parents')]
policies=['alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore'];labels=['Alignment\nE-value','Alignment\nbit score','Envelope\nE-value','Envelope\nbit score'];guides=['profile','mafft'];colors=['#126782','#b04c15']
coverage=[];curves=[]
for guide in guides:
 for policy in policies:
  x=d[d.guide.eq(guide)&d.policy.eq(policy)];total=int(x.targets.sum())
  for metric in ['supported_1_5','focal_supported_1_5']:
   n=int(x[metric].sum());coverage.append(dict(guide=guide,policy=policy,metric=metric,numerator=n,denominator=total,percent=100*n/total))
 x=d[d.guide.eq(guide)&d.policy.eq('alignment_evalue')].sort_values(['supported_1_5','taxon_id'],ascending=[False,True]);total=int(x.supported_1_5.sum());cumulative=0
 for rank,row in enumerate(x.itertuples(),1):
  cumulative+=int(row.supported_1_5);curves.append(dict(guide=guide,rank=rank,taxon_id=row.taxon_id,supported_targets=int(row.supported_1_5),cumulative_supported=cumulative,total_supported=total,percent=100*cumulative/total))
c=pd.DataFrame(coverage);q=pd.DataFrame(curves);a.prefix.parent.mkdir(parents=True,exist_ok=True)
outputs=[]
for suffix,data in [('coverage.tsv',c),('concentration.tsv',q)]:
 path=Path(str(a.prefix)+'_'+suffix);data.to_csv(path,sep='\t',index=False);outputs.append(path)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
fig,ax=plt.subplots(1,2,figsize=(12,5.4),gridspec_kw={'width_ratios':[1.15,1]})
for i,(guide,color) in enumerate(zip(guides,colors)):
 for metric,marker,face in [('supported_1_5','o',color),('focal_supported_1_5','s','white')]:
  v=c[c.guide.eq(guide)&c.metric.eq(metric)].set_index('policy').loc[policies]
  ax[0].scatter(np.arange(4)+(i-.5)*.13,v.percent,s=52,marker=marker,facecolor=face,edgecolor=color,label=guide+' — '+('any taxon' if metric=='supported_1_5' else 'focal taxon'),zorder=3)
 ax[1].plot([0]+q[q.guide.eq(guide)]['rank'].tolist(),[0]+q[q.guide.eq(guide)].percent.tolist(),color=color,linestyle='-' if guide=='profile' else '--',linewidth=2,label=guide)
ax[0].set(title='A  Available controls',ylabel='Targets with ≥1 eligible background (%)',xticks=np.arange(4),xticklabels=labels,ylim=(0,21));ax[0].grid(axis='y',alpha=.2);ax[0].legend(loc='center right',frameon=False,fontsize=9)
n=q[q.guide.eq('profile')]['rank'].max();ax[1].plot([0,n],[0,100],color='#999999',linestyle=':',label='Equal shares across observed taxa')
five=q[q.guide.eq('profile')&q['rank'].eq(5)].iloc[0];ax[1].scatter([5],[five.percent],color=colors[0],s=25,zorder=5);ax[1].annotate(f'Five taxa: {five.percent:.1f}%',xy=(5,five.percent),xytext=(30,27),arrowprops={'arrowstyle':'-','color':'#555555'},fontsize=10)
ax[1].set(title='B  Concentration of supported targets',xlabel='Taxa ranked by supported target count',ylabel='Cumulative supported targets (%)',xlim=(0,n),ylim=(0,103));ax[1].legend(loc='lower right',frameon=False,fontsize=9);ax[1].grid(alpha=.2)
fig.suptitle('Duplication control availability is limited and uneven',fontsize=15,y=.98)
fig.text(.06,.045,'Same family + conservative shared architecture + sequence distance within a factor of 1.5.\nStrict background set; modeled terminal duplicate targets only. Guide alternatives are not independent replicates.',fontsize=9,color='#444444')
fig.subplots_adjust(left=.07,right=.98,bottom=.23,top=.85,wspace=.3)
for ext in ['png','pdf','svg']:
 path=Path(str(a.prefix)+'.'+ext);fig.savefig(path,dpi=180);outputs.append(path)
plt.close(fig)
for path,h in bindings.items():assert sha(path)==h
result=dict(status='complete_background_support_figure',source_bindings=bindings,artifacts={str(p):sha(p) for p in outputs},coverage_rows=len(c),concentration_rows=len(q),scope='Descriptive availability, no effect estimates or confidence intervals. Panel A uses all observed modeled target records as denominators, including identical models; focal means a background includes the duplicate taxon. Panel B ranks taxa by support using alignment E-value policy; zero-support taxa remain. Absent target taxa are outside these denominators. Strict set requires native ortholog assignment and unreported parent duplication in both guide alternatives; this does not establish speciation. Other qualification sets and distance ranges remain in the full source table.')
a.receipt.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(coverage_rows=len(c),concentration_rows=len(q))))
