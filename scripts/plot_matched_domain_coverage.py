#!/usr/bin/env python3
"""Plot verified matched-domain availability across dependent sensitivity settings."""
import csv,json,hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter,MaxNLocator
import pandas as pd

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
root=Path('results/structural_comparisons/selected-domain-coverage-projection-20260927-v2')
proof=Path('metadata/selected_domain_coverage_projection_readback_20260927_v2.json');audit=json.loads(proof.read_text());r=json.loads((root/'receipt.json').read_text())
assert audit['status']=='passed_full_selected_domain_coverage_projection_readback' and audit['source_receipt_sha256']==sha(root/'receipt.json')
source=root/'selected_domain_coverage.tsv';assert sha(source)==r['artifacts'][source.name]
data=pd.read_csv(source,sep='\t');chosen=data[data.scenario_id.isin(['S45','S46']) & data.mask_cohort.eq('both')].copy();assert len(chosen)==192
keys=['guide','policy','scenario_id','boundary','mask_cohort','screen'];assert not chosen.duplicated(keys).any()
chosen['retained_record_fraction']=chosen.usable_selected_records/chosen.selected_records
assert chosen.retained_record_fraction.between(0,1).all()
screens=['n30_c50','n30_c70','n30_c90','n50_c50','n50_c70','n50_c90'];policies=['alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore'];colors={'profile':'#2166ac','mafft':'#b35806'}
fig,axes=plt.subplots(2,2,figsize=(11,7.5),sharex=True)
for col,scenario in enumerate(['S45','S46']):
 subset=chosen[chosen.scenario_id.eq(scenario)]
 for guide in ['profile','mafft']:
  g=subset[subset.guide.eq(guide)]
  x=[screens.index(s)+((policies.index(p)*2+['alignment','envelope'].index(b))-3.5)*.033+(-.018 if guide=='profile' else .018) for s,p,b in zip(g.screen,g.policy,g.boundary)]
  for row,metric in [(0,'retained_record_fraction'),(1,'usable_taxa')]:
   axes[row,col].scatter(x,g[metric],s=24,color=colors[guide],alpha=.65,label=guide+' guide',edgecolors='none')
 axes[0,col].set_title('S45: any eligible background' if scenario=='S45' else 'S46: focal-taxon backgrounds',fontsize=12)
 axes[0,col].set_ylim(0,1);axes[0,col].yaxis.set_major_formatter(PercentFormatter(1))
 axes[1,col].set_ylim(0,160 if scenario=='S45' else 80);axes[1,col].yaxis.set_major_locator(MaxNLocator(integer=True,nbins=5))
 for ax in axes[:,col]:
  ax.set_xticks(range(6),['30 / 50%','30 / 70%','30 / 90%','50 / 50%','50 / 70%','50 / 90%'],rotation=30,ha='right');ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False)
 axes[1,col].set_xlabel('Minimum aligned residues / original-domain coverage')
axes[0,0].set_ylabel('Selected records with a usable domain');axes[1,0].set_ylabel('Taxa with usable matched records')
axes[0,1].legend(frameon=False,loc='lower left')
fig.suptitle('Domain coverage after metadata matching',fontsize=16,y=.99)
fig.text(.06,.025,'Each point: one guide × annotation policy × boundary setting. Both masks must retain the same Pfam domain.\nSensitivity alternatives are dependent; these are coverage counts, not independent replicates or effect estimates.',fontsize=9)
fig.tight_layout(rect=[0,.09,1,.96])
out=Path('docs/figures');stem=out/'matched_domain_coverage_20260927';paths=[]
for suffix in ['png','pdf','svg']:
 p=stem.with_suffix('.'+suffix);assert not p.exists();fig.savefig(p,dpi=180);paths.append(p)
p=stem.with_suffix('.tsv');assert not p.exists();chosen.to_csv(p,sep='\t',index=False);paths.append(p)
# Read every exported identity/count back against the verified complete source.
export=pd.read_csv(p,sep='\t');pd.testing.assert_frame_equal(export.drop(columns='retained_record_fraction').sort_values(keys).reset_index(drop=True),data.merge(chosen[keys],on=keys,validate='one_to_one').sort_values(keys).reset_index(drop=True))
assert all(abs(x-y)<1e-14 for x,y in zip(export.retained_record_fraction,export.usable_selected_records/export.selected_records))
result=dict(status='plotted_verified_domain_coverage_sensitivity',exported_rows=len(export),source_receipt_sha256=sha(root/'receipt.json'),source_readback_sha256=sha(proof),script_sha256=sha(__file__),artifacts={str(x):sha(x) for x in paths},scope='All 192 S45/S46 both-mask guide/policy/boundary/threshold cells plotted. Exact source counts and identities read back; fraction arithmetic checked. Illustrative scenarios previously reported in metadata balance, not selected by structural effects. No uncertainty intervals or independent-replicate claims.')
with stem.with_suffix('.receipt.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
