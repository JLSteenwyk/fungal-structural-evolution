#!/usr/bin/env python3
"""Plot descriptive pooled and within-taxon structural availability."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from catalog_whole_proteome_structures import sha


def main():
    root=Path('results/structures/annotation-stratified-coverage-20260923-v1')
    source=root/'taxon_annotation_coverage.tsv';receipt=json.loads((root/'receipt.json').read_text())
    if sha(source)!=receipt['artifacts'][source.name]:raise ValueError('Changed coverage table')
    t=pd.read_csv(source,sep='\t');w=t.pivot(index='taxon_id',columns='annotation_status',values='modeled_fraction')
    if len(w)!=receipt['taxa'] or w.isna().any().any():raise ValueError('Incomplete taxon strata')
    recomputed=t.modeled/t.proteins
    if not np.allclose(recomputed,t.modeled_fraction,rtol=0,atol=1e-14):raise ValueError('Coverage arithmetic differs')
    pooled=t.groupby('annotation_status')[['proteins','modeled']].sum()
    for key,row in pooled.iterrows():
        if any(int(row[field])!=receipt['pooled'][key][field] for field in ['proteins','modeled']):raise ValueError('Pooled counts differ')
    delta=100*(w.with_pfam_hit-w.without_pfam_hit)
    fig,ax=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    colors=['#28698D','#B16B35'];labels=['With Pfam hit','Without Pfam hit'];keys=['with_pfam_hit','without_pfam_hit']
    values=[100*pooled.loc[k,'modeled']/pooled.loc[k,'proteins'] for k in keys]
    ax[0].bar(labels,values,color=colors,width=0.6)
    for i,v in enumerate(values):ax[0].text(i,v+0.7,f'{v:.2f}%',ha='center')
    ax[0].set(ylim=(0,30),ylabel='Proteins with a catalog model (%)',title='A  Pooled protein records')
    bins=np.linspace(-30,30,25);counts,edges=np.histogram(delta,bins=bins)
    if counts.sum()!=len(w):raise ValueError('Histogram excludes taxa')
    ax[1].stairs(counts,edges,fill=True,color='#477F72');ax[1].axvline(0,color='black',lw=0.9)
    ax[1].set(xlabel='With-hit minus no-hit coverage (percentage points)',ylabel='Number of taxa',title='B  Within-taxon coverage differences')
    ax[1].text(.98,.97,f'Higher: {(delta>0).sum()} taxa\nLower: {(delta<0).sum()} taxa\nEqual: {(delta==0).sum()} taxa',transform=ax[1].transAxes,ha='right',va='top',fontsize=9)
    for a in ax:a.spines[['top','right']].set_visible(False)
    fig.suptitle('Frozen AlphaFold availability by annotation status',fontsize=13)
    out=Path('docs/figures');out.mkdir(exist_ok=True)
    artifacts={}
    for suffix in ['svg','pdf']:
        p=out/('structure_annotation_coverage.'+suffix);fig.savefig(p);artifacts[str(p)]=sha(p)
    fig.savefig('/tmp/structure_annotation_coverage.png',dpi=130);plt.close(fig)
    result={'status':'complete_descriptive_structure_annotation_coverage_figure','taxa':len(w),'taxa_higher_with_hit':int((delta>0).sum()),'taxa_lower_with_hit':int((delta<0).sum()),'taxa_equal':int((delta==0).sum()),'median_within_taxon_difference_percentage_points':float(delta.median()),'range_within_taxon_difference_percentage_points':[float(delta.min()),float(delta.max())],'histogram_counts':counts.tolist(),'histogram_edges':edges.tolist(),'source_sha256':sha(source),'source_receipt_sha256':sha(root/'receipt.json'),'script_sha256':sha(__file__),'artifacts':artifacts,'scope':'Descriptive coverage only; pooled observations weight taxa by protein counts, while histogram counts taxa equally. No confidence-qualified coverage, independent statistical replicates, equivalence, causal or evolutionary inference.'}
    Path('metadata/structure_annotation_coverage_figure_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['histogram_counts','histogram_edges','artifacts']},indent=2))


if __name__=='__main__':main()
