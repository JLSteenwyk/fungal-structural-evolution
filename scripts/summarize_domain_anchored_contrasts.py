#!/usr/bin/env python3
"""Summarize paired A-reference minus B-reference residuals and plot every case."""
import csv,json
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from screen_duplication_domain_alignment_coverage import sha

METRICS=['anchor_rmsd','outside_domain_anchored_rmsd','outside_independent_fit_rmsd']
KEY=['domain_triad','whole_triad','mask','order_ab','order_ar','order_br','mapping_definition']
CASE=['family','gene_a','gene_b','pfam_accession']


def main():
    base=Path('results/structural_comparisons');sources={}
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name
        assert sha(p)==r['artifacts'][name];sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return p
    p=checked(base/'domain-anchored-displacement-20260927-v1','pair_displacements.tsv')
    d=pd.read_csv(p,sep='\t');assert len(d)==2496 and d.status.eq('computed').all()
    assert not d.duplicated(KEY+['pair']).any()
    ar=d[d.pair.eq('ar')].set_index(KEY);br=d[d.pair.eq('br')].set_index(KEY)
    assert ar.index.equals(br.index) and len(ar)==832
    contrasts=(ar[METRICS]-br[METRICS]).reset_index()
    # Independent full CSV reconstruction ensures pairing precedes subtraction.
    scalar=defaultdict(dict)
    with p.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):scalar[tuple(row[k] for k in KEY)][row['pair']]=row
    for row in contrasts.to_dict('records'):
        group=scalar[tuple(str(row[k]) for k in KEY)];assert set(group)=={'ab','ar','br'}
        for metric in METRICS:assert abs(row[metric]-(float(group['ar'][metric])-float(group['br'][metric])))<1e-12
    links=pd.read_csv(checked(base/'whole-domain-case-dossiers-20260927-v1','domain_reference_links.tsv'),sep='\t')
    relation=links[CASE+['triad_key']].drop_duplicates().rename(columns={'triad_key':'domain_triad'})
    linked=relation.merge(contrasts,on='domain_triad',how='left',validate='one_to_many')
    assert len(linked)==832 and linked[METRICS].notna().all().all()
    coverage=pd.read_csv(checked(base/'domain-anchored-coverage-20260927-v1','case_screens.tsv'),sep='\t')
    aggregated=[]
    for key,group in linked.groupby(CASE,sort=True):
        for metric in METRICS:
            values=group[metric].to_numpy()
            aggregated.append(dict(zip(CASE,key),measurement=metric,alternatives=len(values),minimum=float(min(values)),maximum=float(max(values))))
    ranges=pd.DataFrame(aggregated);assert len(ranges)==39
    results=coverage.merge(ranges,on=CASE,how='left',validate='many_to_many');assert len(results)==234
    out=base/'domain-anchored-contrast-summary-20260927-v1';out.mkdir(exist_ok=False);artifacts={}
    for name,data in [('paired_contrasts.tsv',linked),('case_contrast_ranges.tsv',results)]:
        path=out/name;data.to_csv(path,sep='\t',index=False)
        back=pd.read_csv(path,sep='\t');pd.testing.assert_frame_equal(back,data.reset_index(drop=True),check_dtype=False,rtol=1e-12,atol=1e-12)
        artifacts[name]=sha(path)
    plotdata=results[results.screen.eq('n50_c70')];order=list(coverage[coverage.screen.eq('n50_c70')].family)
    assert len(order)==13 and len(set(order))==13
    titles=['Inside domain, domain fit','Outside domain, domain fit','Outside domain, independent fit']
    fig,axes=plt.subplots(1,3,figsize=(14,8),sharey=True)
    labels=[]
    for family in order:
        row=plotdata[plotdata.family.eq(family)].iloc[0]
        labels.append(row.species_name.strip()+'\n'+family+' / '+row.pfam_accession)
    for ax,metric,title in zip(axes,METRICS,titles):
        sub=plotdata[plotdata.measurement.eq(metric)].set_index('family')
        ax.axvline(0,color='black',lw=.8);ax.axvspan(-.1,.1,color='#e8e8e8',zorder=0)
        for y,family in enumerate(order):
            row=sub.loc[family];okay=int(row.all_anchors_and_outside_pass)==1
            color='#007a87' if okay else '#999999'
            ax.plot([row.minimum,row.maximum],[y,y],color=color,lw=2.4)
            ax.plot([row.minimum,row.maximum],[y,y],linestyle='',marker='|',color=color,markersize=8)
        ax.set_title(title,fontsize=11);ax.set_xlabel('RMSD(A, reference) − RMSD(B, reference) (Å)',fontsize=9)
        ax.grid(axis='x',alpha=.2);ax.spines[['top','right']].set_visible(False)
    axes[0].set_yticks(range(13),labels,fontsize=8);axes[0].invert_yaxis()
    fig.suptitle('Domain-anchored contrasts: all 13 inspection cases',fontsize=15,y=.98)
    fig.text(.37,.025,'Ranges cover all alternatives, not confidence intervals.\nTeal: all anchor and outside sets pass n50/c70. Gray: incomplete coverage. Shading: ±0.1 Å descriptive margin.',fontsize=9,ha='left')
    fig.subplots_adjust(left=.26,right=.99,top=.92,bottom=.14,wspace=.25)
    for ext in ['png','pdf','svg']:
        path=out/('domain_anchored_contrasts.'+ext);fig.savefig(path,dpi=180);artifacts[path.name]=sha(path)
    plt.close(fig)
    for path,digest in sources.items():assert sha(path)==digest
    receipt=dict(status='complete_domain_anchored_paired_contrast_summary_pending_visual_review',source_hashes=sources,script_sha256=sha(__file__),paired_configurations=len(contrasts),case_measurement_ranges=len(ranges),case_screen_measurement_rows=len(results),artifacts=artifacts,scope='A-reference minus B-reference residuals paired by the identical anchor/map configuration before range aggregation. All cases and settings retained; n50/c70 used only for figure coverage coloring. Ranges are alternatives, not uncertainty intervals; outside anchored residuals are not independently superposed RMSDs.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
