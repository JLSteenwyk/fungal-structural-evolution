#!/usr/bin/env python3
"""Plot descriptive domain-boundary extension counts with explicit denominators."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--boundaries',type=Path,required=True);ap.add_argument('--pfam-table',type=Path,required=True)
    ap.add_argument('--pfam-receipt',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError('Use fresh output')
    source=a.boundaries/'boundary_sensitivity.tsv.gz';rp=a.boundaries/'receipt.json';readback=a.boundaries/'readback.json'
    r=json.loads(rp.read_text());rb=json.loads(readback.read_text());pr=json.loads(a.pfam_receipt.read_text())
    if (rb['status']!='passed_all_boundary_pair_rows_and_threshold_counts' or rb['producer_receipt_sha256']!=sha(rp)
            or sha(source)!=r['artifacts'][source.name]
            or pr['status']!='complete_pfam_boundary_sensitivity_with_full_sql_aggregation_readback'
            or sha(a.pfam_table)!=pr['artifacts']['pfam_boundary_sensitivity.tsv']
            or pr['source_hashes'][str(rp)]!=sha(rp) or pr['source_hashes'][str(readback)]!=sha(readback)):
        raise ValueError('Missing or mismatched source validation')
    pins={str(p):sha(p) for p in [source,rp,readback,a.pfam_table,a.pfam_receipt]}
    labels=['0 (identical)','1–4','5–9','10–19','20–49','50 or more'];counts=Counter();total=large=0
    with gzip.open(source,'rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            extra=int(row['added_residues']);length=int(row['alignment_residues'])
            if extra<0 or length<=0:raise ValueError('Invalid boundary length')
            index=0 if extra==0 else 1 if extra<5 else 2 if extra<10 else 3 if extra<20 else 4 if extra<50 else 5
            counts[index]+=1;total+=1;large+=5*extra>=length
    if (total!=r['candidate_model_hit_pairs'] or counts[0]!=r['identical_boundary_pairs']
            or sum(counts[i] for i in [3,4,5])!=r['pairs_with_at_least_10_added_residues']
            or large!=r['pairs_with_at_least_20_percent_added_residues']):raise ValueError('Histogram/threshold totals differ')
    with a.pfam_table.open() as f:pfams=list(csv.DictReader(f,delimiter='\t'))
    if len(pfams)!=pr['pfam_accessions'] or sum(int(x['candidate_pairs']) for x in pfams)!=total or sum(int(x['added_ge20pct']) for x in pfams)!=large:raise ValueError('Pfam aggregation differs')
    ranked=sorted(pfams,key=lambda x:(-int(x['added_ge20pct']),x['pfam_accession']))[:10]
    data=[dict(panel='A',label=label,numerator=counts[i],denominator=total,percent=100*counts[i]/total) for i,label in enumerate(labels)]
    data.extend(dict(panel='B',label=x['pfam_accession'],numerator=int(x['added_ge20pct']),denominator=int(x['candidate_pairs']),percent=100*int(x['added_ge20pct'])/int(x['candidate_pairs'])) for x in ranked)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(13,6.4),gridspec_kw={'width_ratios':[1,1.45]})
    fig.subplots_adjust(left=.115,right=.975,bottom=.23,top=.79,wspace=.47)
    fig.suptitle('Domain boundaries: alignment span versus HMM envelope',x=.07,ha='left',y=.97,fontsize=17)
    fig.text(.07,.9,f'{total:,} candidate model/hit pairs · {r["source_models"]:,} models · {pr["pfam_accessions"]:,} Pfam accessions',fontsize=11,color='#444444')
    checks=0
    for ax in axes:
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',color='#dddddd',linewidth=.6);ax.set_axisbelow(True)
    selected=data[:6];values=[x['percent'] for x in selected]
    bars=axes[0].barh(range(6),values,color='#346c91',height=.62)
    axes[0].set_yticks(range(6),labels);axes[0].invert_yaxis();axes[0].set_xlim(0,max(values)*1.36)
    axes[0].set_xlabel('Percentage of all candidate pairs');axes[0].set_ylabel('Envelope residues added')
    axes[0].set_title('A  Extension length',loc='left',fontsize=12,pad=15)
    for b,d in zip(bars,selected):
        assert abs(b.get_width()-d['percent'])<1e-12;checks+=1
        axes[0].text(b.get_width()+.6,b.get_y()+b.get_height()/2,f'{d["percent"]:.1f}%\n({d["numerator"]:,})',va='center',fontsize=9)
    selected=data[6:];values=[x['numerator'] for x in selected]
    bars=axes[1].barh(range(len(selected)),values,color='#bd7035',height=.62)
    axes[1].set_yticks(range(len(selected)),[x['label'] for x in selected]);axes[1].invert_yaxis();axes[1].set_xlim(0,max(values)*1.6)
    axes[1].set_xlabel('Model/hit pairs with ≥20% added length')
    axes[1].set_title('B  Top 10 Pfams by count of ≥20% extensions',loc='left',fontsize=12,pad=15)
    for b,d in zip(bars,selected):
        assert b.get_width()==d['numerator'];checks+=1
        axes[1].text(b.get_width()+12,b.get_y()+b.get_height()/2,f'{d["numerator"]:,}/{d["denominator"]:,} ({d["percent"]:.1f}%)',va='center',fontsize=9)
    fig.text(.07,.12,'Panel A uses disjoint length bins. Panel B labels show affected / all candidate pairs within each Pfam.',fontsize=10)
    fig.text(.07,.07,'Counts are model/hit observations, not independent species or evolutionary changes. Thresholds are descriptive.',fontsize=10,color='#555555')
    a.output.mkdir(parents=True)
    with (a.output/'plot_data.tsv').open('w') as f:
        w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    for ext in ['svg','pdf','png']:fig.savefig(a.output/('domain_boundary_sensitivity.'+ext),dpi=170,facecolor='white')
    plt.close(fig)
    svg=a.output/'domain_boundary_sensitivity.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    for path,h in pins.items():
        if sha(path)!=h:raise ValueError('Source changed during plotting')
    receipt=dict(status='complete_descriptive_domain_boundary_figure_pending_visual_review',candidate_pairs=total,pfam_accessions=len(pfams),bars_checked=checks,source_hashes=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in a.output.iterdir()},scope='Full source-table aggregation and all 16 bar lengths checked. Displayed Pfams selected by absolute affected-pair count; proportions descriptive. Not enrichment, structural-confidence, cluster-stability or evolutionary-effect inference.')
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
