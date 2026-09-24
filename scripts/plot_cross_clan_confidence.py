#!/usr/bin/env python3
"""Plot descriptive cross-clan candidate retention across confidence thresholds."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--summary',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();rp=args.summary/'receipt.json';table=args.summary/'candidate_confidence_summary.tsv'
    pins={str(p):sha(p) for p in [rp,table,Path(__file__)]};receipt=json.loads(rp.read_text())
    if receipt['status']!='complete_descriptive_cross_clan_confidence_summary' or receipt['artifacts'][table.name]!=sha(table):raise ValueError('Incomplete or changed summary')
    with table.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    thresholds=[0,.5,.8,.9];data=[];seen=set()
    for r in rows:
        key=tuple(r[k] for k in ['boundary','representative','pfam_left','pfam_right','confidence_fraction'])
        if key in seen:raise ValueError('Duplicate summary row')
        seen.add(key)
    for t in thresholds:
        total_screen=0
        for boundary in ['alignment','envelope']:
            group=[r for r in rows if r['boundary']==boundary and float(r['confidence_fraction'])==t]
            n=len(group);k=sum(int(r['screen_interval_pairs'])>0 for r in group)
            if n!={'alignment':32,'envelope':36}[boundary]:raise ValueError('Candidate scope differs')
            data.append(dict(confidence_fraction=t,boundary=boundary,candidates=n,screened=k));total_screen+=k
        if total_screen!=receipt['screen_entries_by_threshold'][str(t)]:raise ValueError('Screen count differs')
    plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(10,6.5));fig.subplots_adjust(left=.10,right=.98,bottom=.30,top=.78)
    fig.suptitle('Cross-clan domain candidates: confidence sensitivity',fontsize=16,y=.96)
    fig.text(.10,.885,'Entries with ≥1 interval pair meeting TM-score ≥0.5 and coverage ≥0.8',fontsize=11)
    fig.text(.10,.84,'Criteria must hold for both alignment directions and both length normalizations.',fontsize=10)
    for boundary,offset,color in [('alignment',-.19,'#346c91'),('envelope',.19,'#b97239')]:
        subset=[r for r in data if r['boundary']==boundary]
        bars=ax.bar([i+offset for i in range(4)],[r['screened'] for r in subset],width=.34,color=color,label=boundary.capitalize()+' boundaries')
        for bar,row in zip(bars,subset):
            if bar.get_height()!=row['screened']:raise ValueError('Bar height differs')
            ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+.5,f"{row['screened']}/{row['candidates']}",ha='center',fontsize=10)
    ax.set_xticks(range(4),['None','50%','80%','90%']);ax.set_ylim(0,40)
    ax.set_ylabel('Candidate entries passing the screen')
    ax.set_xlabel('Minimum fraction with C-alpha pLDDT ≥70')
    ax.spines[['top','right']].set_visible(False);ax.set_axisbelow(True);ax.grid(axis='y',alpha=.2)
    ax.legend(frameon=False,loc='upper right')
    fig.text(.10,.18,'Confidence applies to both whole domains and both matched-residue alignments.',fontsize=10)
    fig.text(.10,.125,'Labels: screened / all entries in that boundary view; empty strata remain in the denominator.',fontsize=9)
    fig.text(.10,.07,'Views and models overlap. Descriptive prioritization only; no homology or evolutionary inference.',fontsize=9)
    args.output.mkdir(parents=True,exist_ok=False)
    with (args.output/'plot_data.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    for ext in ['svg','pdf','png']:fig.savefig(args.output/('cross_clan_confidence.'+ext),dpi=160,facecolor='white')
    plt.close(fig)
    svg=args.output/'cross_clan_confidence.svg';svg.write_text('\n'.join(l.rstrip() for l in svg.read_text().splitlines())+'\n')
    for p,h in pins.items():
        if sha(p)!=h:raise ValueError('Input changed')
    result=dict(status='complete_cross_clan_confidence_figure_pending_visual_review',source_hashes=pins,
                bars_checked=len(data),artifacts={p.name:sha(p) for p in args.output.iterdir()},
                scope='Descriptive screen counts; 68 overlapping boundary-specific entries, not 68 independent biological discoveries.')
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
