#!/usr/bin/env python3
"""Plot checked lineage coverage with sequence and model gaps separated."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True)
    p.add_argument('--readback',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    receipt=json.loads(a.receipt.read_text());readback=json.loads(a.readback.read_text())
    if (receipt['status']!='complete_frozen_marker_gap_characterization'
            or readback['status']!='passed_gap_sequence_and_lineage_summary_readback'
            or readback['source_receipt_sha256']!=sha(a.receipt)
            or receipt['artifacts']['lineage_marker_availability.tsv']!=sha(a.summary)):
        raise ValueError('Unbound or unverified coverage table')
    with a.summary.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    rows.sort(key=lambda r:(r['study_role']!='ingroup',r['manifest_lineage']))
    fields=['taxa','marker_slots','recovered_records','unrecovered_sequences','model_linked_records','uncovered_recovered_records']
    for row in rows:
        for field in fields:row[field]=int(row[field])
        if (row['model_linked_records']+row['uncovered_recovered_records']!=row['recovered_records']
                or row['recovered_records']+row['unrecovered_sequences']!=row['marker_slots']
                or row['marker_slots']!=125*row['taxa']):raise ValueError('Invalid coverage partition')
    assert sum(r['taxa'] for r in rows)==526
    assert sum(r['model_linked_records'] for r in rows)==56171
    a.output.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.hashsalt':'current-marker-coverage-20260922'})
    fig,(left,right)=plt.subplots(1,2,figsize=(13,12),sharey=True,gridspec_kw={'width_ratios':[4.6,1.55]})
    y=np.arange(len(rows),dtype=float);y[[r['study_role']=='outgroup' for r in rows]]+=.6
    offset=np.zeros(len(rows));segments=[]
    categories=[('model_linked_records','#26817e','Sequence and model available'),
                ('uncovered_recovered_records','#e4be74','Sequence available; model missing'),
                ('unrecovered_sequences','#dedede','Marker sequence not recovered')]
    for field,color,label in categories:
        values=np.array([r[field]/r['marker_slots'] for r in rows])
        bars=left.barh(y,values,left=offset,height=.72,color=color,label=label,edgecolor='white',linewidth=.25)
        for i,patch in enumerate(bars):
            if not np.isclose(patch.get_width(),values[i],rtol=0,atol=1e-12) or not np.isclose(patch.get_x(),offset[i],rtol=0,atol=1e-12):
                raise ValueError('Plotted coverage interval differs')
            segments.append(dict(study_role=rows[i]['study_role'],manifest_lineage=rows[i]['manifest_lineage'],measure=field,
                                 count=rows[i][field],denominator=rows[i]['marker_slots'],left_fraction=float(offset[i]),width_fraction=float(values[i])))
        offset+=values
    if not np.allclose(offset,1,rtol=0,atol=1e-12):raise ValueError('Incomplete plotted bars')
    labels=[r['manifest_lineage']+(' [outgroup]' if r['study_role']=='outgroup' else '')+f"  (n={r['taxa']})" for r in rows]
    left.set_yticks(y,labels);left.invert_yaxis();left.set_xlim(0,1);left.xaxis.set_major_formatter(PercentFormatter(1))
    left.set_xlabel('Fraction of all 125 marker slots per sampled entry')
    left.set_title('Sequence recovery and model availability',loc='left',fontsize=12,pad=15)
    conditional=np.array([r['model_linked_records']/r['recovered_records'] for r in rows])
    dots=right.scatter(conditional,y,color='#26817e',s=28,zorder=3)
    if not np.allclose(np.asarray(dots.get_offsets())[:,0],conditional,rtol=0,atol=1e-12):raise ValueError('Plotted conditional coverage differs')
    right.set_xlim(.85,1.01);right.set_xticks([.85,.90,.95,1]);right.xaxis.set_major_formatter(PercentFormatter(1,decimals=0))
    right.set_xlabel('Models among recovered\nmarker sequences')
    right.set_title('Conditional coverage',loc='left',fontsize=12,pad=15)
    for ax in [left,right]:
        ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
        ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True)
    fig.suptitle('Marker coverage across the full sampling design',x=.04,ha='left',fontsize=17,y=.985)
    fig.text(.04,.953,'526 sampled entries · 125 markers · September 22 frozen model catalogs',fontsize=11,color='#444444')
    handles,legend_labels=left.get_legend_handles_labels()
    fig.legend(handles,legend_labels,loc='lower left',bbox_to_anchor=(.035,.066),frameon=False,ncol=1,fontsize=10)
    fig.text(.04,.043,'93.9% of recovered marker records have models; 85.4% of all marker slots have both sequence and model.',fontsize=10)
    fig.text(.04,.019,'Availability precedes confidence qualification. Lineage bins are not equivalent taxonomic ranks; n counts sampled entries.\nUnrecovered sequences do not establish gene absence. The right panel uses a restricted 85–101% axis.',fontsize=8.5,color='#444444')
    fig.subplots_adjust(left=.33,right=.98,top=.90,bottom=.19,wspace=.18)
    for suffix in ['svg','pdf','png']:
        metadata={'Date':None} if suffix=='svg' else {'CreationDate':None,'ModDate':None} if suffix=='pdf' else None
        fig.savefig(a.output/('current_marker_coverage.'+suffix),dpi=170,metadata=metadata)
    plt.close(fig)
    svg=a.output/'current_marker_coverage.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    with (a.output/'plotted_segments.tsv').open('w') as f:
        w=csv.DictWriter(f,list(segments[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(segments)
    result=dict(status='complete_checked_marker_coverage_figure',lineage_rows=len(rows),plotted_segments=len(segments),
                conditional_points=len(rows),input_hashes={str(path):sha(path) for path in [a.summary,a.receipt,a.readback]},
                script_sha256=sha(__file__),artifacts={path.name:sha(path) for path in a.output.iterdir()},
                scope='Descriptive frozen-catalog availability for all sampled entries, with distinct sequence/model gaps and declared denominators. No confidence-qualified coverage, model accuracy or evolutionary effect claims.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
