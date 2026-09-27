#!/usr/bin/env python3
"""Plot every eligible candidate's original and sequence-locked contrast ranges."""
import argparse,csv,json,hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text())
    for f,h in p['pins'].items():assert sha(f)==h
    proof=json.loads(Path(p['readback']).read_text());assert proof['status']=='passed_full_sequence_locked_reference_core_readback'
    root=Path(p['source']);receipt=json.loads((root/'receipt.json').read_text());assert sha(root/'receipt.json')==proof['producer_receipt_sha256']
    table=root/'candidate_control_summary.tsv';assert sha(table)==receipt['artifacts'][table.name]
    with table.open() as f:source=list(csv.DictReader(f,delimiter='\t'))
    assert len(source)==48 and all(x['control_complete']=='1' for x in source)
    source.sort(key=lambda x:(int(x['control_retains_original_direction']),x['family'],x['gene_a'],x['gene_b'],x['pfam_accession']))
    rows=[]
    for i,x in enumerate(source,1):
        rows.append(dict(display_row=i,**{k:x[k] for k in ['family','gene_a','gene_b','pfam_accession','pfam_name','candidate_class','study_role','control_direction_margin_0_1','control_retains_original_direction']},original_min=float(x['all_alternatives_min_rmsd_ar_minus_br']),original_max=float(x['all_alternatives_max_rmsd_ar_minus_br']),control_min=float(x['control_contrast_min']),control_max=float(x['control_contrast_max'])))
    fig,axes=plt.subplots(1,2,figsize=(12,14.8),sharex=True,sharey=True)
    palette={'a_farther_from_reference':'#187d75','b_farther_from_reference':'#187d75','direction_flip_beyond_margin':'#bf3b3b','touches_or_enters_margin_band':'#a96700'}
    labels=[f"{r['display_row']:02d}  {r['family']}  {r['pfam_name']}" for r in rows]
    artists=[]
    for panel,ax in enumerate(axes):
        prefix='original' if panel==0 else 'control'
        ax.axvspan(-.1,.1,color='#ededed',zorder=0);ax.axvline(0,color='#666666',lw=.6,zorder=1)
        for y,row in enumerate(rows):
            color='#5b6570' if panel==0 else palette[row['control_direction_margin_0_1']]
            bounds=[row[prefix+'_min'],row[prefix+'_max']]
            line,=ax.plot(bounds,[y,y],color=color,lw=1.8,marker='|',markersize=5)
            artists.append((line,bounds,[y,y]))
        ax.set_title('Original residue correspondence' if panel==0 else 'Sequence-locked duplicate correspondence',fontsize=12,pad=12)
        ax.set_xlabel('RMSD(A, reference) − RMSD(B, reference) [Å]',fontsize=10)
        ax.grid(axis='x',alpha=.15);ax.spines[['top','right']].set_visible(False);ax.tick_params(axis='y',length=0)
    axes[0].set_yticks(range(len(rows)),labels,fontsize=7.5);axes[0].set_ylim(len(rows)-.3,-.7)
    axes[0].set_xlim(min(r[k] for r in rows for k in ['original_min','control_min'])-.15,max(r[k] for r in rows for k in ['original_max','control_max'])+.15)
    fig.suptitle('Residue correspondence changes four candidate classifications',fontsize=16,y=.988)
    fig.text(.5,.965,'All 48 identical-domain candidates • 77 reference triads • 4,928 constrained fits',ha='center',fontsize=11)
    legend=[Line2D([0],[0],color=c,lw=2,label=l) for c,l in [('#187d75','Original direction retained (44)'),('#bf3b3b','Both signs beyond margin (3)'),('#a96700','Enters margin band (1)')]]
    fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.56,.025),ncol=3,frameon=False,fontsize=9)
    fig.text(.54,.013,'Ranges across alternative fits; shaded band = ±0.1 Å (descriptive). Lines are not confidence intervals.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.30,right=.985,top=.938,bottom=.08,wspace=.12)
    # Match every drawn range to the underlying plot table before export.
    for line,x,y in artists:assert list(line.get_xdata())==x and list(line.get_ydata())==y
    base=Path(p['output_prefix']);base.parent.mkdir(parents=True,exist_ok=True)
    for suffix in ['.png','.pdf','.svg','.tsv']:
        assert not base.with_suffix(suffix).exists()
    for suffix in ['.png','.pdf','.svg']:fig.savefig(base.with_suffix(suffix),dpi=180)
    with base.with_suffix('.tsv').open('w') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    plt.close(fig)
    result=dict(status='complete_sequence_locked_reference_control_figure',plan_sha256=sha(a.plan),source_receipt_sha256=sha(root/'receipt.json'),readback_sha256=sha(p['readback']),candidate_rows=len(rows),range_endpoints_checked=4*len(rows),script_sha256=sha(__file__),artifacts={str(base.with_suffix(s)):sha(base.with_suffix(s)) for s in ['.png','.pdf','.svg','.tsv']},scope='All 48 eligible candidates, both original and controlled extrema, shared axis scale. Descriptive alternative ranges, not uncertainty intervals or biological effect validation.')
    Path(p['receipt']).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
