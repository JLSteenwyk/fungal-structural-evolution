#!/usr/bin/env python3
"""Plot descriptive within-group codon divergence without independence claims."""
import argparse,csv,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_busco_gene_copies import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable figure directory')
    r=json.loads((a.input/'receipt.json').read_text());table=a.input/'group_pair_summary.tsv'
    if r['status']!='complete_coverage_group_observed_codon_divergence' or sha(table)!=r['artifacts'][table.name]:raise ValueError('Invalid pair summary')
    with table.open() as f:rows=[x for x in csv.DictReader(f,delimiter='\t') if x['policy']=='exclude_recorded_annotation_gene_label_flags' and x['median_amino_acid_difference_fraction']!='']
    genera=sorted({x['genus_label'] for x in rows});fig,ax=plt.subplots(figsize=(9,5.3));colors=plt.get_cmap('tab20')
    for i,g in enumerate(genera):
        rr=[x for x in rows if x['genus_label']==g]
        ax.scatter([float(x['median_amino_acid_difference_fraction']) for x in rr],[float(x['median_third_position_difference_fraction']) for x in rr],s=13,alpha=.6,color=colors(i),label=g,edgecolors='none')
    ax.set_xlabel('Median amino-acid difference fraction')
    ax.set_ylabel('Median third-position difference fraction')
    ax.set_title('Observed divergence within fungal genus/code groups',loc='left',fontsize=12,pad=28)
    ax.text(0,1.025,'Each point: one marker/group; pairs share ≥100 called codons',transform=ax.transAxes,fontsize=9,color='#555555')
    ax.set_xlim(left=0);ax.set_ylim(bottom=0);ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.16);ax.set_axisbelow(True)
    ax.legend(bbox_to_anchor=(1.02,1),loc='upper left',frameon=False,fontsize=8,markerscale=1.7,title='Genus label',title_fontsize=9)
    fig.text(.10,.015,'Recorded annotation, gene and label flags excluded. Uncorrected differences; shared ancestry is not modeled.',fontsize=8,color='#555555')
    fig.tight_layout(rect=(0,.045,1,1));a.output.mkdir(parents=True)
    for suffix in ['svg','png']:fig.savefig(a.output/('codon_group_divergence.'+suffix),dpi=180)
    plt.close(fig)
    result={'status':'complete_descriptive_codon_divergence_figure','source_receipt_sha256':sha(a.input/'receipt.json'),'source_summary_sha256':sha(table),'points':len(rows),'genus_labels':genera,'script_sha256':sha(Path(__file__)),'interpretation':'Each point is a marker/genus/code group median across eligible overlapping taxon pairs. Points/pairs are not independent; no fit, significance test, dN/dS or saturation estimate is shown. Genus labels are not guaranteed clades.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
