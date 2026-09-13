#!/usr/bin/env python3
"""Plot per-protein experimental agreement, deposition multiplicity and coverage exclusions."""
import argparse,json
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['comparisons','summary','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable figure output')
    cr=checked_receipt(a.comparisons);sr=checked_receipt(a.summary)
    if sr['comparison_receipt_sha256']!=sha(a.comparisons/'receipt.json'):raise ValueError('Summary source mismatch')
    rows=read_table(a.comparisons/'comparisons.tsv');excluded=read_table(a.comparisons/'exclusions.tsv');proteins=read_table(a.summary/'protein_summary.tsv')
    baseline=[r for r in proteins if r['predicted_plddt_cutoff']=='0'];baseline.sort(key=lambda r:float(r['median_CA_RMSD_angstrom']))
    order={r['predicted_model_id']:i for i,r in enumerate(baseline)}
    if len(order)!=len(baseline):raise ValueError('Duplicate protein summaries')
    table=[]
    for cut in [0,70,90]:
        included=[r for r in rows if int(r['predicted_plddt_cutoff'])==cut];ex=[r for r in excluded if int(r['predicted_plddt_cutoff'])==cut]
        medians=[float(r['median_CA_RMSD_angstrom']) for r in proteins if int(r['predicted_plddt_cutoff'])==cut]
        for protein in proteins:
            if int(protein['predicted_plddt_cutoff'])!=cut:continue
            values=[float(r['ca_superposition_rmsd_angstrom']) for r in included if r['predicted_model_id']==protein['predicted_model_id']]
            if len(values)!=int(protein['chain_model_comparisons']) or not np.isclose(np.median(values),float(protein['median_CA_RMSD_angstrom']),rtol=1e-12,atol=1e-12):raise ValueError('Protein summary mismatch')
        table.append({'predicted_plddt_cutoff':cut,'accepted_comparisons':len(included),'excluded_comparisons':len(ex),'accepted_proteins':len(medians),'comparison_weighted_median_RMSD':float(np.median([float(r['ca_superposition_rmsd_angstrom']) for r in included])),'median_of_protein_median_RMSD':float(np.median(medians))})
    a.output.mkdir(parents=True);write_table(a.output/'plotted_threshold_summary.tsv',table)
    fig,(left,right)=plt.subplots(1,2,figsize=(12,10),gridspec_kw={'width_ratios':[1.4,1]})
    colors={0:'#335c81',70:'#b55c29',90:'#258777'}
    for cut,offset in [(0,-.20),(70,0),(90,.20)]:
        subset=[r for r in proteins if int(r['predicted_plddt_cutoff'])==cut]
        if any(r['predicted_model_id'] not in order for r in subset):raise ValueError('Filtered protein absent from baseline')
        left.scatter([float(r['median_CA_RMSD_angstrom']) for r in subset],[order[r['predicted_model_id']]+offset for r in subset],color=colors[cut],s=24,label=f'Predicted pLDDT ≥{cut}',zorder=3)
    left.set_yticks(range(len(baseline)),[r['predicted_model_id'].removeprefix('AF-').removesuffix('-F1')+'  (n='+r['chain_model_comparisons']+')' for r in baseline],fontsize=8)
    left.set(xscale='log',xlabel='Within-protein median Cα RMSD (Å; log scale)',title='Each protein shown once per threshold')
    left.grid(axis='x',alpha=.2);left.legend(loc='lower right',fontsize=8)
    ticks=np.arange(3);width=.33
    right.bar(ticks-width/2,[r['comparison_weighted_median_RMSD'] for r in table],width,color='#808c98',label='Median across comparisons')
    right.bar(ticks+width/2,[r['median_of_protein_median_RMSD'] for r in table],width,color='#335c81',label='Median across protein medians')
    right.set_xticks(ticks,[f"≥{r['predicted_plddt_cutoff']}\n{r['accepted_proteins']} proteins\n{r['accepted_comparisons']} accepted\n{r['excluded_comparisons']} excluded" for r in table],fontsize=8)
    right.set(ylabel='Cα RMSD (Å)',xlabel='Predicted focal pLDDT threshold',title='Deposition multiplicity changes the summary')
    right.legend(fontsize=8);right.grid(axis='y',alpha=.2)
    for ax in [left,right]:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Initial prediction–experiment agreement: partial reference coverage',fontsize=14)
    fig.text(.5,.04,'n = accepted chain/model comparisons at pLDDT ≥0; points are dependent descriptive summaries.\nThresholds change residue subsets and protein coverage. One reference protein fails coverage at every threshold.\nExperimental quality, biological context and training independence remain unresolved.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.10,1,.96])
    for ext in ['svg','png','pdf']:fig.savefig(a.output/('experimental_agreement.'+ext),dpi=180)
    plt.close(fig)
    receipt={'status':'complete_experimental_agreement_figure','comparison_receipt_sha256':sha(a.comparisons/'receipt.json'),'summary_receipt_sha256':sha(a.summary/'receipt.json'),'script_sha256':sha(Path(__file__)),'plotted_proteins':len(baseline),'recomputed_threshold_summary':table,'interpretation':'Every plotted per-protein count and median independently recomputed from accepted comparison rows. Descriptive partial snapshot; no inferential independence or accuracy claim. Baseline cohort includes 33 of 34 reference proteins; no point is drawn for absent threshold coverage.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
