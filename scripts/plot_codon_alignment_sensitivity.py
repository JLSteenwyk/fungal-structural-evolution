#!/usr/bin/env python3
"""Export full-case alignment sensitivity distributions and retention plot."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root = Path('results/cds/codon-divergence-sensitivity-summary-20260927-v1')
    proof_path = Path('metadata/codon_divergence_sensitivity_summary_readback_20260927.json')
    proof = json.loads(proof_path.read_text())
    assert proof['status'] == 'passed_full_codon_divergence_sensitivity_summary_readback'
    assert proof['source_receipt_sha256'] == sha(root/'receipt.json')
    receipt = json.loads((root/'receipt.json').read_text())
    for name,digest in receipt['artifacts'].items():assert sha(root/name)==digest
    tree_root=Path('results/cds/codon-alignment-tree-comparison-20260927-v1')
    assert receipt['topology_receipt_sha256']==sha(tree_root/'receipt.json')
    tr=json.loads((tree_root/'receipt.json').read_text())
    assert tr['artifacts']['cases.tsv']==sha(tree_root/'cases.tsv')
    frame=pd.read_csv(root/'matched_cases.tsv',sep='\t',keep_default_na=False).set_index('case_id')
    trees=pd.read_csv(tree_root/'cases.tsv',sep='\t',keep_default_na=False).set_index('case_id')
    data=frame[['topology_changed','historical_review_flags','selection_eligibility','global_omega_ratio_disposition','global_omega_log2_local_over_original','tree_ds_equal_alternative_log2_local_over_original','tree_dn_equal_alternative_log2_local_over_original']].copy()
    data['median_pair_retention']=trees.loc[data.index,'alignment_median_pair_retention'].astype(float)
    assert len(data)==1625 and data.index.is_unique
    prefix=Path('docs/figures/codon_alignment_sensitivity_20260927')
    for suffix in ['.png','.pdf','.svg','.tsv','.receipt.json']:assert not prefix.with_suffix(suffix).exists()
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,(left,right)=plt.subplots(1,2,figsize=(12.8,6.3))
    fig.subplots_adjust(left=.08,right=.98,bottom=.30,top=.78,wspace=.28)
    fig.suptitle('Codon divergence depends on alignment choices',fontsize=17,y=.96)
    fig.text(.08,.88,'1,625 matched gene groups • full distributions and all cases retained',fontsize=11)
    columns=[('tree_ds_equal_alternative_log2_local_over_original','Tree-total dS','#2166ac'),('tree_dn_equal_alternative_log2_local_over_original','Tree-total dN','#d95f02'),('global_omega_log2_local_over_original','Fitted global omega','#4d9221')]
    for col,label,color in columns:
        values=pd.to_numeric(data[col],errors='coerce').dropna().sort_values().to_numpy()
        assert np.isfinite(values).all()
        left.step(values,np.arange(1,len(values)+1)/len(values),where='post',color=color,label=f'{label} (n={len(values):,})')
    left.axvline(0,color='.4',lw=.8,ls='--')
    left.set(xlabel='log₂(local / original)',ylabel='Cumulative fraction of groups',ylim=(0,1.02))
    left.set_title('A  Divergence sensitivity',loc='left',pad=12)
    left.legend(frameon=False,fontsize=9,loc='lower right')
    for changed,color,label in [(False,'#2166ac','Unchanged topology'),(True,'#d95f02','Changed topology')]:
        subset=data[data.topology_changed==changed]
        right.scatter(subset.median_pair_retention,subset.tree_dn_equal_alternative_log2_local_over_original.astype(float),s=16,c=color,alpha=.5,edgecolors='none',label=f'{label} (n={len(subset):,})',rasterized=True)
    right.axhline(0,color='.4',lw=.8,ls='--')
    right.set(xlabel='Median retained original residue-pair fraction',ylabel='Tree-total dN: log₂(local / original)')
    right.set_xlim(min(data.median_pair_retention)-.025,1.025)
    right.set_title('B  Alignment retention and divergence',loc='left',pad=12)
    right.legend(frameon=False,fontsize=9,loc='lower left')
    for ax in (left,right):ax.spines[['top','right']].set_visible(False)
    fig.text(.08,.18,'0 = unchanged; +1 = doubled; −1 = halved. Tree-total dS/dN use the declared opportunity normalization.',fontsize=10)
    fig.text(.08,.125,'One zero/zero global-omega case has no log ratio; it remains in the data export and panel B.',fontsize=10)
    fig.text(.08,.07,'Descriptive sensitivity, not selection evidence. Coverage, alignment, topology and fitted parameters can all contribute.',fontsize=10)
    prefix.parent.mkdir(parents=True,exist_ok=True)
    data.to_csv(prefix.with_suffix('.tsv'),sep='\t')
    for suffix in ['.png','.pdf','.svg']:
        path=prefix.with_suffix(suffix)
        fig.savefig(path,dpi=180)
        if suffix=='.svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)
    report=dict(status='complete_codon_alignment_sensitivity_figure_pending_readback',cases=len(data),summary_receipt_sha256=sha(root/'receipt.json'),summary_readback_sha256=sha(proof_path),topology_receipt_sha256=sha(tree_root/'receipt.json'),script_sha256=sha(__file__),artifacts={prefix.with_suffix(s).name:sha(prefix.with_suffix(s)) for s in ['.png','.pdf','.svg','.tsv']},scope='All matched groups shown; ECDFs use defined positive-value log ratios without clipping or pseudocounts. Retention is median across taxon pairs per case; no causal regression or significance claim. Historical flags retained in export.')
    prefix.with_suffix('.receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
