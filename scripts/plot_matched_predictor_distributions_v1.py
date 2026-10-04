#!/usr/bin/env python3
"""Plot validated descriptive predictor ranges, preserving their conditional scope."""
from collections import Counter
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import verify


def main():
    rp=Path('metadata/matched_predictor_distribution_readback_20261004_v1.json')
    tp=Path('metadata/matched_predictor_distribution_readback_transport_20261004_v1.json')
    r,t=[json.loads(p.read_text()) for p in [rp,tp]]
    assert r['status']=='passed_full_independent_paired_predictor_distribution_summary_readback'
    assert t['validation_sha256']==sha(rp) and t['original_tool_terminal_exit_code']==0
    assert t['whole_wrapper_initial_and_terminal_payloads_matched']
    verify(r['source_hashes']);verify(t['source_hashes'])
    table=Path('results/phylogeny/matched-predictor-distribution-summary-20261004-v1/paired_branch_distributions.tsv')
    counts=Counter();n=0
    with table.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            assert row['calibrated_effect_interval']==row['scientific_eligibility']=='False'
            counts[row['mode']+':'+row['model']+':'+row['empirical_direction']]+=1;n+=1
    assert n==13170 and dict(counts)==r['empirical_direction_counts']
    fig,axes=plt.subplots(1,2,figsize=(10,5.8),sharey=True)
    model_names=['AF+G4','AF+F+G4','LLM+G4'];models=['af','af_empirical','llm']
    labels=['Range includes zero','Range entirely below zero','Range entirely above zero']
    directions=['includes_zero','negative','positive'];colors=['#b8c2cc','#315b83','#c17635']
    for ax,mode,title in zip(axes,['iid_sites','circular_blocks10'],['Independent-site draws','Circular blocks of 10 retained columns']):
        bottom=np.zeros(3)
        for direction,label,color in zip(directions,labels,colors):
            totals=[sum(counts[mode+':'+model+':'+d] for d in directions) for model in models]
            assert totals==[2195]*3
            value=np.array([counts[mode+':'+m+':'+direction]/2195 for m in models])
            ax.bar(np.arange(3),value,bottom=bottom,color=color,label=label,width=.65)
            bottom+=value
        assert np.allclose(bottom,1)
        ax.set_xticks(np.arange(3),model_names);ax.set_title(title,fontsize=11)
        ax.set_ylim(0,1);ax.set_yticks(np.arange(0,1.01,.2),[f'{int(x)}%' for x in np.arange(0,101,20)])
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Fraction of branch/input combinations')
    fig.suptitle('Matched AlphaFold–ESMFold predictor sensitivity',fontsize=14,y=.96)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.15),ncol=3,frameon=False,fontsize=9)
    fig.text(.5,.045,'ESMFold minus AlphaFold branch length; empirical 2.5th–97.5th percentile ranges from 200 draws.\n'
             '71 markers, 21 fungal taxa; 2,195 branch/input combinations per model and mode.\n'
             'Fixed predictions and topology; overlapping topologies/branches are dependent. These are not calibrated discoveries.',
             ha='center',fontsize=9)
    fig.subplots_adjust(left=.09,right=.98,top=.84,bottom=.30,wspace=.12)
    paths=[Path('docs/figures/matched_predictor_paired_ranges_20261004.'+ext) for ext in ['png','pdf']]
    for path in paths:
        assert not path.exists();fig.savefig(path,dpi=180)
    plt.close(fig)
    proof=dict(status='plotted_verified_full_paired_predictor_empirical_ranges',summary_rows=13170,
        branch_input_combinations_per_model_mode=2195,markers=71,fungal_taxa=21,
        counts=dict(counts),source_hashes={str(p):sha(p) for p in [rp,tp,table,Path(__file__)]},
        artifacts={str(p):sha(p) for p in paths},scientific_eligibility=False,
        calibrated_discoveries=False,gpu=False,new_native_fits=0)
    out=Path('metadata/matched_predictor_distribution_figure_20261004_v1.json')
    with out.open('x') as f:json.dump(proof,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in proof.items() if k not in ['source_hashes','counts']},indent=2))


if __name__=='__main__':main()
