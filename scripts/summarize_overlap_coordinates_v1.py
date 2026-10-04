#!/usr/bin/env python3
"""Describe the complete coordinate benchmark without inferential or accuracy claims."""
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


METRICS=['ca_superposition_rmsd_angstrom','all_distance_rms_change_angstrom',
         'local_distance_rms_change_angstrom','sequence_local_distance_rms_change_angstrom']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['producer','reader','reader-transport','output','figure-prefix','receipt']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();assert not a.output.exists() and not a.receipt.exists()
    producer=json.loads(a.producer.read_text());reader=json.loads(a.reader.read_text());transport=json.loads(a.reader_transport.read_text())
    assert reader['status']=='passed_full_independent_matched_predictor_coordinate_readback'
    assert reader['producer_receipt_sha256']==sha(a.producer)
    assert transport['validation_sha256']==sha(a.reader) and transport['original_tool_terminal_exit_code']==0
    pins=dict(reader['source_hashes'])
    for q in [a.producer,a.reader,a.reader_transport,Path(__file__)]:bind(pins,q)
    table=Path(producer['table']);assert sha(table)==producer['table_sha256'];verify(pins)
    records=list(csv.DictReader(table.open(),delimiter='\t'));assert len(records)==7716
    groups={}
    for cutoff in ['0','70','90']:
        for pae in ['unfiltered','5','10','15']:
            group=[r for r in records if r['plddt_cutoff']==cutoff and r['pae_cutoff']==pae]
            assert len(group)==643;groups[cutoff,pae]=group
    summaries=[];correlations=[]
    for (cutoff,pae),group in groups.items():
        good=[r for r in group if r['status']=='coordinates_compared']
        values={metric:np.array([float(r[metric]) for r in good]) for metric in METRICS}
        summary=dict(plddt_cutoff=int(cutoff),pae_cutoff=pae,declared_model_pairs=643,compared_model_pairs=len(good),
            insufficient_common_coverage=sum(r['status']=='insufficient_common_coverage' for r in group),
            coordinate_validation_rejected=sum(r['status']=='coordinate_validation_rejected' for r in group),
            retained_residues=sum(int(r['retained_residues']) for r in good),
            state_mismatches=sum(int(r['state_mismatches']) for r in good))
        for metric,v in values.items():
            for name,q in [('p025',.025),('median',.5),('p975',.975)]:
                summary[metric+'_'+name]=float(np.quantile(v,q,method='linear')) if len(v) else ''
        summaries.append(summary)
        state=np.array([int(r['state_mismatches'])/int(r['retained_residues']) for r in good])
        partner=np.array([int(r['partner_changes'])/int(r['retained_residues']) for r in good])
        for field,axis in [('state_mismatch_fraction',state),('partner_change_fraction',partner)]:
            for metric,v in values.items():
                rho=float(spearmanr(axis,v).statistic) if len(v)>1 and np.ptp(axis)>0 and np.ptp(v)>0 else ''
                correlations.append(dict(plddt_cutoff=int(cutoff),pae_cutoff=pae,model_pairs=len(good),
                    x=field,y=metric,descriptive_spearman_rho=rho,scope='Dependent selected predictor pairs; no p-value or evolutionary attribution'))
    a.output.mkdir()
    for name,data in [('threshold_geometry_summary.tsv',summaries),('descriptive_geometry_correlations.tsv',correlations)]:
        with (a.output/name).open('x') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
    primary=[r for r in groups['70','10'] if r['status']=='coordinates_compared']
    assert len(primary)==545
    state=np.array([int(r['state_mismatches'])/int(r['retained_residues']) for r in primary])
    contact=np.array([float(r['local_distance_rms_change_angstrom']) for r in primary])
    seq=np.array([float(r['sequence_local_distance_rms_change_angstrom']) for r in primary])
    ca=np.array([float(r['ca_superposition_rmsd_angstrom']) for r in primary])
    allpair=np.array([float(r['all_distance_rms_change_angstrom']) for r in primary])
    counts=[len([r for r in g if r['status']=='coordinates_compared']) for g in groups.values()]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(11,8.5),layout='constrained')
    ax=axes[0,0];im=ax.scatter(contact,state*100,c=ca,s=17,cmap='viridis',alpha=.75,edgecolor='none')
    ax.set(xscale='log',xlabel='Spatial-contact distance RMS change (Å)',ylabel='3Di state disagreement (%)',title='A  States versus local geometry')
    fig.colorbar(im,ax=ax,label='Cα superposition RMSD (Å)')
    ax=axes[0,1];ax.scatter(seq,allpair,s=17,color='#2073a1',alpha=.65,edgecolor='none')
    lo=min(seq.min(),allpair.min());hi=max(seq.max(),allpair.max());ax.plot([lo,hi],[lo,hi],color='0.5',lw=1,ls='--')
    ax.set(xscale='log',yscale='log',xlabel='Sequence-local distance RMS change (Å)',ylabel='All-pair distance RMS change (Å)',title='B  Local and whole-protein geometry')
    ax=axes[1,0]
    data=[seq,contact,allpair,ca];colors=['#277da1','#43aa8b','#f8961e','#f94144']
    for i,(v,color) in enumerate(zip(data,colors),1):
        sorted_values=np.sort(v);ax.plot(sorted_values,np.arange(1,len(v)+1)/len(v),color=color,lw=2,
            label=['Sequence-local distances','Spatial-contact distances','All-pair distances','Cα superposition'][i-1])
    ax.set(xscale='log',xlabel='RMS difference (Å)',ylabel='Fraction of 545 compared pairs',title='C  Conditional empirical distributions')
    ax.legend(fontsize=8,loc='lower right')
    ax=axes[1,1];grid=np.array(counts).reshape(3,4)
    ax.imshow(grid,vmin=0,vmax=643,cmap='Blues',aspect='auto')
    for i in range(3):
        for j in range(4):ax.text(j,i,f'{grid[i,j]} / 643',ha='center',va='center',color='white' if grid[i,j]>400 else 'black')
    ax.set(xticks=range(4),xticklabels=['Unfiltered','≤5 Å','≤10 Å','≤15 Å'],yticks=range(3),yticklabels=['0','70','90'],
        xlabel='Maximum feature-context PAE',ylabel='Minimum feature-context pLDDT',title='D  Coverage at every original mask')
    fig.suptitle('Matched AlphaFold–ESMFold coordinate benchmark\nPanels A–C: pLDDT ≥70 / PAE ≤10 Å; same complete sequences',fontsize=14)
    fig.supxlabel('Selected overlap control · 21 fungi · 78 linked markers · 643 model pairs\nWhole-protein predictor differences; no experimental accuracy or evolutionary-rate interpretation',fontsize=9)
    artifacts=[]
    for suffix in ['.png','.pdf']:
        q=Path(str(a.figure_prefix)+suffix);assert not q.exists();q.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(q,dpi=180);artifacts.append(q)
    plt.close(fig)
    for q in [*a.output.iterdir(),*artifacts]:bind(pins,q)
    verify(pins)
    primary_summary=next(s for s in summaries if s['plddt_cutoff']==70 and s['pae_cutoff']=='10')
    result=dict(status='complete_descriptive_full_overlap_coordinate_summary_and_figure',checked_utc=datetime.now(timezone.utc).isoformat(),
        model_pairs=643,confidence_settings=12,source_rows=7716,threshold_summary_rows=12,descriptive_correlation_rows=96,
        primary_mask_summary=primary_summary,figure_paths=[str(q) for q in artifacts],source_hashes=pins,
        scientific_eligibility=False,all_eight_aims_incomplete=True,
        scope='Full independently replayed coordinate benchmark summarized at every original threshold, with explicit coverage. '
              'Linear empirical quantiles and descriptive rank correlations only; neither calibrated confidence intervals nor '
              'p-values. Model pairs, marker/taxon links and alternative thresholds are not independent evolutionary replications. '
              'Panels A-C condition on pLDDT70/PAE10; panel D retains every cohort count. Whole-protein RMSD includes domain '
              'orientation effects. Source model pairs link78markers/21fungi; ready phylogenetic branch controls cover71markers, '
              'a different cohort. Full501fungi+25outgroup project and original branch/path calibration remain unfinished.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
