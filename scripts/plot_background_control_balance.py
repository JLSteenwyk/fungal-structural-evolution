#!/usr/bin/env python3
"""Show within-match balance separately from selection shifts in the original target pool."""
import argparse,csv,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_ortholog_pair_guide_comparison import sha

FEATURES=['sequence_distance','log_positive_sequence_distance','mean_log_length','log_length_asymmetry','mean_plddt','minimum_plddt','mean_lowconf_fraction','maximum_lowconf_fraction']
LABELS=['Sequence distance','Log positive sequence distance','Mean log protein length','Log length asymmetry','Mean pLDDT','Minimum pLDDT','Mean low-confidence fraction','Maximum low-confidence fraction']
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--proof',type=Path,required=True);p.add_argument('--prefix',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    rp=a.source/'receipt.json';r=json.loads(rp.read_text());proof=json.loads(a.proof.read_text());assert proof['status']=='passed_full_background_control_balance_readback' and proof['producer_receipt_sha256']==sha(rp)
    paths=[rp,a.proof]
    for name,h in r['artifacts'].items():assert sha(a.source/name)==h;paths.append(a.source/name)
    bindings={str(p):sha(p) for p in paths};data={};coverage={}
    for row in csv.DictReader(open(a.source/'covariate_balance.tsv'),delimiter='\t'):
        if row['scenario_id']=='S45' and row['policy']=='alignment_evalue':data[row['guide'],row['feature']]=row
    for row in csv.DictReader(open(a.source/'selection_coverage.tsv'),delimiter='\t'):
        if row['scenario_id']=='S45' and row['policy']=='alignment_evalue':coverage[row['guide']]=row
    assert set(data)=={(g,f) for g in ['profile','mafft'] for f in FEATURES}
    assert all(r['smd_status']=='estimable' and r['selection_shift_status']=='estimable' for r in data.values())
    points=[]
    for guide in ['profile','mafft']:
        for feature in FEATURES:
            row=data[guide,feature]
            for panel,field in [('matched_balance','standardized_mean_difference'),('target_selection_shift','selection_shift_in_baseline_sd')]:points.append(dict(panel=panel,guide=guide,feature=feature,value=row[field],pairs=row['pairs'],baseline_targets=row['baseline_targets']))
    a.prefix.parent.mkdir(parents=True,exist_ok=True);table=Path(str(a.prefix)+'.tsv')
    with table.open('w') as f:
        w=csv.DictWriter(f,list(points[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(points)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(12,6),sharey=True)
    for ax,panel,title in zip(axes,['matched_balance','target_selection_shift'],['A  Selected targets versus controls','B  Selected targets versus original pool']):
        ax.axvline(0,color='#999999',linewidth=1)
        for i,(guide,color,marker) in enumerate([('profile','#126782','o'),('mafft','#b04c15','s')]):
            vals=[float(next(r['value'] for r in points if r['panel']==panel and r['guide']==guide and r['feature']==feature)) for feature in FEATURES]
            ax.scatter(vals,[j+(i-.5)*.16 for j in range(len(FEATURES))],color=color,marker=marker,s=38,label=guide,zorder=3)
        ax.set_title(title,fontsize=12);ax.set_xlim(-1.05,1.05);ax.grid(axis='x',alpha=.2)
    axes[0].set_yticks(range(len(FEATURES)),LABELS);axes[0].invert_yaxis();axes[0].set_xlabel('Mean difference / pooled matched SD')
    axes[1].set_xlabel('Target mean shift / original-target SD');axes[1].legend(frameon=False,loc='lower right')
    c=coverage['profile'];fig.suptitle('Close matches represent a selected target subset',fontsize=15,y=.97)
    fig.text(.04,.08,f"S45 · strict background set · factor-1.5 sequence distance · moderate metadata tolerances · alignment-E-value policy\n{int(c['matched_target_records']):,} matched targets across {int(c['matched_taxa'])} taxa per guide; original pools: {int(coverage['profile']['all_target_records']):,} profile / {int(coverage['mafft']['all_target_records']):,} mafft targets.\nDescriptive standardized differences, not effects or confidence intervals. Positive-log distance excludes {int(c['zero_sequence_distance_targets']):,} zero-distance pairs.",fontsize=9,color='#444444')
    fig.subplots_adjust(left=.26,right=.98,top=.85,bottom=.25,wspace=.18);outputs=[table]
    for ext in ['png','pdf','svg']:
        path=Path(str(a.prefix)+'.'+ext);fig.savefig(path,dpi=180);outputs.append(path)
    plt.close(fig)
    for path,h in bindings.items():assert sha(path)==h
    result=dict(status='complete_background_control_balance_figure',source_bindings=bindings,artifacts={str(p):sha(p) for p in outputs},points=len(points),scenario='S45',policy='alignment_evalue',scope='Illustrative scenario from the full verified sensitivity grid, not chosen as a confirmatory primary test. Identical horizontal scales, but matched pooled SD and original-target SD denominators differ explicitly. All points descriptive; guides are sensitivity alternatives, not independent replicates. Figures show metadata balance and target selection, not structural effects or phylogenetic correction.')
    a.receipt.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
