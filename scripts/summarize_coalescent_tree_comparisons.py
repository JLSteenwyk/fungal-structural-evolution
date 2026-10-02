#!/usr/bin/env python3
"""Render all 315 closed candidate comparisons without selecting a preferred tree."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

COHORTS=['full_primary','exclude_sparse_boundary_taxon','exclude_below10_both_alignments',
         'exclude_curated_hybrids','exclude_hybrids_and_uncertain_labels']
COAL=['coal:'+a+':'+s for a in ['profile','mafft'] for s in ['uncontracted','sh_alrt10','sh_alrt80']]
REF=['concat:'+r+':'+k for r in ['profile_profile','profile_mafft','mafft_profile','mafft_mafft'] for k in ['ml','consensus']]
COAL_LABELS=['Profile raw','Profile SH10','Profile SH80','MAFFT raw','MAFFT SH10','MAFFT SH80']
REF_LABELS=[r+'\n'+k for r in ['P/P','P/M','M/P','M/M'] for k in ['ML','cons.']]
COHORT_LABELS=['Primary','Exclude sparse boundary taxon','Exclude low occupancy',
               'Exclude curated hybrids','Exclude hybrids / incomplete labels']


def sources(plan,plan_path):
    completion=json.loads(Path(plan['completion']).read_text())
    assert completion['status']=='complete_verified_full_coalescent_reference_tree_comparisons'
    assert completion['cohorts']==5 and completion['tree_views']==70 and completion['comparison_rows']==315
    assert completion['role_boundary_rows']==70 and completion['exact_process_journals_checked']==2
    bindings=dict(plan['pins']);bind(bindings,plan_path);bind(bindings,plan['completion'])
    bind(bindings,completion['full_hash_archive'],completion['full_hash_archive_sha256'])
    archive=json.loads(Path(completion['full_hash_archive']).read_text());assert len(archive['services'])==2
    for p,d in archive['source_hashes'].items():bind(bindings,p,d)
    producer=json.loads(Path(completion['producer_receipt']).read_text())
    reader=json.loads(Path(completion['independent_readback']).read_text())
    bind(bindings,completion['producer_receipt'],completion['producer_receipt_sha256'])
    bind(bindings,completion['independent_readback'],completion['independent_readback_sha256'])
    assert producer['status']=='complete_coalescent_reference_tree_comparisons_pending_independent_readback'
    assert reader['status']=='passed_coalescent_reference_full_raw_tree_comparison_readback'
    assert reader['producer_receipt_sha256']==sha(completion['producer_receipt'])
    assert producer['comparison_rows']==reader['comparison_rows']==315
    assert producer['cohort_summaries']==reader['cohort_summaries']==completion['cohort_summaries']
    for name,digest in producer['artifacts'].items():
        bind(bindings,Path(completion['producer_receipt']).parent/name,digest)
    verify(bindings)
    path=Path(completion['producer_receipt']).parent/'comparisons.tsv'
    rows=list(csv.DictReader(path.open(),delimiter='\t'))
    assert len(rows)==315
    return completion,rows,bindings


def panel_data(completion,rows):
    index={}
    for r in rows:
        key=r['cohort'],tuple(sorted((r['view_a'],r['view_b'])))
        assert key not in index;index[key]=r
    panels=[];placements=[];used=set()
    summaries={r['cohort']:r for r in completion['cohort_summaries']}
    assert set(summaries)==set(COHORTS)
    for cohort in COHORTS:
        cross=[[None]*8 for _ in range(6)];internal=[[None]*6 for _ in range(6)]
        for i,a in enumerate(COAL):
            for j,b in enumerate(REF):
                key=cohort,tuple(sorted((a,b)));r=index[key];used.add(key)
                assert r['comparison_scope']=='coalescent_vs_concatenated'
                value=float(r['normalized_rf']);assert 0<=value<=1
                expected=int(r['rf_distance'])/(2*(summaries[cohort]['taxa']-3))
                assert abs(value-expected)<1e-12
                cross[i][j]=value
                placements.append(dict(**r,panel='coalescent_reference',row=i,column=j))
            for j in range(i):
                key=cohort,tuple(sorted((a,COAL[j])));r=index[key];used.add(key)
                assert r['comparison_scope']=='coalescent_alignment_support_sensitivity'
                value=float(r['normalized_rf']);assert 0<=value<=1
                assert abs(value-int(r['rf_distance'])/(2*(summaries[cohort]['taxa']-3)))<1e-12
                internal[i][j]=value
                placements.append(dict(**r,panel='coalescent_internal',row=i,column=j))
        panels.append(dict(cohort=cohort,taxa=summaries[cohort]['taxa'],roles=summaries[cohort]['roles'],
            coalescent_views=COAL,reference_views=REF,cross_rf=cross,internal_rf_lower_triangle=internal,
            shared_all_coalescent=summaries[cohort]['shared_all_coalescent'],
            shared_all_references=summaries[cohort]['shared_all_references'],shared_all_views=summaries[cohort]['shared_all_views']))
    assert used==set(index) and len(placements)==315
    return panels,placements


def render(panels,output,software_fixture=False):
    output=Path(output)
    values=[v for p in panels for row in p['cross_rf']+p['internal_rf_lower_triangle'] for v in row if v is not None]
    vmax=max(.05,np.ceil(max(values)/.05)*.05)
    paths=[]
    with PdfPages(output/'coalescent_reference_sensitivity.pdf') as pdf:
        for internal in [False,True]:
            fig,axes=plt.subplots(3,2,figsize=(12.2,13.2),constrained_layout=True)
            cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#eeeeee')
            for k,p in enumerate(panels):
                ax=axes.flat[k]
                matrix=np.array([[np.nan if v is None else v for v in row] for row in
                    (p['internal_rf_lower_triangle'] if internal else p['cross_rf'])])
                im=ax.imshow(matrix,vmin=0,vmax=vmax,cmap=cmap,aspect='auto')
                ax.set_xticks(range(matrix.shape[1]),COAL_LABELS if internal else REF_LABELS,
                              fontsize=7,rotation=45 if internal else 0,ha='right' if internal else 'center')
                ax.set_yticks(range(6),COAL_LABELS,fontsize=8)
                for i,j in np.argwhere(np.isfinite(matrix)):
                    value=matrix[i,j]
                    ax.text(j,i,f'{value:.3f}',ha='center',va='center',fontsize=6.5,
                            color='white' if value/vmax<.45 else '#111111')
                ax.set_title(COHORT_LABELS[k]+'\n'+str(p['taxa'])+' entries; '+str(p['roles']['outgroup'])+' outgroups',fontsize=10)
                ax.set_xlabel('Coalescent candidate' if internal else 'Original/projected concatenated reference',fontsize=8)
                ax.set_ylabel('Coalescent candidate',fontsize=8)
            axes.flat[5].axis('off')
            axes.flat[5].text(0,.98,
                'Normalized Robinson–Foulds distance\n\n'
                'Different splits / [2 × (retained taxa − 3)].\n'
                'Every comparison uses identical retained taxa.\n\n'
                'P = profile; M = MAFFT. Reference labels:\n'
                'alignment / guide, with ML or consensus view.\n'
                'SH10 / SH80 = gene SH-aLRT contraction setting.\n\n'
                +('All 75 unique candidate pairs are shown once;\n'
                  'gray diagonal/upper cells are not additional data.\n\n' if internal else
                  'All 240 coalescent/reference pairs are shown.\n\n')+
                'Point-tree comparisons are descriptive.\n'
                'Support metrics retain separate definitions.\n'
                'No preferred tree, biological root or cause\n'
                'of discordance is assigned.',va='top',fontsize=10,linespacing=1.4)
            fig.colorbar(im,ax=list(axes.flat[:5]),shrink=.65,label='Normalized RF distance',
                         ticks=np.linspace(0,vmax,5),format='%.2f')
            title=('Coalescent alignment/support sensitivity' if internal else
                   'Coalescent versus concatenated species relationships')
            fig.suptitle(('SOFTWARE FIXTURE: ' if software_fixture else '')+title,fontsize=15)
            pdf.savefig(fig)
            name='coalescent_internal_rf.png' if internal else 'coalescent_reference_rf.png'
            path=output/name;fig.savefig(path,dpi=180);paths.append(path);plt.close(fig)
    paths.append(output/'coalescent_reference_sensitivity.pdf')
    return paths


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',required=True,type=Path);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());completion,rows,bindings=sources(plan,args.plan)
    panels,placements=panel_data(completion,rows);out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    with (out/'panel_data.json').open('x') as f:f.write(json.dumps(panels,indent=2,allow_nan=False)+'\n')
    with (out/'figure_pair_cells.tsv').open('x') as f:
        writer=csv.DictWriter(f,list(placements[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(placements)
    paths=render(panels,out)
    verify(bindings)
    result=dict(status='complete_full_coalescent_reference_comparison_figure_pending_cell_readback',
        plan_sha256=sha(args.plan),comparison_completion_sha256=sha(plan['completion']),comparisons=315,
        cross_reference_cells=240,internal_candidate_cells=75,panels=10,pdf_pages=2,tree_views=70,cohorts=5,
        source_hashes=bindings,artifacts={p.name:sha(p) for p in [out/'panel_data.json',out/'figure_pair_cells.tsv']+paths},
        scientific_eligibility=False,scope=plan['scope'])
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,bool))}),flush=True)


if __name__=='__main__':main()
