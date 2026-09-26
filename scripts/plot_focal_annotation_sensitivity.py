#!/usr/bin/env python3
"""Plot focal profile scores beside fixed-correspondence geometry and coverage."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from compare_marker_structures import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable output')
    paths = [Path('metadata/cross_clan_focal_profile_sensitivity_receipt_20260926.json'),
             Path('metadata/cross_clan_focal_profile_sensitivity_readback_20260926.json'),
             Path('metadata/cross_clan_focal_geometry_20260924.json')]
    profile, audit, geometry = [json.loads(p.read_text()) for p in paths]
    if audit['search_receipt_sha256'] != sha(paths[0]) or audit['status'] != 'passed_focal_profile_hit_readback':
        raise ValueError('Profile readback binding failed')
    if geometry['status'] != 'complete_focal_fixed_mapping_geometry':
        raise ValueError('Incomplete geometry')
    hits = [r for run in profile['runs'] if run['label'] == 'max_permissive'
            for r in run['hits'] if r['profile'] == 'PF26973.1']
    order = ['F281847_EPZ34216.1', 'F2606893_XP_031856394.1',
             'F42068_XP_018228094.1', 'F4754_XP_018226490.1', 'F263815_XP_007871995.1']
    names = ['EPZ34216.1 (focal)', 'XP_031856394.1', 'XP_018228094.1 (focal)',
             'XP_018226490.1', 'XP_007871995.1']
    scores = []
    for protein in order:
        rows = [r for r in hits if r['protein'] == protein]
        if not rows or len({r['sequence_score'] for r in rows}) != 1:
            raise ValueError('Missing or inconsistent profile scores')
        scores.append(max(rows, key=lambda r: r['domain_score']))
    geo = [r for r in geometry['results'] if r['boundary'] == 'alignment']
    pairs = [('F281847_EPZ34216.1', 'F2606893_XP_031856394.1', 'Rozella pair', '#0072B2'),
             ('F42068_XP_018228094.1', 'F4754_XP_018226490.1', 'Pneumocystis pair', '#D55E00')]
    grid = {(r['focal'], r['sister'], r['mode'], r['minimum_joint_plddt']): r for r in geo}
    expected = {(f,s,m,c) for f,s,_,_ in pairs for m in ['localpair','globalpair'] for c in [0,70,90]}
    if set(grid) != expected or len(geo) != len(expected):
        raise ValueError('Incomplete geometry plotting grid')
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.7), gridspec_kw={'width_ratios':[1.35,1,1]})
    y = np.arange(len(order))
    axes[0].barh(y-.16, [r['sequence_score'] for r in scores], height=.29, color='#555555', label='Sequence score')
    axes[0].barh(y+.16, [r['domain_score'] for r in scores], height=.29, color='#56B4E9', label='Best domain score')
    axes[0].axvline(26.2, color='#AA3377', linestyle='--', linewidth=1.5, label='GA threshold (both)')
    axes[0].set_yticks(y, names, fontsize=9)
    axes[0].invert_yaxis()
    axes[0].set_xlabel('PF26973 profile score (bits)')
    axes[0].set_title('A  Annotation threshold', loc='left', fontsize=12)
    axes[0].legend(loc='lower right', fontsize=8)
    axes[0].set_xlim(0,48)
    for focal,sister,label,color in pairs:
        for mode,style,marker in [('localpair','-','o'),('globalpair','--','x')]:
            rows=[grid[focal,sister,mode,c] for c in [0,70,90]]
            legend=label+' / '+('local' if mode=='localpair' else 'global')
            axes[1].plot(range(3),[r['ca_rmsd_angstrom'] for r in rows],color=color,linestyle=style,marker=marker,label=legend)
            axes[2].plot(range(3),[r['matched_residues']/r['domain_residues'] for r in rows],color=color,linestyle=style,marker=marker)
    axes[1].set_title('B  Mapped-region geometry',loc='left',fontsize=12)
    axes[1].set_ylabel('Cα RMSD (Å)')
    axes[1].set_ylim(bottom=0)
    axes[2].set_title('C  Residues retained',loc='left',fontsize=12)
    axes[2].set_ylabel('Fraction of focal domain residues')
    axes[2].set_ylim(0,1.02)
    for ax in axes[1:]:
        ax.set_xticks(range(3),['All','≥70','≥90'])
        ax.set_xlabel('Joint pLDDT threshold')
        ax.grid(axis='y',alpha=.2)
    for ax in axes:
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Exploratory annotation sensitivity in a thioredoxin-like family',fontsize=15,y=.97)
    handles,labels=axes[1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.66,.12),ncol=2,fontsize=9,frameon=False)
    fig.text(.02,.055,'A: best domain score anywhere in the protein; positions need not correspond. B–C: PF26973 alignment boundaries, fixed sequence mappings; no PAE mask.',fontsize=9)
    fig.text(.02,.02,'Low RMSD at stricter confidence uses fewer residues. XP_007871995.1 lacks a model in the frozen bridge. These results do not establish domain loss or structural acceleration.',fontsize=9)
    fig.subplots_adjust(left=.145,right=.985,bottom=.29,top=.86,wspace=.46)
    args.output.mkdir(parents=True)
    plot_data={'profile_scores':scores,'geometry':[{k:v for k,v in r.items() if k!='residue_pairs'} for r in geo],
               'missing_models':geometry['missing_from_frozen_bridge']}
    (args.output/'plot_data.json').write_text(json.dumps(plot_data,indent=2)+'\n')
    for ext in ['png','svg','pdf']:
        fig.savefig(args.output/('focal_annotation_sensitivity.'+ext),dpi=180)
    plt.close(fig)
    receipt={'status':'complete_exploratory_focal_annotation_figure',
             'source_sha256':{str(p):sha(p) for p in paths},'script_sha256':sha(Path(__file__)),
             'profile_proteins':len(scores),'geometry_points':len(geo),
             'scope':'Descriptive reuse of completed profile and fixed-mapping geometry results. No new biological tests. Full figure requires visual review.',
             'artifacts':{p.name:sha(p) for p in args.output.iterdir()}}
    (args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(receipt['status'])


if __name__=='__main__':
    main()
