#!/usr/bin/env python3
"""Plot the complete original grid and all largest-family retry contrasts."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    pins={}
    def audited_table(root, name):
        root=Path(root);rp=root/'receipt.json';r=json.loads(rp.read_text())
        pins[str(rp)]=sha(rp)
        assert sha(root/name)==r['artifacts'][name]
        with (root/name).open() as h:return list(csv.DictReader(h,delimiter='\t'))
    original=audited_table('results/ancestral/historian-capacity-audit-20260927-final-v1','paired_sequences.tsv')
    retries=audited_table('results/ancestral/historian-largest-family-comparison-20260927-v1','comparisons.tsv')
    assert len(original)==5680 and all(int(r['edit_distance'])==0 for r in original)
    assert len(retries)==16
    counts=Counter(r['comparison'] for r in original)
    matrix=np.full((4,4),np.nan)
    for r in retries:
        if r['comparison']=='branch_floor':row=0 if '-mafft-' in r['job_a'] else 1
        else:
            assert r['comparison']=='alignment'
            row=2 if '-floor1e-09-' in r['job_a'] else 3
        level=int(r['level']);assert np.isnan(matrix[row,level]);matrix[row,level]=int(r['edit_distance'])
        assert (r['assumed_root']=='True')==(level==3)
    assert np.isfinite(matrix).all()
    out=Path('results/ancestral/historian-sensitivity-figure-20260927-v1');out.mkdir(exist_ok=False)
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,(a,b)=plt.subplots(1,2,figsize=(12,5.5),gridspec_kw={'width_ratios':[1,1.5]})
    a.set_title('A  Original runs passing output checks',loc='left',fontsize=11)
    labels=['Polytomy resolution','Minimum branch length']
    for y,(label,kind) in enumerate(zip(labels,['star_resolution','branch_floor'])):
        a.scatter(0,y,s=70,color='#0072B2',zorder=3)
        a.text(1,y,f"{counts[kind]:,} comparisons\nall zero edits",va='center',fontsize=10)
    a.set_yticks([0,1],labels);a.set_ylim(1.7,-.7);a.set_xlim(-1,23)
    a.set_xticks([0,5,10,15,20]);a.set_xlabel('Ungapped sequence edit distance')
    a.spines[['top','right']].set_visible(False);a.grid(axis='x',alpha=.2)
    b.set_title('B  Largest family: four successful memory retries',loc='left',fontsize=11)
    im=b.imshow(matrix,cmap='YlOrRd',vmin=0,vmax=25,aspect='auto')
    for row in range(4):
        for col in range(4):b.text(col,row,str(int(matrix[row,col])),ha='center',va='center',color='white' if matrix[row,col]>=18 else 'black')
    b.set_yticks(range(4),['Floor: MAFFT','Floor: FAMSA','Alignment: floor 10⁻⁹','Alignment: floor 10⁻⁷'])
    b.set_xticks(range(4),['Candidate 0','Candidate 1','Candidate 2','Assumed root'])
    b.axvline(2.5,color='black',linestyle='--',linewidth=1)
    b.tick_params(axis='x',labelsize=9)
    fig.colorbar(im,ax=b,label='Ungapped sequence edit distance',fraction=.05,pad=.04)
    fig.suptitle('Ancestral sequence sensitivity under fixed Historian settings',fontsize=14)
    fig.text(.5,.055,'324 configurations checked: 320 original successes + 4 successful retries; original failures retained.\nComparisons share inputs and are not independent biological events or posterior samples.\nPanel A excludes failed original attempts. Panel B shows all 16 contrasts for OG0000972 (622 proteins).',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.20,1,.94),w_pad=2)
    for ext in ['png','pdf','svg']:fig.savefig(out/('historian_sensitivity.'+ext),dpi=180)
    plt.close(fig)
    result=dict(status='complete_full_grid_and_retry_sensitivity_figure',pins=pins,script_sha256=sha(__file__),
        original_comparisons=dict(counts),retry_comparisons=16,
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='All original-grid paired comparisons and all four-alternative retry contrasts shown separately; no posterior reliability or biological event claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
