#!/usr/bin/env python3
"""Describe all-family context disagreement and local-coverage limitations."""
import csv,hashlib,json
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path('results/ancestral/ancestral-context-coverage-20260927-v1')
    receipt=json.loads((root/'receipt.json').read_text());table=root/'coverage_summary.tsv'
    assert sha(table)==receipt['artifacts'][table.name]
    with table.open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
    families=sorted({r['family'] for r in rows});assert len(families)==13
    data=[]
    for family in families:
        rs=[r for r in rows if r['family']==family]
        n=sum(int(r['comparisons']) for r in rs);changes=sum(int(r['map_disagreements']) for r in rs)
        high={c:sum(int(r['opposing_at_least_090']) for r in rs if r['coverage_class']==c) for c in receipt['coverage_classes']}
        assert high['neither_has_known_descendants']==0
        data.append(dict(family=family,comparisons=n,map_disagreements=changes,percent_disagreement=100*changes/n,
            opposing_with_known_descendants=high['both_have_known_descendants'],opposing_without_domain_descendants=high['whole_only_has_known_descendants']))
    assert sum(r['comparisons'] for r in data)==269154 and sum(r['map_disagreements'] for r in data)==10320
    assert sum(r['opposing_with_known_descendants'] for r in data)==56 and sum(r['opposing_without_domain_descendants'] for r in data)==48
    output=Path('results/ancestral/ancestral-context-sensitivity-figure-20260927-v1');output.mkdir(exist_ok=False)
    with (output/'family_summary.tsv').open('w') as h:
        w=csv.DictWriter(h,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,(a,b)=plt.subplots(1,2,figsize=(12,7),sharey=True,gridspec_kw={'width_ratios':[1.2,1]})
    y=np.arange(13);vals=[r['percent_disagreement'] for r in data]
    a.barh(y,vals,color='#666666',height=.65)
    for i,r in enumerate(data):a.text(vals[i]+.2,i,f"{r['map_disagreements']:,}/{r['comparisons']:,}",va='center',fontsize=8)
    a.set_xlim(0,max(vals)+6);a.set_xlabel('Different most-probable amino acids (%)')
    a.set_title('A  All matched-coordinate comparisons',loc='left',fontsize=11)
    a.set_yticks(y,families);a.invert_yaxis()
    v=np.array([r['opposing_with_known_descendants'] for r in data]);u=np.array([r['opposing_without_domain_descendants'] for r in data])
    b.barh(y,v,color='#0072B2',height=.65,label='Known descendants in both contexts')
    b.barh(y,u,left=v,color='#D55E00',height=.65,label='No known descendants at domain position')
    for i in range(13):
        if v[i]:b.text(v[i]/2,i,str(v[i]),ha='center',va='center',color='white')
        if u[i]:b.text(v[i]+u[i]/2,i,str(u[i]),ha='center',va='center',color='white')
    b.set_xlim(0,75);b.set_xlabel('Opposing calls supported ≥0.9 in both fits')
    b.set_title('B  High-confidence conflicts by local coverage',loc='left',fontsize=11)
    for ax in [a,b]:ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
    fig.suptitle('Whole-protein versus domain ancestral estimates',fontsize=15)
    fig.legend(*b.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.925),ncol=2,frameon=False,fontsize=9)
    fig.text(.5,.04,'13 families; 624 fit pairs; 269,154 dependent ancestor/site comparisons. Labels in A are changed/compared counts.\nThe 104 high-confidence conflicts represent three exact coordinate sets across five ancestral nodes.\nObserved descendant coverage does not establish ancestral residue presence; indel uncertainty remains unresolved.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.15,1,.89))
    for ext in ['png','pdf','svg']:fig.savefig(output/('ancestral_context_sensitivity.'+ext),dpi=180)
    plt.close(fig)
    result=dict(status='complete_all_family_context_sensitivity_figure',input_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),families=13,comparisons=269154,artifacts={p.name:sha(p) for p in output.iterdir()},scope='Descriptive dependent comparisons only; no independence, evolutionary event counts, posterior mixture weighting or statistical significance claim.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
