#!/usr/bin/env python3
"""Plot the complete audited refinement/baseline and gamma-bound comparisons."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    closure_path = Path('metadata/whole_optimization_probability_comparison_completed_20260927.json')
    closure = json.loads(closure_path.read_text())
    receipt_path = Path(closure['completed_receipt_path'])
    assert sha(receipt_path) == closure['completed_receipt_sha256']
    receipt = json.loads(receipt_path.read_text())
    table = receipt_path.parent / 'comparison_summary.tsv'
    assert sha(table) == receipt['artifacts'][table.name]
    rows = list(csv.DictReader(table.open(), delimiter='\t'))
    assert len(rows) == 702
    kinds = ['refined_versus_baseline', 'refined_lower_versus_original_bound']
    families = sorted({r['base_job_id'].split('-')[0] for r in rows})
    assert len(families) == 13
    data = []
    for family in families:
        for kind in kinds:
            subset = [r for r in rows if r['base_job_id'].startswith(family+'-') and r['comparison'] == kind]
            assert len(subset) == (36 if kind == kinds[0] else 18)
            data.append(dict(family=family, comparison=kind,
                node_sites=sum(int(r['columns']) for r in subset),
                state_changes=sum(int(r['map_disagreements']) for r in subset),
                opposing_high_support=sum(int(r['both_at_least_090_disagreements']) for r in subset),
                maximum_total_variation=max(float(r['maximum_total_variation']) for r in subset)))
    for kind in kinds:
        subset = [r for r in data if r['comparison'] == kind]
        expected = receipt['descriptive_aggregates'][kind]
        assert sum(r['state_changes'] for r in subset) == expected['map_disagreements']
        assert sum(r['node_sites'] for r in subset) == expected['sites']
        assert max(r['maximum_total_variation'] for r in subset) == expected['maximum_total_variation']
    out = Path('results/ancestral/whole-refinement-sensitivity-figure-20260927-v2')
    out.mkdir(parents=True, exist_ok=False)
    with (out/'family_summary.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,list(data[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(data)
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig, axes = plt.subplots(1,2,figsize=(11,6),sharey=True,gridspec_kw={'width_ratios':[1,1.4]})
    y=np.arange(len(families));colors=['#0072B2','#D55E00']
    for index,kind in enumerate(kinds):
        subset=[r for r in data if r['comparison']==kind]
        offset=-.15 if index==0 else .15
        label=['Refined vs original fit','Lower vs original gamma bound'][index]
        axes[0].barh(y+offset,[r['state_changes'] for r in subset],height=.28,color=colors[index])
        values=[r['maximum_total_variation'] for r in subset]
        assert min(values)>0
        axes[1].scatter(values,y+offset,color=colors[index],s=28,label=label)
    axes[0].set_yticks(y, families);axes[0].invert_yaxis()
    axes[0].set_xlabel('Changes in most probable amino acid')
    axes[0].set_xticks(range(0,15,2));axes[0].set_xlim(0,14)
    axes[1].set_xscale('log');axes[1].set_xlabel('Maximum probability-distribution change\n(total variation; log scale)')
    for ax in axes:
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
    fig.legend(*axes[1].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.91),ncol=2,frameon=False,fontsize=9)
    fig.suptitle('Ancestral amino-acid estimates: refinement sensitivity',fontsize=14)
    fig.text(.5,.04,'13 families; 156 refined fits. All columns and three identifiable ancestors retained.\nCounts include dependent model/bound comparisons; they are not independent substitutions.\nNo opposing amino-acid calls have ≥90% support in both fits. Indel uncertainty is separate.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.16,1,.87))
    for ext in ['png','pdf','svg']:
        fig.savefig(out/('whole_refinement_sensitivity.'+ext),dpi=180)
    plt.close(fig)
    result=dict(status='complete_audited_refinement_sensitivity_figure',families=13,
        input_receipt=str(receipt_path),input_receipt_sha256=sha(receipt_path),
        script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='Descriptive optimization sensitivity from all702 audited comparison rows. No uncertainty interval, independence or biological significance claim. Alternate starts and indels are not represented.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'output':str(out),'families':13,'state_changes':sum(r['state_changes'] for r in data)}))

if __name__=='__main__':
    main()
