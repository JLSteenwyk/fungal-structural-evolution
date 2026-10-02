#!/usr/bin/env python3
"""Exercise full 315-cell placement and two-page rendering on declared fixtures."""
import argparse
import copy
import csv
import json
from pathlib import Path

from summarize_coalescent_tree_comparisons import panel_data,render
from readback_coalescent_comparison_figure import inspect
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    prior=Path('results/software-checks/coalescent-tree-comparison-cases-20261002-v1')
    completion=json.loads((prior/'receipt.json').read_text())
    assert completion['status']=='passed_full_coalescent_reference_comparison_software_grid'
    bindings={}
    for path in [str(Path(__file__)), 'scripts/summarize_coalescent_tree_comparisons.py',
                 'scripts/readback_coalescent_comparison_figure.py',str(prior/'receipt.json'),
                 str(prior/'full-grid/comparisons.tsv')]:bind(bindings,path)
    for path,digest in completion['source_hashes'].items():bind(bindings,path,digest)
    # Rename the old generic fixture identities only in memory, keeping all original output intact.
    names=['coal:'+a+':'+s for a in ['profile','mafft'] for s in ['uncontracted','sh_alrt10','sh_alrt80']]
    names+=['concat:'+a+':'+k for a in ['profile_profile','profile_mafft','mafft_profile','mafft_mafft'] for k in ['ml','consensus']]
    mapping={('coal:' if i<6 else 'concat:')+str(i).zfill(2):name for i,name in enumerate(names)}
    with (prior/'full-grid/comparisons.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    for row in rows:
        for field in ['view_a','view_b']:row[field]=mapping[row[field]]
    panels,placements=panel_data(completion,rows)
    paths=render(panels,args.output,software_fixture=True)
    with (args.output/'figure_pair_cells.tsv').open('x') as f:
        writer=csv.DictWriter(f,list(placements[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(placements)
    with (args.output/'figure_pair_cells.tsv').open() as f:placements=list(csv.DictReader(f,delimiter='\t'))
    summary=inspect(completion,rows,panels,placements,args.output/'coalescent_reference_sensitivity.pdf')
    (args.output/'panel_data.json').write_text(json.dumps(panels,indent=2,allow_nan=False)+'\n')
    mutations=[]
    def altered(name,change):
        ps,cs=copy.deepcopy(panels),copy.deepcopy(placements);change(ps,cs);mutations.append((name,ps,cs))
    altered('cross_cell_value',lambda p,c:p[0]['cross_rf'][0].__setitem__(0,.123456))
    altered('internal_cell_value',lambda p,c:p[0]['internal_rf_lower_triangle'][1].__setitem__(0,.123456))
    altered('unused_upper_cell',lambda p,c:p[0]['internal_rf_lower_triangle'][0].__setitem__(1,0.))
    altered('cohort_taxon_count',lambda p,c:p[0].__setitem__('taxa',526))
    altered('cohort_order',lambda p,c:p.reverse())
    altered('reference_order',lambda p,c:p[0]['reference_views'].reverse())
    altered('placement_row',lambda p,c:c[0].__setitem__('row','5'))
    altered('placement_column',lambda p,c:c[0].__setitem__('column','7'))
    altered('placement_panel',lambda p,c:c[0].__setitem__('panel','coalescent_internal'))
    altered('source_rf',lambda p,c:c[0].__setitem__('rf_distance','999'))
    altered('missing_pair',lambda p,c:c.pop())
    altered('duplicate_pair',lambda p,c:c.__setitem__(1,c[0].copy()))
    rejected=0
    for name,ps,cs in mutations:
        try:inspect(completion,rows,ps,cs)
        except (AssertionError,ValueError,KeyError,TypeError):rejected+=1
        else:raise AssertionError('Altered figure accepted: '+name)
    verify(bindings)
    result=dict(status='passed_full_315_cell_coalescent_figure_software_contract',**summary,
        altered_exports_rejected=rejected,source_hashes=bindings,
        artifacts={p.name:sha(p) for p in paths+[args.output/'panel_data.json',args.output/'figure_pair_cells.tsv']},
        scientific_eligibility=False,scope='All 315 synthetic pairs from the already declared five-cohort/70-view software grid, renamed only in memory to production view identities. Full independently reconstructed cell placements, RF and normalized RF, unused cells, all 315 PDF numeric labels and two-page rendering checked. False cell/placement/scope/missing/duplicate exports rejected. No production inference, source/journal closure or biological pilot evidence.')
    with (args.output/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,bool))}),flush=True)


if __name__=='__main__':main()
