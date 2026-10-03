#!/usr/bin/env python3
"""Describe the complete coverage geometry only after full independent readback."""
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from audit_selected_taxon_identity_snapshot_v2 import rows,sha,table
from structural_marker_tree_projection import STATUSES,SOURCES


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--producer',type=Path,required=True)
    p.add_argument('--readback',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists() and not a.receipt.exists()
    producer=json.loads(a.producer.read_text());reader=json.loads(a.readback.read_text())
    assert reader['status']=='passed_full_structural_marker_raw_tree_projection_readback'
    assert reader['producer_sha256']==sha(a.producer)
    assert reader['projection_cases']==17500 and reader['branch_projection_cells']==9047500
    bindings={str(a.producer):sha(a.producer),str(a.readback):sha(a.readback),str(Path(__file__)):sha(__file__)}
    for record in [producer,reader]:
        for field in ['source_hashes','artifacts']:
            for path,digest in record.get(field,{}).items():
                assert sha(path)==digest,path
                assert path not in bindings or bindings[path]==digest,path
                bindings[path]=digest
    roots={Path(path).parent for path in producer['artifacts']};assert len(roots)==1
    root=roots.pop();cases=rows(root/'projection_summary.tsv')
    taxa=rows(root/'taxon_identity_source_coverage.tsv')
    axes=json.loads((root/'input_axes.json').read_text())
    assert len(cases)==17500 and len(taxa)==526 and len(axes['views'])==70
    assert sum(r['study_role']=='outgroup' for r in taxa)==25
    grouped=defaultdict(list)
    for row in cases:grouped[int(row['view_index']),row['source']].append(row)
    view_rows=[]
    for (vi,source),entries in sorted(grouped.items()):
        assert len(entries)==125 and {r['marker'] for r in entries}==set(axes['markers'])
        view=axes['views'][vi];n=len(view['branches']);denominator=125*n
        totals={status:sum(int(r[status]) for r in entries) for status in STATUSES}
        assert sum(totals.values())==denominator
        for entry in entries:assert sum(int(entry[s]) for s in STATUSES)==n
        view_rows.append(dict(view_index=vi,cohort=view['cohort'],view=view['view'],source=source,
            cohort_taxa=len(view['taxa']),cohort_outgroups=view['roles']['outgroup'],
            original_internal_branches=n,marker_slots=125,branch_marker_cells=denominator,
            markers_with_four_or_more_observed_taxa=sum(int(r['retained_tips'])>=4 for r in entries),
            minimum_retained_tips=min(int(r['retained_tips']) for r in entries),
            maximum_retained_tips=max(int(r['retained_tips']) for r in entries),
            **totals,**{status+'_fraction':totals[status]/denominator for status in STATUSES}))
    assert len(view_rows)==140 and sum(r['branch_marker_cells'] for r in view_rows)==9047500
    bins=defaultdict(list)
    for row in taxa:bins[row['study_role'],row['manifest_lineage_group']].append(row)
    lineage_rows=[]
    for (role,lineage),entries in sorted(bins.items()):
        for source in SOURCES:
            counts=[int(r[source+'_usable_markers']) for r in entries]
            lineage_rows.append(dict(study_role=role,manifest_lineage_group=lineage,source=source,
                panel_entries=len(entries),usable_marker_taxon_cells=sum(counts),
                **{'entries_at_least_'+str(k)+'_markers':sum(v>=k for v in counts) for k in [1,10,50,100]},
                entries_with_identity_review_flags=sum(bool(r['identity_review_flags']) for r in entries)))
    cohorts=sorted({r['cohort'] for r in view_rows});assert len(cohorts)==5
    source_summary={}
    for source in SOURCES:
        selected=[r for r in view_rows if r['source']==source and r['cohort']=='full_primary']
        assert len(selected)==14
        source_summary[source]={status+'_fraction_range_across_primary_views':
            [min(r[status+'_fraction'] for r in selected),max(r[status+'_fraction'] for r in selected)]
            for status in STATUSES}
    a.output.mkdir(parents=True,exist_ok=False)
    table(a.output/'all_view_source_branch_coverage.tsv',view_rows)
    table(a.output/'all_manifest_lineage_source_coverage.tsv',lineage_rows)
    colors=['#737373','#d95f02','#e6ab02','#1b9e77','#7570b3']
    labels=['Fewer than four taxa','One side unobserved','Becomes terminal','Distinct internal branch','Shared internal path']
    fig,panels=plt.subplots(5,2,figsize=(14,19),layout='constrained')
    for ci,cohort in enumerate(cohorts):
        for si,source in enumerate(SOURCES):
            ax=panels[ci,si];selected=[r for r in view_rows if r['source']==source and r['cohort']==cohort]
            assert len(selected)==14
            y=np.arange(14);left=np.zeros(14)
            for status,color,label in zip(STATUSES,colors,labels):
                values=np.array([r[status+'_fraction'] for r in selected])*100
                ax.barh(y,values,left=left,color=color,label=label,height=.8);left+=values
            ax.set_yticks(y,[r['view'].replace('coal:','').replace('ref:','') for r in selected],fontsize=7)
            ax.invert_yaxis();ax.set_xlim(0,100)
            ax.set_title(source+' | '+cohort.replace('_',' ')+'\n'+str(selected[0]['cohort_taxa'])+
                         ' taxa, '+str(selected[0]['cohort_outgroups'])+' outgroups',fontsize=9)
            ax.set_xlabel('Original internal branch × marker entries (%)',fontsize=8)
    handles,legend_labels=panels[0,0].get_legend_handles_labels()
    fig.legend(handles,legend_labels,loc='outside upper center',ncol=3,fontsize=9)
    for extension in ['pdf','png']:fig.savefig(a.output/('all_tree_view_coverage.'+extension),dpi=160)
    plt.close(fig)
    for path,digest in bindings.items():assert sha(path)==digest,path
    result=dict(status='complete_descriptive_full_tree_view_structural_coverage_summary',
        checked_utc=datetime.now(timezone.utc).isoformat(),tree_views=70,view_source_summaries=140,
        full_panel_entries=526,full_panel_outgroups=25,primary_view_fraction_ranges=source_summary,
        lineage_bin_source_summaries=len(lineage_rows),source_hashes=bindings,
        artifacts={str(path):sha(path) for path in a.output.iterdir()},scientific_eligibility=False,
        scope='Descriptive coverage of every qualified view/source/marker projection after full raw-tree readback. All five cohorts,70views,125slots,526taxon rows and25primaryoutgroups retained. Repeated views and marker entries are not independent observations; bins are manifest labels, not established clades/species counts. No pooled predictor effects, fitted rates, statistical estimability, accepted tree/root/model or biological aim completion.')
    with a.receipt.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts','scope']},indent=2))


if __name__=='__main__':main()
