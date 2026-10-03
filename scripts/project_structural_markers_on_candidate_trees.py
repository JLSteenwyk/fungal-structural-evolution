#!/usr/bin/env python3
"""Project complete source-specific marker memberships onto all70closed views."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
from pathlib import Path

import numpy as np

from audit_selected_taxon_identity_snapshot_v2 import sha,table
from structural_marker_tree_projection import load,project,mask_for,SOURCES,STATUSES,DTYPES


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True)
    p.add_argument('--resources',type=Path,required=True)
    p.add_argument('--software',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists() and not a.receipt.exists()
    resources=json.loads(a.resources.read_text())
    assert resources['tree_views']==70 and resources['projection_cases']==17500
    software=json.loads(a.software.read_text())
    assert software['status']=='passed_full_structural_marker_tree_projection_software_and_source_contracts'
    assert software['synthetic_retained_sets']==512 and software['all_five_status_classes_observed']
    assert len(software['altered_cases_rejected'])==11
    assert software['full_real_source_views']==70 and software['future_projection_cases']==17500
    assert software['future_branch_projection_cells']==9047500
    for field in ['source_hashes','artifacts']:
        for path,digest in software[field].items():assert sha(path)==digest,path
    assert software['source_hashes']['scripts/project_structural_markers_on_candidate_trees.py']==sha(__file__)
    views,positions,markers,datasets,taxa,bindings=load()
    bindings.update(software['source_hashes']);bindings.update(software['artifacts'])
    for path in [Path(__file__),Path('scripts/structural_marker_tree_projection.py'),a.resources,a.software]:
        bindings[str(path)]=sha(path)
    a.output.mkdir(parents=True,exist_ok=False)
    table(a.output/'taxon_identity_source_coverage.tsv',taxa)
    inputs=dict(positions=positions,markers=markers,sources=list(SOURCES),statuses=STATUSES,views=views,
        eligible_taxa={s:{m:sorted(d['eligible'][m]) for m in markers} for s,d in datasets.items()})
    with (a.output/'input_axes.json').open('x') as h:h.write(json.dumps(inputs,separators=(',',':'))+'\n')
    projections=[];totals=Counter()
    for view_index,view in enumerate(views):
        n=len(view['branches']);shape=(len(SOURCES),len(markers),n)
        arrays={k:np.zeros(shape,dtype=dtype) for k,dtype in DTYPES.items()}
        arrays['path_length_sum'].fill(np.nan)
        for source_index,(source,dataset) in enumerate(datasets.items()):
            for marker_index,marker in enumerate(markers):
                kept=set(dataset['eligible'][marker])&set(view['taxa'])
                projected,split_count=project(view['branches'],mask_for(kept,positions))
                for key,value in projected.items():arrays[key][source_index,marker_index]=value
                counts=Counter(int(v) for v in projected['status'])
                totals.update(counts)
                projections.append(dict(view_index=view_index,cohort=view['cohort'],view=view['view'],
                    family=view['family'],source=source,marker=marker,retained_tips=len(kept),
                    retained_outgroups=sum(t in kept for t in [r['taxon_id'] for r in taxa if r['study_role']=='outgroup']),
                    original_internal_branches=n,projected_internal_splits=split_count,
                    **{status:counts[i] for i,status in enumerate(STATUSES)}))
        output=a.output/(f'view_{view_index:02d}.npz')
        with output.open('xb') as h:np.savez_compressed(h,**arrays)
        with np.load(output,allow_pickle=False) as saved:
            assert set(saved.files)==set(DTYPES)
            for key,value in arrays.items():
                assert saved[key].dtype==value.dtype and saved[key].shape==shape
                assert np.array_equal(saved[key],value,equal_nan=True)
        print('all_marker_source_projections_written',view_index,view['cohort'],view['view'],flush=True)
    table(a.output/'projection_summary.tsv',projections)
    assert len(projections)==17500 and sum(totals.values())==9047500
    for path,digest in bindings.items():assert sha(path)==digest,path
    receipt=dict(status='complete_full_structural_marker_tree_projection_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),tree_views=70,cohorts=5,
        sources=2,markers=125,full_panel_entries=526,full_panel_outgroups=25,
        original_view_internal_branches=36190,projection_cases=17500,branch_projection_cells=9047500,
        status_counts={STATUSES[i]:totals[i] for i in range(len(STATUSES))},source_hashes=bindings,
        artifacts={str(path):sha(path) for path in a.output.iterdir()},
        resources_path=str(a.resources),scientific_eligibility=False,
        scope='All70closednative/coalescent/original/projectedviews acrossfivepredeclaredcohorts,all125marker slots andbothqualifiedpredictor inputs. Alloriginalinternalbranches retained in complete arrays, including insufficient-family, lost-side and terminal projections. Grouped internal paths are arithmetic sums in each original view unit, never physical displacement or a new fitted rate. Source eligibility/masks remain separate; taxonomy identity overlay applied only to display/provenance. No terminal-path length claim, root/model/mixing/estimability/ancestral/biological acceptance, taxon replacement, inference/GPU restart or charges. Independent raw-tree/pruning and full array readback remain required.')
    with a.receipt.open('x') as h:h.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_hashes','artifacts','scope']},indent=2))


if __name__=='__main__':main()
