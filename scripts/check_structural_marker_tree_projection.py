#!/usr/bin/env python3
"""Exhaustive synthetic pruning and serialized-array contracts before execution."""
import argparse
import copy
from datetime import datetime,timezone
import json
from pathlib import Path

import dendropy
import numpy as np

from audit_selected_taxon_identity_snapshot_v2 import sha
from structural_marker_tree_projection import load,project,mask_for,DTYPES,STATUSES
from readback_structural_marker_tree_projection import tree_edges,check_slice,read_arrays


def reject(action):
    try: action()
    except (AssertionError,ValueError,KeyError): return
    raise AssertionError('Altered projection accepted')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists() and not a.receipt.exists()
    a.output.mkdir(parents=True,exist_ok=False)
    newick='((a:1,b:1):2,(c:1,d:1):3,((e:1,f:1):4,(g:1,(h:1,i:1):5):6):7);'
    tree=dendropy.Tree.get(data=newick,schema='newick',rooting='force-unrooted')
    positions={t:i for i,t in enumerate('abcdefghi')}
    tips,full,edges=tree_edges(tree,positions)
    assert len(tips)==9 and len(edges)==6
    branches=[dict(split_mask_hex=hex(mask),branch_length=length,branch_length_unit='synthetic_units')
              for mask,length in sorted(edges.items())]
    seen=set();cases=[]
    for retained_mask in range(512):
        retained={t for t,index in positions.items() if retained_mask&(1<<index)}
        values,split_count=project(branches,retained_mask)
        counts,independent_count=check_slice(branches,retained,positions,tree,values)
        assert split_count==independent_count
        seen.update(counts)
        cases.append(dict(retained_mask=retained_mask,retained_tips=len(retained),
            actual_internal_splits=split_count,status_counts={STATUSES[k]:v for k,v in counts.items()}))
    assert seen==set(range(5))
    values,_=project(branches,full)
    negatives=[]
    for field in DTYPES:
        bad={k:v.copy() for k,v in values.items()}
        bad[field][0]+=1
        reject(lambda:check_slice(branches,tips,positions,tree,bad));negatives.append('altered_'+field)
    none,_=project(branches,0)
    bad={k:v.copy() for k,v in none.items()};bad['path_length_sum'][0]=0
    reject(lambda:check_slice(branches,set(),positions,tree,bad));negatives.append('invented_missing_path_length')
    serialized={key:np.tile(value,(2,125,1)) for key,value in values.items()}
    good=a.output/'good.npz'
    with good.open('xb') as h:np.savez_compressed(h,**serialized)
    decoded=read_arrays(good,(2,125,6))
    assert all(np.array_equal(decoded[k],v,equal_nan=True) for k,v in serialized.items())
    for case in ['missing_key','extra_key','wrong_dtype','wrong_shape','object_dtype']:
        bad={k:v.copy() for k,v in serialized.items()}
        if case=='missing_key':bad.pop('status')
        elif case=='extra_key':bad['invented']=bad['status']
        elif case=='wrong_dtype':bad['status']=bad['status'].astype('uint16')
        elif case=='wrong_shape':bad['status']=bad['status'][:,:124]
        else:bad['status']=bad['status'].astype(object)
        path=a.output/(case+'.npz')
        with path.open('xb') as h:np.savez_compressed(h,**bad)
        reject(lambda:read_arrays(path,(2,125,6)));negatives.append('serialized_'+case)
    views,positions,markers,datasets,taxa,bindings=load()
    assert len(views)==70 and len(taxa)==526 and len(markers)==125
    cases_path=a.output/'all_synthetic_pruning_cases.json'
    with cases_path.open('x') as h:h.write(json.dumps(cases,indent=2)+'\n')
    for path in [Path(__file__),Path('scripts/structural_marker_tree_projection.py'),
        Path('scripts/project_structural_markers_on_candidate_trees.py'),
        Path('scripts/readback_structural_marker_tree_projection.py')]:bindings[str(path)]=sha(path)
    for path,digest in bindings.items():assert sha(path)==digest,path
    result=dict(status='passed_full_structural_marker_tree_projection_software_and_source_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(),synthetic_retained_sets=512,
        synthetic_original_branches=6,synthetic_original_branch_cells=3072,
        all_five_status_classes_observed=True,altered_cases_rejected=negatives,
        full_real_source_views=70,full_real_source_cohorts=5,full_real_marker_slots=125,
        full_real_taxon_entries=526,full_original_internal_branches=36190,
        future_projection_cases=17500,future_branch_projection_cells=9047500,
        source_hashes=bindings,artifacts={str(path):sha(path) for path in a.output.iterdir()},
        scientific_eligibility=False,
        scope='All512subsets of a nine-tip unrooted synthetic tree exhaustively compared with independent DendroPy actual pruning/path lengths; allfiveprojectionclasses and11altered value/serialization cases checked. Complete actual70view/two-source/125marker source metadata, original branch identities and every emitted paired-alignment mask/membership verified. Synthetic checks are not production fungal projection, a biological pilot, inference, estimability/model/root/predictor qualification or any biological aim completion.')
    with a.receipt.open('x') as h:h.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2))


if __name__=='__main__':main()
