#!/usr/bin/env python3
"""Independently prune raw candidate trees and check every projected array cell."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
import math
from pathlib import Path

import dendropy
import numpy as np

from audit_selected_taxon_identity_snapshot_v2 import rows,sha
from structural_marker_tree_projection import load,SOURCES,STATUSES,DTYPES


def tree_edges(tree,positions):
    """Raw DendroPy graph traversal; no producer split/project routine."""
    tree.deroot();tree.encode_bipartitions()
    tips={node.taxon.label for node in tree.leaf_node_iter()}
    full=sum(1<<positions[t] for t in tips)
    descendants={}
    for node in tree.postorder_node_iter():
        descendants[node]=(1<<positions[node.taxon.label] if node.is_leaf()
                           else sum(descendants[c] for c in node.child_node_iter()))
    edges={}
    for node in tree.preorder_node_iter():
        if node is tree.seed_node or node.is_leaf():continue
        a,b=descendants[node],full^descendants[node]
        key=a if (a.bit_count(),a)<=(b.bit_count(),b) else b
        assert key not in edges and min(a.bit_count(),b.bit_count())>=2
        assert node.edge.length is not None and math.isfinite(node.edge.length) and node.edge.length>=0
        edges[key]=node.edge.length
    assert len(edges)==len(tips)-3
    return tips,full,edges


def check_slice(branches,retained,positions,base_tree,saved):
    """Check all original branches against actual graph pruning/path lengths."""
    n=len(retained)
    mask=sum(1<<positions[t] for t in retained)
    projected_edges={}
    if n>=4:
        pruned=base_tree.clone(depth=2)
        pruned.retain_taxa_with_labels(sorted(retained))
        actual_tips,actual_mask,projected_edges=tree_edges(pruned,positions)
        assert actual_tips==retained and actual_mask==mask
    contributors=Counter()
    keys=[]
    for branch in branches:
        side=int(branch['split_mask_hex'],16)&mask
        other=mask^side
        if min(side.bit_count(),other.bit_count())>=2:
            key=side if (side.bit_count(),side)<=(other.bit_count(),other) else other
            contributors[key]+=1;keys.append(key)
        else:keys.append(None)
    if n>=4:assert set(contributors)==set(projected_edges)
    counts=Counter()
    for index,branch in enumerate(branches):
        a=(int(branch['split_mask_hex'],16)&mask).bit_count();b=n-a
        assert int(saved['observed_side'][index])==a and int(saved['observed_complement'][index])==b
        if n<4:status=0
        elif min(a,b)==0:status=1
        elif min(a,b)==1:status=2
        else:status=3 if contributors[keys[index]]==1 else 4
        assert int(saved['status'][index])==status
        counts[status]+=1
        if status in [3,4]:
            assert int(saved['path_internal_branch_count'][index])==contributors[keys[index]]
            assert math.isclose(float(saved['path_length_sum'][index]),projected_edges[keys[index]],
                                rel_tol=1e-10,abs_tol=1e-10)
        else:
            assert int(saved['path_internal_branch_count'][index])==0
            assert np.isnan(saved['path_length_sum'][index])
    return counts,len(projected_edges)


def read_arrays(path,shape):
    with np.load(path,allow_pickle=False) as data:
        assert set(data.files)==set(DTYPES)
        arrays={key:data[key] for key in DTYPES}
        for key,dtype in DTYPES.items():
            assert arrays[key].dtype==np.dtype(dtype) and arrays[key].shape==shape
    return arrays


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--producer',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    producer=json.loads(a.producer.read_text());producer_sha=sha(a.producer)
    assert producer['status']=='complete_full_structural_marker_tree_projection_pending_independent_readback'
    for path,digest in producer['source_hashes'].items():assert sha(path)==digest,path
    for path,digest in producer['artifacts'].items():assert sha(path)==digest,path
    roots={Path(path).parent for path in producer['artifacts']};assert len(roots)==1
    root=roots.pop();axes=json.loads((root/'input_axes.json').read_text())
    views,positions,markers,datasets,taxa,bindings=load()
    assert axes['positions']==positions and axes['markers']==markers and axes['views']==views
    assert axes['sources']==list(SOURCES) and axes['statuses']==STATUSES
    assert rows(root/'taxon_identity_source_coverage.tsv')==[{k:str(v) for k,v in r.items()} for r in taxa]
    for source in SOURCES:
        assert axes['eligible_taxa'][source]=={m:sorted(datasets[source]['eligible'][m]) for m in markers}
    summary=rows(root/'projection_summary.tsv');assert len(summary)==17500
    index={(int(r['view_index']),r['source'],r['marker']):r for r in summary}
    assert len(index)==len(summary)
    totals=Counter();cases=0
    for view_index,view in enumerate(views):
        base=dendropy.Tree.get(path=view['source_tree'],schema='newick',rooting='force-unrooted',preserve_underscores=True)
        before={node.taxon.label for node in base.leaf_node_iter()}
        assert before==(set(positions) if view['prune'] else set(view['taxa']))
        if view['prune']:base.retain_taxa_with_labels(view['taxa'])
        actual_tips,full,actual_edges=tree_edges(base,positions)
        assert actual_tips==set(view['taxa']) and hex(full)==view['full_mask_hex']
        assert set(actual_edges)=={int(b['split_mask_hex'],16) for b in view['branches']}
        for branch in view['branches']:
            assert math.isclose(actual_edges[int(branch['split_mask_hex'],16)],branch['branch_length'],
                                rel_tol=1e-12,abs_tol=1e-12)
        arrays=read_arrays(root/(f'view_{view_index:02d}.npz'),(2,125,len(view['branches'])))
        for si,(source,dataset) in enumerate(datasets.items()):
            for mi,marker in enumerate(markers):
                kept=set(dataset['eligible'][marker])&actual_tips
                counts,split_count=check_slice(view['branches'],kept,positions,base,
                    {key:value[si,mi] for key,value in arrays.items()})
                row=index[view_index,source,marker];cases+=1;totals.update(counts)
                for key,value in [('cohort',view['cohort']),('view',view['view']),('family',view['family'])]:assert row[key]==value
                assert int(row['retained_tips'])==len(kept)
                assert int(row['retained_outgroups'])==sum(t in kept for t in [r['taxon_id'] for r in taxa if r['study_role']=='outgroup'])
                assert int(row['original_internal_branches'])==len(view['branches'])
                assert int(row['projected_internal_splits'])==split_count
                for code,status in enumerate(STATUSES):assert int(row[status])==counts[code]
        print('independent_full_raw_tree_projection_verified',view_index,view['cohort'],view['view'],flush=True)
    assert cases==17500 and sum(totals.values())==9047500
    assert producer['status_counts']=={STATUSES[i]:totals[i] for i in range(len(STATUSES))}
    for path,digest in producer['source_hashes'].items():assert sha(path)==digest,path
    for path,digest in producer['artifacts'].items():assert sha(path)==digest,path
    assert sha(a.producer)==producer_sha
    bindings[str(a.producer)]=producer_sha
    bindings[str(Path(__file__))]=sha(__file__)
    bindings['scripts/structural_marker_tree_projection.py']=sha('scripts/structural_marker_tree_projection.py')
    result=dict(status='passed_full_structural_marker_raw_tree_projection_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),producer=str(a.producer),producer_sha256=producer_sha,
        tree_views=70,cohorts=5,projection_cases=cases,branch_projection_cells=sum(totals.values()),
        status_counts=producer['status_counts'],source_hashes=bindings,dendropy_version=dendropy.__version__,
        scientific_eligibility=False,
        scope='Every source-qualified eligibility grid and emitted paired-alignment membership bound; all raw candidate trees independently parsed, pruned to exact declared cohort then source-marker memberships, and their actual projected split/path lengths checked against every one of9,047,500array cells. All125slots and17500cases retained. Unit-specific arithmetic path sums are not new fitted rates; unique projection is not statistical estimability. No accepted root/tree/model/likelihood, independent parental purity, predictor calibration, native restart/GPU/new charges or biological aim completion.')
    with a.output.open('x') as h:h.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','scope']},indent=2))


if __name__=='__main__':main()
