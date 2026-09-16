#!/usr/bin/env python3
"""Read back all branch-frame lengths by undirected graph edge deletion."""
import argparse
import json
from pathlib import Path
from collections import defaultdict
from Bio import Phylo
import numpy as np
import pandas as pd
from audit_joint_path_uncertainty import checked,sha


def graph_edges(path):
    tree=Phylo.read(path,'newick');nodes=list(tree.find_clades());number={c:i for i,c in enumerate(nodes)}
    graph=defaultdict(list);edges=[];leaves={number[c]:c.name for c in tree.get_terminals()};universe=set(leaves.values())
    for c in nodes:
        for child in c.clades:
            u,v=number[c],number[child];graph[u].append(v);graph[v].append(u);edges.append((u,v,child.branch_length))
    result=defaultdict(float)
    for u,v,length in edges:
        visited={u};pending=[v];side=set()
        while pending:
            n=pending.pop()
            if n in visited:continue
            visited.add(n)
            if n in leaves:side.add(leaves[n])
            pending.extend(graph[n])
        alternatives=[tuple(sorted(side)),tuple(sorted(universe-side))]
        canonical=sorted(alternatives,key=lambda x:(len(x),x))[0]
        result[','.join(canonical)]+=length
    return dict(result)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--frame',type=Path,required=True);ap.add_argument('--projection',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=checked(a.frame);checked(a.projection)
    for path,digest in r['source_pins'].items():
        if sha(Path(path))!=digest:raise ValueError('Changed pinned source '+path)
    d=pd.read_csv(a.frame/'branch_frame.tsv',sep='\t',keep_default_na=False)
    keys=['cohort','marker','marker_split_taxa'];assert not d.duplicated(keys).any()
    expected=pd.read_csv(a.projection/'guide_sensitivity.tsv',sep='\t',keep_default_na=False).set_index(keys).sort_index()
    actual=d.set_index(keys).sort_index();assert actual.index.equals(expected.index)
    np.testing.assert_array_equal(actual.projection_status,expected.status)
    np.testing.assert_array_equal(actual.marker_terminal,expected.marker_terminal)
    projections=pd.read_csv(a.projection/'edge_projection.tsv',sep='\t',keep_default_na=False)
    for row in projections.to_dict('records'):
        item=actual.loc[tuple(row[k] for k in keys)]
        assert json.loads(item.guide_edge_sets_json)[row['guide']]==json.loads(row['full_edge_ids_json'])
    lengths=0;trees=0
    for source in r['sources']:
        subset=d[d.cohort==source['cohort']]
        for marker,group in subset.groupby('marker'):
            for label in ['aa','3di_af','3di_af_empirical','3di_llm']:
                edges=graph_edges(Path(source['rates'])/marker/(label+'.treefile'))
                assert set(edges)==set(group.marker_split_taxa)
                np.testing.assert_allclose(group[label+'_branch_length'],[edges[x] for x in group.marker_split_taxa],rtol=1e-12,atol=1e-12)
                lengths+=len(edges);trees+=1
    selected=d[d.unique_full_edge_id!=''];coverage=pd.read_csv(a.frame/'unique_edge_coverage.tsv',sep='\t')
    grouped=selected.groupby(['cohort','unique_full_edge_id'])
    assert len(coverage)==len(grouped)
    for row in coverage.to_dict('records'):
        group=grouped.get_group((row['cohort'],row['full_edge_id']))
        assert group.marker.nunique()==len(group)==row['markers']
        assert sorted(group.marker)==json.loads(row['marker_ids_json'])
        for mapping in group.guide_edge_sets_json:
            assert all(v==[row['full_edge_id']] for v in json.loads(mapping).values())
    result={'status':'passed_projected_branch_frame_readback','frame_rows':len(d),'fitted_trees_checked':trees,'branch_lengths_checked':lengths,'projection_rows_checked':len(projections),'unique_coverage_rows_checked':len(coverage),'frame_receipt_sha256':sha(a.frame/'receipt.json'),'script_sha256':sha(Path(__file__)),'scope':'Every numerical branch length reconstructed by independent undirected graph traversal; full projection key/status/mapping preservation and unique-edge marker coverage checked. Does not validate model biology, supported topology or acceleration.'}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
