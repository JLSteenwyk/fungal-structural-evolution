#!/usr/bin/env python3
"""Compare every marker topology on explicitly shared alignment-eligible taxa."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from Bio import Phylo


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def table(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def projected_splits(tree, common):
    result=set()
    for node in tree.get_nonterminals():
        side={t.name for t in node.get_terminals()} & common
        other=common-side
        if min(len(side),len(other))<2:continue
        result.add(min(tuple(sorted(side)),tuple(sorted(other)),key=lambda x:(len(x),x)))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for label in ['first','second']:
        for key in ['trees','support','audit']:p.add_argument('--'+label+'-'+key,type=Path,required=True)
        p.add_argument('--'+label+'-coverage',choices=['profile','mafft'],required=True)
    p.add_argument('--coverage',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    coverage_receipt=json.loads((a.coverage/'receipt.json').read_text())
    for name,d in coverage_receipt['artifacts'].items():assert sha(a.coverage/name)==d
    coverage={r['marker']:r for r in table(a.coverage/'marker_summary.tsv')};assert len(coverage)==125
    sources={};collections={}
    for label in ['first','second']:
        trees=getattr(a,label+'_trees');support=getattr(a,label+'_support');audit=getattr(a,label+'_audit')
        sr=json.loads((support/'receipt.json').read_text());ar=json.loads((audit/'receipt.json').read_text())
        assert ar['status']=='passed_all_snapshot_input_and_graph_split_readbacks' and ar['markers']==125
        assert ar['snapshot_receipt_sha256']==sha(support/'receipt.json')
        assert sr['completed_markers']==sr['planned_markers']==125 and not sr['pending_markers']
        for name,d in sr['artifacts'].items():assert sha(support/name)==d
        inputs={r['marker']:r for r in sr['inputs']};assert len(inputs)==len(sr['inputs'])==125 and set(inputs)==set(coverage)
        sources[label]=dict(support_receipt_sha256=sha(support/'receipt.json'),audit_receipt_sha256=sha(audit/'receipt.json'))
        collections[label]={}
        for marker,proof in inputs.items():
            folder=trees/marker
            assert sha(folder/'receipt.json')==proof['receipt_sha256']
            assert sha(folder/'tree.treefile')==proof['tree_sha256']
            tree=Phylo.read(folder/'tree.treefile','newick');tips=[t.name for t in tree.get_terminals()]
            c=coverage[marker];method=getattr(a,label+'_coverage')
            expected=set(filter(None,c['common_taxa'].split(';'))) | set(filter(None,c[method+'_only_taxa'].split(';')))
            assert len(tips)==len(set(tips)) and set(tips)==expected
            collections[label][marker]=(tree,set(tips))
    summaries=[];splits=[]
    for marker in sorted(coverage):
        first,ft=collections['first'][marker];second,st=collections['second'][marker];common=ft&st
        fs=projected_splits(first,common);ss=projected_splits(second,common)
        assert len(fs)==len(ss)==len(common)-3
        summaries.append(dict(marker=marker,first_taxa=len(ft),second_taxa=len(st),shared_taxa=len(common),first_only_taxa=';'.join(sorted(ft-st)),second_only_taxa=';'.join(sorted(st-ft)),shared_internal_splits=len(fs&ss),first_only_internal_splits=len(fs-ss),second_only_internal_splits=len(ss-fs),unrooted_rf=len(fs^ss),topology_changed=fs!=ss))
        for side in sorted(fs|ss):splits.append(dict(marker=marker,split_taxa_json=json.dumps(side),in_first=side in fs,in_second=side in ss))
    a.output.mkdir(parents=True)
    for name,rows in [('marker_comparison.tsv',summaries),('split_comparison.tsv',splits)]:
        with (a.output/name).open('w') as f:
            w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    r=dict(status='complete_full_marker_common_taxon_topology_comparison_pending_independent_readback',markers=125,changed_topologies=sum(r['topology_changed'] for r in summaries),split_rows=len(splits),sources=sources,coverage_receipt_sha256=sha(a.coverage/'receipt.json'),script_sha256=sha(Path(__file__)),artifacts={p.name:sha(p) for p in a.output.iterdir()},scope='Unrooted topology comparison after restricting each split to shared taxa. Pruning is not refitting and does not remove inference effects of different taxon sampling. Zero-length edges remain represented. No likelihood, support significance, preferred alignment or biological-cause claim.')
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))


if __name__=='__main__':main()
