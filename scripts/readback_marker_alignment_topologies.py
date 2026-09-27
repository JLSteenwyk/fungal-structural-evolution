#!/usr/bin/env python3
"""Independently prune DendroPy trees and verify every marker comparison row."""
import argparse,csv,hashlib,json
from pathlib import Path
import dendropy


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for side in ['first','second']:
        for key in ['trees','support','audit']:ap.add_argument('--'+side+'-'+key,type=Path,required=True)
        ap.add_argument('--'+side+'-coverage',choices=['profile','mafft'],required=True)
    for key in ['coverage','comparison','output']:ap.add_argument('--'+key,type=Path,required=True)
    a=ap.parse_args();assert not a.output.exists()
    receipt=json.loads((a.comparison/'receipt.json').read_text())
    assert receipt['status']=='complete_full_marker_common_taxon_topology_comparison_pending_independent_readback'
    for name,h in receipt['artifacts'].items():assert sha(a.comparison/name)==h
    assert sha(a.coverage/'receipt.json')==receipt['coverage_receipt_sha256']
    cr=json.loads((a.coverage/'receipt.json').read_text())
    for name,h in cr['artifacts'].items():assert sha(a.coverage/name)==h
    cov={r['marker']:r for r in rows(a.coverage/'marker_summary.tsv')}
    summaries=rows(a.comparison/'marker_comparison.tsv');summary={r['marker']:r for r in summaries}
    assert len(summary)==len(summaries)==len(cov)==receipt['markers']==125 and set(summary)==set(cov)
    actual={}
    for row in rows(a.comparison/'split_comparison.tsv'):
        key=row['marker'],tuple(json.loads(row['split_taxa_json']))
        assert key not in actual
        assert row['in_first'] in ['True','False'] and row['in_second'] in ['True','False']
        actual[key]=(row['in_first']=='True',row['in_second']=='True')
    trees={}
    for side in ['first','second']:
        support=getattr(a,side+'_support');audit=getattr(a,side+'_audit');folder=getattr(a,side+'_trees')
        assert sha(support/'receipt.json')==receipt['sources'][side]['support_receipt_sha256']
        assert sha(audit/'receipt.json')==receipt['sources'][side]['audit_receipt_sha256']
        sr=json.loads((support/'receipt.json').read_text());ar=json.loads((audit/'receipt.json').read_text())
        assert ar['status']=='passed_all_snapshot_input_and_graph_split_readbacks' and ar['markers']==125
        assert ar['snapshot_receipt_sha256']==sha(support/'receipt.json')
        assert sr['completed_markers']==sr['planned_markers']==125 and not sr['pending_markers']
        index={r['marker']:r for r in sr['inputs']};assert len(index)==len(sr['inputs'])==125 and set(index)==set(cov)
        trees[side]={}
        for marker,proof in index.items():
            path=folder/marker/'tree.treefile'
            assert sha(path)==proof['tree_sha256'] and sha(folder/marker/'receipt.json')==proof['receipt_sha256']
            tree=dendropy.Tree.get(path=str(path),schema='newick',rooting='force-unrooted',preserve_underscores=True)
            tips=[n.taxon.label for n in tree.leaf_node_iter()]
            c=cov[marker];method=getattr(a,side+'_coverage')
            expected=set(filter(None,c['common_taxa'].split(';')))|set(filter(None,c[method+'_only_taxa'].split(';')))
            assert len(tips)==len(set(tips)) and set(tips)==expected
            trees[side][marker]=(tree,set(tips))
    expected_rows={};changed=0
    for marker in sorted(cov):
        ft,st=trees['first'][marker][1],trees['second'][marker][1];common=ft&st;sets=[]
        for side in ['first','second']:
            tree=trees[side][marker][0];tree.retain_taxa_with_labels(common);tree.encode_bipartitions();splits=set()
            for edge in tree.postorder_edge_iter():
                left={n.taxon.label for n in edge.head_node.leaf_iter()};right=common-left
                if min(len(left),len(right))<2:continue
                splits.add(min(tuple(sorted(left)),tuple(sorted(right)),key=lambda x:(len(x),x)))
            assert len(splits)==len(common)-3;sets.append(splits)
        fs,ss=sets;changed+=fs!=ss
        expected=dict(marker=marker,first_taxa=len(ft),second_taxa=len(st),shared_taxa=len(common),first_only_taxa=';'.join(sorted(ft-st)),second_only_taxa=';'.join(sorted(st-ft)),shared_internal_splits=len(fs&ss),first_only_internal_splits=len(fs-ss),second_only_internal_splits=len(ss-fs),unrooted_rf=len(fs^ss),topology_changed=fs!=ss)
        assert summary[marker]=={k:str(v) for k,v in expected.items()}
        for split in fs|ss:expected_rows[marker,split]=(split in fs,split in ss)
    assert actual==expected_rows and len(actual)==receipt['split_rows'] and changed==receipt['changed_topologies']
    result=dict(status='passed_full_marker_topology_comparison_independent_pruning_readback',markers=125,split_rows=len(actual),changed_topologies=changed,source_receipt_sha256=sha(a.comparison/'receipt.json'),script_sha256=sha(Path(__file__)),scope='All source tree identities, coverage memberships, pruned split rows and marker summaries independently reconstructed with DendroPy. No refitting, preferred-method, branch-support significance or biological-cause claim.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
