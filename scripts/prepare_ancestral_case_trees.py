#!/usr/bin/env python3
"""Prepare local full/domain tree inputs and verify every induced edge length."""
import copy,csv,json,math
from collections import defaultdict
from pathlib import Path
from Bio import Phylo,SeqIO
from prepare_case_ancestral_neighborhoods import sha,read


def descendants(tree):
    result={}
    for node in tree.find_clades(order='postorder'):
        result[node]=set().union(*(result[c] for c in node.clades)) if node.clades else {node.name}
    return result


def projected_edges(tree,keep):
    desc=descendants(tree);lengths=defaultdict(float)
    for node,genes in desc.items():
        selected=genes & keep
        if node is tree.root or not selected or selected==keep:continue
        assert node.branch_length is not None and math.isfinite(node.branch_length) and node.branch_length>=0
        lengths[tuple(sorted(selected))]+=node.branch_length
    return dict(lengths)


def main():
    root=Path('results/ancestral/case-neighborhoods-20260927-v1');pins={}
    def checked(base,name):
        rp=base/'receipt.json';r=json.loads(rp.read_text());p=base/name
        assert sha(p)==r['artifacts'][name];pins[str(rp)]=sha(rp);pins[str(p)]=sha(p);return p
    proof=Path('metadata/case_ancestral_neighborhoods_completed_20260927.json')
    assert json.loads(proof.read_text())['source_receipt_sha256']==sha(root/'receipt.json');pins[str(proof)]=sha(proof)
    neighborhoods=read(checked(root,'neighborhoods.tsv'))
    seqroot=Path('results/ancestral/case-sequences-20260927-v1');domainroot=Path('results/ancestral/case-domain-sequences-20260927-v1')
    out=Path('results/ancestral/case-local-trees-20260927-v1');out.mkdir(exist_ok=False)
    reports=[];nodes=[];edge_rows=[];tips_rows=[]
    for row in neighborhoods:
        if row['level']!='3':continue
        family,guide=row['family'],row['guide']
        source=Phylo.read(checked(root,row['tree_file']),'newick')
        original_nodes={n.name:n for n in source.find_clades()};local=Phylo.BaseTree.Tree(root=copy.deepcopy(original_nodes[row['node']]),rooted=True)
        local.root.branch_length=0.0
        full_path=checked(seqroot,family+'.faa');full={r.id for r in SeqIO.parse(full_path,'fasta')}
        assert {t.name for t in local.get_terminals()}==full
        domain={r.id for r in SeqIO.parse(checked(domainroot,family+'-alignment.faa'),'fasta')}
        assert domain=={r.id for r in SeqIO.parse(checked(domainroot,family+'-envelope.faa'),'fasta')}
        assert {row['gene_a'],row['gene_b']}<=domain<=full
        levels=[r for r in neighborhoods if r['guide']==guide and r['family']==family]
        source_desc=descendants(local)
        for dataset,keep in [('whole',full),('domain',domain)]:
            tree=copy.deepcopy(local)
            for tip in list(tree.get_terminals()):
                if tip.name not in keep:tree.prune(tip)
            tree.root.branch_length=0.0
            path=out/(guide+'-'+family+'-'+dataset+'.nwk')
            Phylo.write(tree,path,'newick',format_branch_length='%1.17g')
            recovered=Phylo.read(path,'newick');tiplist=[n.name for n in recovered.get_terminals()]
            assert len(tiplist)==len(keep)==len(set(tiplist)) and set(tiplist)==keep
            expected=projected_edges(local,keep);observed=projected_edges(recovered,keep)
            assert set(expected)==set(observed)
            error=max(abs(expected[k]-observed[k]) for k in expected)
            assert error<1e-12
            observed_desc=descendants(recovered)
            byset=defaultdict(list)
            for node,genes in observed_desc.items():byset[tuple(sorted(genes))].append(node.name)
            for level in levels:
                node=next(n for n in source_desc if n.name==level['node']);projected=source_desc[node]&keep
                matches=byset[tuple(sorted(projected))]
                assert matches
                nodes.append(dict(guide=guide,family=family,dataset=dataset,level=int(level['level']),source_node=node.name,source_descendants=len(source_desc[node]),retained_descendants=len(projected),output_nodes_json=json.dumps(matches),unique_output_node=int(len(matches)==1),retained_set_json=json.dumps(sorted(projected))))
            focal=recovered.common_ancestor(row['gene_a'],row['gene_b'])
            assert observed_desc[focal]=={row['gene_a'],row['gene_b']}
            reports.append(dict(guide=guide,family=family,dataset=dataset,proteins=len(keep),taxa=len({g.split('_',1)[0] for g in keep}),removed_proteins=len(full-keep),source_local_root=local.root.name,output_root=recovered.root.name,root_changed=int(local.root.name!=recovered.root.name),focal_node=focal.name,source_projected_edges=len(expected),maximum_edge_length_error=error,tree_file=path.name,sha256=sha(path)))
            for genes,value in sorted(expected.items()):edge_rows.append(dict(guide=guide,family=family,dataset=dataset,descendants_json=json.dumps(genes),expected_length=value,observed_length=observed[genes]))
            tips_rows.extend(dict(guide=guide,family=family,dataset=dataset,gene=g,retained=int(g in keep)) for g in sorted(full))
    assert len(reports)==52 and len(nodes)==208
    for filename,data in [('tree_inputs.tsv',reports),('ancestral_node_mapping.tsv',nodes),('edge_length_readback.tsv',edge_rows),('tip_dispositions.tsv',tips_rows)]:
        with (out/filename).open('w') as f:w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
        assert read(out/filename)==[{k:str(v) for k,v in row.items()} for row in data]
    for p,h in pins.items():assert sha(p)==h
    result=dict(status='complete_local_ancestral_tree_inputs_with_projected_edge_readback',trees=len(reports),node_mappings=len(nodes),edge_records=len(edge_rows),tip_dispositions=len(tips_rows),changed_roots=sum(r['root_changed'] for r in reports),ambiguous_node_mappings=sum(not r['unique_output_node'] for r in nodes),maximum_edge_length_error=max(r['maximum_edge_length_error'] for r in reports),source_hashes=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Local level-three clades under both reconciliation guides, full and domain tip sets. All induced rooted descendant sets and summed edge lengths verified after serialization. Root stem reset to zero; no contribution to within-clade distances. Source lengths are inherited gene-tree estimates, not refitted domain branch lengths, dates or posterior support. No ancestral sequence inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2))


if __name__=='__main__':main()
