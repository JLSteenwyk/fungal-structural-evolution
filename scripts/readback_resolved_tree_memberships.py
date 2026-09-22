#!/usr/bin/env python3
"""Check complete native resolved-tree membership without recursion or pair expansion."""
import argparse
import hashlib
from io import StringIO
import json
import math
from pathlib import Path
from Bio import Phylo
from assess_small_family_output_exposure import groups, sha


def tree_membership(newick, labels):
    tree=Phylo.read(StringIO(newick),'newick')
    stack=[tree.root];genes=[];nodes=0
    while stack:
        node=stack.pop();nodes+=1
        if node is not tree.root and node.branch_length is None:
            raise ValueError('Missing resolved branch length')
        if node.branch_length is not None and (not math.isfinite(node.branch_length) or node.branch_length<0):
            raise ValueError('Invalid resolved branch length')
        if node.clades:stack.extend(node.clades)
        else:
            if node.name not in labels:raise ValueError('Unrecognized resolved leaf label')
            genes.append(labels[node.name])
    if len(set(genes))!=len(genes):raise ValueError('Repeated resolved leaf')
    return len(genes),hashlib.sha256(('\n'.join(sorted(genes))+'\n').encode()).hexdigest(),nodes


def audit(source, result, output):
    if output.exists():raise FileExistsError(output)
    inputs={str(source/name):sha(source/name) for name in ['SpeciesIDs.txt','SequenceIDs.txt','clusters_OrthoFinder.txt_id_pairs.txt']}
    species={}
    with (source/'SpeciesIDs.txt').open() as f:
        for line in f:
            native,name=line.rstrip('\n').split(': ',1)
            if native in species:raise ValueError('Repeated species ID')
            species[native]=name.rsplit('.',1)[0].replace('.','_').replace(' ','_')
    if len(set(species.values()))!=len(species):raise ValueError('Ambiguous native species labels')
    expected={};needed=set();families=proteins=0
    for family,genes in groups(source/'clusters_OrthoFinder.txt_id_pairs.txt'):
        families+=1;proteins+=len(genes)
        if len(genes)>=4:
            digest=hashlib.sha256(('\n'.join(sorted(genes))+'\n').encode()).hexdigest()
            expected[family]=(len(genes),digest)
            if needed.intersection(genes):raise ValueError('Gene appears in multiple tree families')
            needed.update(genes)
    labels={}
    with (source/'SequenceIDs.txt').open() as f:
        for line in f:
            native,name=line.rstrip('\n').split(': ',1)
            if native not in needed:continue
            if len(name.split())!=1:raise ValueError('Full-header fallback requires separate naming review')
            for c in ':,()':name=name.replace(c,'_')
            label=species[native.split('_')[0]]+'_'+name
            if label in labels:raise ValueError('Ambiguous resolved leaf label')
            labels[label]=native
    if set(labels.values())!=needed:raise ValueError('Incomplete label mapping')
    del needed
    path=result/'Resolved_Gene_Trees/Resolved_Gene_Trees.txt'
    observed=set();tip_count=node_count=0;digest=hashlib.sha256();before=path.stat()
    with path.open('rb') as f:
        for raw in f:
            digest.update(raw)
            family,newick=raw.decode().rstrip('\r\n').split(': ',1)
            if family in observed or family not in expected:raise ValueError('Unexpected or repeated resolved family')
            tips,membership,nodes=tree_membership(newick,labels)
            if (tips,membership)!=expected[family]:raise ValueError('Resolved tree changes family membership: '+family)
            observed.add(family);tip_count+=tips;node_count+=nodes
            if len(observed)%5000==0:print('Checked',len(observed),'resolved trees',flush=True)
    after=path.stat()
    if (before.st_size,before.st_mtime_ns,before.st_ino)!=(after.st_size,after.st_mtime_ns,after.st_ino):
        raise ValueError('Resolved output changed during readback')
    if observed!=set(expected):raise ValueError('Missing resolved families: '+str(len(set(expected)-observed)))
    for path,value in inputs.items():
        if sha(path)!=value:raise ValueError('Source changed during audit')
    receipt=dict(status='passed_complete_resolved_tree_membership_readback',source_families=families,source_proteins=proteins,
                 resolved_trees=len(observed),resolved_tips=tip_count,resolved_nodes=node_count,input_hashes=inputs,
                 resolved_tree_file_sha256=digest.hexdigest(),script_sha256=sha(__file__),cluster_parser_sha256=sha(Path(__file__).with_name('assess_small_family_output_exposure.py')),
                 scope='Every native resolved tree for families with at least four genes checked for exact source membership, unique labels and finite nonnegative branch lengths. Native min_seq=4 output contract; three-gene input trees are not expected as resolved outputs. Does not independently infer rooting/topology, reconciliation events, HOGs, ortholog pairs or biological homology.')
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(receipt,indent=2)+'\n');return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','results','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();print(json.dumps(audit(a.source,a.results,a.output),indent=2))


if __name__=='__main__':main()
