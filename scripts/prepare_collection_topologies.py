#!/usr/bin/env python3
"""Export exact native AA topologies from a verified per-marker fit collection."""
import argparse
import json
import os
from pathlib import Path
from Bio import Phylo, SeqIO
from paired_collection_resampling_sources import load_sources
from audit_busco_gene_copies import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','models','collection','inventory','output']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    ready,native,executable,bindings=load_sources(a.inputs,a.models,a.collection,a.inventory)
    trees={};sites=0;observations=0
    for row in ready:
        marker=row['marker'];path=native[marker]/'aa.treefile'
        seqs=list(SeqIO.parse(a.inputs/marker/'aa.faa','fasta'))
        tips=[x.name for x in Phylo.read(path,'newick').get_terminals()]
        if len(tips)!=len(set(tips)) or len(seqs)!=len({x.id for x in seqs}) or set(tips)!={x.id for x in seqs}:raise ValueError('Tree/sequence identity differs')
        if any(len(x.seq)!=int(row['retained_columns']) for x in seqs):raise ValueError('Alignment width differs')
        if sha(path)!=bindings[marker]['topology_sha256']:raise ValueError('Tree differs')
        trees[marker]=dict(path=str(path),sha256=sha(path),taxa=len(tips),sites=int(row['retained_columns']),source=bindings[marker])
        sites+=int(row['retained_columns']);observations+=sum(sum(c!='?' for c in str(x.seq)) for x in seqs)
    a.output.mkdir(parents=True)
    for marker,t in trees.items():
        folder=a.output/marker;folder.mkdir();(folder/'aa.treefile').symlink_to(os.path.relpath(t['path'],folder.resolve()))
    (a.output/'tree_sources.json').write_text(json.dumps(trees,indent=2)+'\n')
    artifacts={'tree_sources.json':sha(a.output/'tree_sources.json')}
    artifacts.update({marker+'/aa.treefile':t['sha256'] for marker,t in trees.items()})
    receipt=dict(status='complete_verified_collection_topologies',markers=len(ready),sites=sites,taxon_site_observations=observations,source_receipts={k:sha(getattr(a,k)/'receipt.json') for k in ['inputs','models','collection']},inventory_sha256=sha(a.inventory),script_sha256=sha(Path(__file__)),source_helper_sha256=sha(Path(__file__).with_name('paired_collection_resampling_sources.py')),artifacts=artifacts,scope='Relative symlinks to exact native AA point-estimate topologies selected by the audited collection inventory. No topology inference or uncertainty integration performed.')
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='artifacts'},indent=2))


if __name__=='__main__':main()
