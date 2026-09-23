#!/usr/bin/env python3
"""Verify every inherited geometry field and every tree path by graph traversal."""
import argparse
from collections import Counter
import csv
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
from Bio import Phylo, SeqIO


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def rows(path):
    with Path(path).open() as f:yield from csv.DictReader(f,delimiter='\t')


def fingerprint(row):
    return hashlib.sha256(json.dumps(row,sort_keys=True,separators=(',',':')).encode()).digest()


def graph_paths(path,taxa):
    tree=Phylo.read(path,'newick');nodes=list(tree.find_clades());adj={id(n):[] for n in nodes}
    leaves={n.name:id(n) for n in tree.get_terminals()}
    if set(leaves)!=set(taxa) or len(leaves)!=len(tree.get_terminals()):raise ValueError('Tree tip scope differs')
    for parent in nodes:
        for child in parent.clades:
            length=float(child.branch_length) if child.branch_length is not None else 0.
            if not math.isfinite(length) or length<0:raise ValueError('Invalid tree edge')
            a,b=id(parent),id(child);adj[a].append((b,length));adj[b].append((a,length))
    distances={}
    for i,a in enumerate(taxa):
        reached={};stack=[(leaves[a],None,0.)]
        while stack:
            node,parent,distance=stack.pop()
            if node in reached:raise ValueError('Cyclic tree graph')
            reached[node]=distance
            stack.extend((other,node,distance+length) for other,length in adj[node] if other!=parent)
        if len(reached)!=len(nodes):raise ValueError('Disconnected tree')
        for b in taxa[i+1:]:distances[a,b]=reached[leaves[b]]
    return distances


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ['inputs','fits','geometry','benchmark','review','output']:ap.add_argument('--'+key,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError('Use fresh output')
    pins={};receipts={}
    for key in ['inputs','fits','geometry','benchmark']:
        root=getattr(a,key);rp=root/'receipt.json';pins[str(rp)]=sha(rp);r=json.loads(rp.read_text());receipts[key]=r
        for name,h in r.get('artifacts',{}).items():pins[str(root/name)]=h
    def verify():
        for path,h in pins.items():
            if sha(path)!=h:raise ValueError('Changed artifact: '+path)
    pins[str(a.review)]=sha(a.review);verify()
    br=receipts['benchmark'];gr=receipts['geometry'];fr=receipts['fits'];ir=receipts['inputs']
    if (br['status']!='complete_paired_site_tree_path_point_benchmark'
            or fr['status']!='complete_matched_topology_point_estimates'
            or gr['status']!='complete_paired_site_geometry'
            or any(br['source_receipts'][k]!=sha(getattr(a,k)/'receipt.json') for k in ['inputs','fits','geometry'])
            or br['marker_review_sha256']!=sha(a.review)):
        raise ValueError('Benchmark/source lineage differs')
    review={r['marker']:r['status'] for r in json.loads(a.review.read_text())['records']}
    ready={r['marker'] for r in rows(a.inputs/'marker_summary.tsv') if r['status']=='ready_for_inference'}
    if len(ready)!=ir['ready_markers'] or len(ready)!=br['markers']:raise ValueError('Wrong marker universe')
    source={};expected_by_marker=Counter();accepted_by_marker=Counter()
    for file,accepted in [('paired_site_geometry.tsv',True),('coverage_exclusions.tsv',False)]:
        for r in rows(a.geometry/file):
            key=r['marker'],r['taxon_a'],r['taxon_b']
            if key in source or key[0] not in ready:raise ValueError('Unexpected or repeated geometry pair')
            source[key]=(fingerprint(r),accepted);expected_by_marker[key[0]]+=1;accepted_by_marker[key[0]]+=accepted
    labels=['aa','3di_af','3di_af_empirical','3di_llm'];paths={};tips={}
    for marker in sorted(ready):
        taxa=sorted(r.id for r in SeqIO.parse(a.inputs/marker/'aa.faa','fasta'));tips[marker]=len(taxa)
        if len(taxa)!=len(set(taxa)):raise ValueError('Repeated aligned taxon')
        expected={(marker,x,y) for x,y in combinations(taxa,2)}
        if len(expected)!=expected_by_marker[marker] or any(k not in source for k in expected):raise ValueError('Incomplete source pair grid')
        paths[marker]={}
        for label in labels:
            tree=a.fits/marker/(label+'.treefile');rp=a.fits/marker/(label+'.receipt.json');r=json.loads(rp.read_text());pins[str(rp)]=sha(rp);pins[str(tree)]=r['artifacts'][tree.name]
            if sha(tree)!=pins[str(tree)]:raise ValueError('Changed fitted tree')
            paths[marker][label]=graph_paths(tree,taxa)
    counts=Counter();checked=0;max_error=0.
    for file,accepted in [('path_geometry_points.tsv',True),('geometry_exclusions_with_paths.tsv',False)]:
        for r in rows(a.benchmark/file):
            key=r['marker'],r['taxon_a'],r['taxon_b'];marker,x,y=key
            if key not in source:raise ValueError('Repeated or extra benchmark pair')
            digest,source_accepted=source.pop(key)
            if source_accepted!=accepted:raise ValueError('Acceptance differs')
            if r.pop('marker_review_status')!=review.get(marker,'no_current_copy_review_flag') or r.pop('uncertainty_status')!='point_only_resampling_pending':raise ValueError('Incorrect review or uncertainty label')
            for label in labels:
                reported=float(r.pop(label+'_tree_path_point'));expected=paths[marker][label][x,y]
                if not math.isfinite(reported) or not math.isclose(reported,expected,rel_tol=1e-10,abs_tol=1e-10):raise ValueError('Independently traversed tree path differs')
                max_error=max(max_error,abs(reported-expected));checked+=1
            if fingerprint(r)!=digest:raise ValueError('Inherited geometry field changed')
            counts['accepted' if accepted else 'excluded']+=1
    if source or counts['accepted']!=br['accepted_pairs'] or counts['excluded']!=br['excluded_pairs']:raise ValueError('Incomplete output grid')
    summary={r['marker']:r for r in rows(a.benchmark/'marker_summary.tsv')}
    if set(summary)!=ready:raise ValueError('Summary marker scope differs')
    for marker,r in summary.items():
        values=dict(taxa=tips[marker],pairs=expected_by_marker[marker],accepted_geometry_pairs=accepted_by_marker[marker],excluded_geometry_pairs=expected_by_marker[marker]-accepted_by_marker[marker])
        if any(int(r[k])!=v for k,v in values.items()) or r['marker_review_status']!=review.get(marker,'no_current_copy_review_flag'):raise ValueError('Marker summary differs')
    verify();a.output.mkdir(parents=True)
    result=dict(status='passed_all_tree_path_geometry_point_rows',markers=len(ready),counts=dict(counts),tree_paths_checked=checked,maximum_absolute_path_difference=max_error,source_hashes=pins,script_sha256=sha(__file__),scope='All accepted/excluded pair identities, inherited geometry strings, review/uncertainty labels, summaries and four path values independently checked using adjacency-graph traversal. Same fitted trees and BioPython parser; no independent inference, numerical geometry recomputation, rank-summary audit or biological validation.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
