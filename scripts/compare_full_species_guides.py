#!/usr/bin/env python3
"""Compare audited unrooted species guides without interpreting support or dates."""
import argparse,csv,json,hashlib
from pathlib import Path
import numpy as np
from Bio import Phylo
from scipy.stats import spearmanr
from audit_species_guide import edges


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def graph_splits(tree,taxa):
    adjacent={n:[] for n in tree.find_clades()}
    for n in adjacent:
        for child in n.clades:adjacent[n].append(child);adjacent[child].append(n)
    found=set()
    for parent in adjacent:
        for child in parent.clades:
            seen={parent};stack=[child];side=set()
            while stack:
                n=stack.pop()
                if n in seen:continue
                seen.add(n)
                if n.is_terminal():side.add(n.name)
                stack.extend(adjacent[n])
            left,right=tuple(sorted(side)),tuple(sorted(taxa-side));found.add(min(left,right,key=lambda x:(len(x),x)))
    return found


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['profile','mafft','manifest','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable comparison')
    with a.manifest.open() as f:manifest={x['taxon_id']:x for x in csv.DictReader(f,delimiter='\t')}
    taxa=sorted(manifest);taxaset=set(taxa);maps={};matrices={};receipts={};source_checks={};traversal_checks=0
    sample_pairs=sorted([(i,j) for i in range(len(taxa)) for j in range(i+1,len(taxa))],key=lambda pair:hashlib.sha256(('|'.join(taxa[i] for i in pair)).encode()).digest())[:100]
    for label in ['profile','mafft']:
        folder=getattr(a,label);r=json.loads((folder/'receipt.json').read_text());receipts[label]=r
        if r['status']!='passed_full_species_guide_readback':raise ValueError('Guide audit missing')
        for name,digest in r['artifacts'].items():
            if sha(folder/name)!=digest:raise ValueError('Changed guide artifact')
        tree=Phylo.read(folder/'guide.treefile','newick');mapping=edges(tree,taxaset)
        if graph_splits(tree,taxaset)!=set(mapping):raise ValueError('Independent graph splits differ')
        maps[label]=mapping;matrix=np.zeros((len(taxa),len(taxa)))
        for split,length in mapping.items():
            inside=np.array([t in split for t in taxa]);matrix+=np.logical_xor(inside[:,None],inside[None,:])*length
        for i,j in sample_pairs:
            if not np.isclose(matrix[i,j],tree.distance(taxa[i],taxa[j]),atol=1e-10,rtol=1e-10):raise ValueError('Path sum/traversal differs')
            traversal_checks+=1
        matrices[label]=matrix;source_checks[label]=sha(folder/'receipt.json')
    internal={label:{s for s in mapping if len(s)>1} for label,mapping in maps.items()};shared=internal['profile']&internal['mafft'];union=internal['profile']|internal['mafft'];rf=len(union-shared)
    splits=[]
    for side in sorted(union,key=lambda s:(len(s),s)):
        splits.append({'split_smaller_side_taxa':';'.join(side),'smaller_side_taxa':len(side),'status':'shared' if side in shared else 'profile_only' if side in internal['profile'] else 'mafft_only','profile_branch_length':maps['profile'].get(side,''),'mafft_branch_length':maps['mafft'].get(side,'')})
    nearest=[]
    for label in matrices:
        for focal in ['F610337','F27376','O400682']:
            i=taxa.index(focal);ranked=sorted((float(matrices[label][i,j]),taxa[j]) for j in range(len(taxa)) if j!=i)[:10]
            for rank,(distance,other) in enumerate(ranked,1):nearest.append({'guide':label,'focal_taxon':focal,'focal_species':manifest[focal]['species_name'],'rank':rank,'neighbor_taxon':other,'neighbor_species':manifest[other]['species_name'],'tree_path_length':distance})
    upper=np.triu_indices(len(taxa),1);corr=float(spearmanr(matrices['profile'][upper],matrices['mafft'][upper]).statistic)
    a.output.mkdir(parents=True)
    for name,data in [('internal_split_comparison.tsv',splits),('flagged_taxon_nearest_guide_paths.tsv',nearest)]:
        with (a.output/name).open('w',newline='') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    result={'status':'complete_audited_full_species_guide_comparison','taxa':len(taxa),'internal_splits_per_guide':{k:len(v) for k,v in internal.items()},'shared_internal_splits':len(shared),'unrooted_robinson_foulds_distance':rf,'normalized_rf_by_total_internal_splits':rf/sum(map(len,internal.values())),'tip_pair_distances_compared':len(upper[0]),'descriptive_tip_pair_path_spearman':corr,'independent_graph_split_checks':sum(map(len,maps.values())),'independent_sampled_path_traversals':traversal_checks,'source_audit_receipts':source_checks,'manifest_sha256':sha(a.manifest),'script_sha256':sha(Path(__file__)),'artifacts':{p.name:sha(p) for p in a.output.iterdir()},'interpretation':'Unrooted split disagreement and descriptive path rank association between guides from different alignments. No support, uncertainty, topology significance, calibrated dates, independent pair observations or likelihood comparison across alignments. Nearest tree paths are descriptive and cannot confirm contamination or taxonomic identity; homogeneous-model composition concerns remain.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
