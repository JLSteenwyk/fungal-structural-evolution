#!/usr/bin/env python3
"""Document an exploratory flagged-marker placement in a completed gene tree."""
import argparse,csv,json,hashlib,math
from pathlib import Path
from Bio import Phylo,SeqIO


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use immutable output')
    marker='129234at2759';focal='F610337';folder=Path('results/phylogeny/marker-gene-trees-v2')/marker;rp=folder/'receipt.json';r=json.loads(rp.read_text());manifest=Path('metadata/analysis_manifest.tsv')
    with manifest.open() as f:taxa={x['taxon_id']:x for x in csv.DictReader(f,delimiter='\t')}
    command=r['command'];alignment=Path(command[command.index('-s')+1]);tp=Path(r['tree_path'])
    if r['status']!='inferred' or r['returncode']!=0 or sha(tp)!=r['tree_sha256'] or sha(alignment)!=r['input_sha256']:raise ValueError('Incomplete or changed gene tree')
    if command[command.index('--alrt')+1]!='1000':raise ValueError('Unexpected support method')
    records=list(SeqIO.parse(alignment,'fasta'));expected={x.id for x in records};tree=Phylo.read(tp,'newick');tips={x.name for x in tree.get_terminals()}
    if tips!=expected or len(tips)!=r['taxa'] or focal not in tips:raise ValueError('Tree/input tips differ')
    if any(len(x.seq)!=r['columns'] for x in records):raise ValueError('Input column count differs')
    for node in tree.find_clades():
        if node is not tree.root and (node.branch_length is None or not math.isfinite(node.branch_length) or node.branch_length<0):raise ValueError('Invalid branch')
        if node.confidence is not None and not 0<=node.confidence<=100:raise ValueError('Invalid support')
    # Find all unrooted split sides with the focal plus >=3 Ascomycota entries.
    candidates=[]
    for node in tree.find_clades():
        if node is tree.root:continue
        side={x.name for x in node.get_terminals()}
        for group in [side,tips-side]:
            if focal in group and len(group)>=4 and all(taxa[x]['lineage'].split(';')[0]=='Ascomycota' for x in group-{focal}):
                candidates.append({'taxa':sorted(group),'sh_alrt_support':node.confidence,'separating_branch_length':node.branch_length})
    if not candidates:raise ValueError('Exploratory source-review split not present')
    selected=max(candidates,key=lambda r:len(r['taxa']));neighbors=sorted((tree.distance(focal,t),t) for t in tips if t!=focal)[:10]
    fcs=Path('results/qc/fcs-cds-overlap-v1/marker_overlap_review.tsv')
    with fcs.open() as f:flag=next(x for x in csv.DictReader(f,delimiter='\t') if x['marker']==marker and x['taxon_id']==focal)
    mapping=json.loads((fcs.parent/'receipt.json').read_text())
    if sha(fcs)!=mapping['artifacts'][fcs.name] or flag['fcs_overlap_status']!='cds_overlaps_fcs_region':raise ValueError('FCS flag lineage differs')
    a.output.mkdir(parents=True)
    rows=[{'taxon_id':t,'species_name':taxa[t]['species_name'],'manifest_lineage':taxa[t]['lineage'],'is_focal':t==focal} for t in selected['taxa']]
    with (a.output/'reported_split_taxa.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result={'status':'complete_exploratory_flagged_marker_tree_source_review','marker':marker,'focal_taxon':focal,'focal_species':taxa[focal]['species_name'],'focal_protein':flag['protein_id'],'fcs_actions':flag['fcs_actions'],'tree_taxa':len(tips),'columns':r['columns'],'largest_ascomycota_plus_focal_split':selected,'nearest_tree_path_taxa':[{'taxon_id':t,'species_name':taxa[t]['species_name'],'path_length':d} for d,t in neighbors],'source_gene_tree_receipt_sha256':sha(rp),'tree_sha256':sha(tp),'alignment_sha256':sha(alignment),'manifest_sha256':sha(manifest),'fcs_overlap_receipt_sha256':sha(fcs.parent/'receipt.json'),'script_sha256':sha(Path(__file__)),'artifacts':{'reported_split_taxa.tsv':sha(a.output/'reported_split_taxa.tsv')},'interpretation':'Exploratory marker-specific discordance: focal manifest Basidiomycota entry lies on an unrooted split with Ascomycota entries in this completed gene tree. Support is reported SH-aLRT, not bootstrap probability or a recomputed test. This is consistent with the FCS concern but does not distinguish mixed-source assembly, misannotation, paralogy, gene history or model error. No genome-wide taxonomic reassignment.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['nearest_tree_path_taxa','artifacts']},indent=2))


if __name__=='__main__':main()
