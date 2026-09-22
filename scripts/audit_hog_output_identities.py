#!/usr/bin/env python3
"""Check all HOG output identities and species-clade containment; not HOG inference."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
from Bio import Phylo
from assess_small_family_output_exposure import groups, sha


def audit(source, result, output):
    if output.exists():raise FileExistsError(output)
    csv.field_size_limit(16 * 1024 * 1024)
    tree_path=result/'Species_Tree/SpeciesTree_rooted_node_labels.txt'
    paths=[source/name for name in ('SpeciesIDs.txt','SequenceIDs.txt','clusters_OrthoFinder.txt_id_pairs.txt')]+[tree_path]
    pins={str(p):sha(p) for p in paths}
    species={}
    for line in (source/'SpeciesIDs.txt').read_text().splitlines():
        native,name=line.split(': ',1)
        if native in species:raise ValueError('Repeated source species')
        species[native]=name.rsplit('.',1)[0]
    if len(set(species.values()))!=len(species):raise ValueError('Ambiguous source species')
    family_by_gene={};family_count=0
    for family,genes in groups(source/'clusters_OrthoFinder.txt_id_pairs.txt'):
        family_count+=1
        for gene in genes:
            if gene in family_by_gene:raise ValueError('Repeated source gene membership')
            family_by_gene[gene]=family
    count=len(family_by_gene);lookup={taxon:{} for taxon in species.values()}
    with (source/'SequenceIDs.txt').open() as handle:
        for line in handle:
            gene,protein=line.rstrip('\n').split(': ',1)
            if len(protein.split())!=1:raise ValueError('Full-header naming needs separate review')
            taxon=species[gene.split('_')[0]]
            if protein in lookup[taxon]:raise ValueError('Ambiguous source protein identity')
            lookup[taxon][protein]=(gene,family_by_gene.pop(gene))
    if family_by_gene:raise ValueError('Source genes missing sequence labels')
    del family_by_gene
    tree=Phylo.read(tree_path,'newick');desc={};tips=[]
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal():tips.append(node.name)
        else:
            if node.name in desc or not node.name:raise ValueError('Ambiguous tree node label')
            desc[node.name]={tip.name for tip in node.get_terminals()}
    if len(tips)!=len(set(tips)) or set(tips)!=set(lookup):raise ValueError('Tree/source taxa differ')
    folder=result/'Phylogenetic_Hierarchical_Orthogroups'
    files={p.stem:p for p in folder.glob('*.tsv')}
    if set(files)!=set(desc):raise ValueError('Missing or extra HOG node tables')
    records=[];output.mkdir(parents=True)
    expected_header=['HOG','OG','Gene Tree Parent Clade']+[species[k] for k in sorted(species,key=int)]
    for node in sorted(files,key=lambda name:int(name[1:])):
        path=files[node];before=path.stat();digest=hashlib.sha256();seen=set();hogs=set();families=set();members=0
        def lines():
            with path.open('rb') as handle:
                for raw in handle:
                    digest.update(raw);yield raw.decode()
        reader=csv.reader(lines(),delimiter='\t')
        if next(reader)!=expected_header:raise ValueError('Unexpected HOG header: '+node)
        for row in reader:
            if len(row)!=len(expected_header):raise ValueError('Malformed HOG row')
            hog,family,parent=row[:3]
            prefix=node+'.HOG'
            if not hog.startswith(prefix) or not hog[len(prefix):].isdigit() or hog in hogs:
                raise ValueError('Invalid or repeated HOG ID')
            hogs.add(hog);families.add(family);n=0
            for taxon,cell in zip(expected_header[3:],row[3:]):
                if not cell:continue
                if taxon not in desc[node]:raise ValueError('Gene outside HOG species clade')
                for protein in cell.split(', '):
                    if protein not in lookup[taxon]:raise ValueError('Unknown HOG protein')
                    gene,expected_family=lookup[taxon][protein]
                    if family!=expected_family:raise ValueError('HOG changes source family')
                    if gene in seen:raise ValueError('Gene repeated within HOG level')
                    seen.add(gene);n+=1
            if not n:raise ValueError('Empty HOG')
            members+=n
        if {int(h[len(node+'.HOG'):]) for h in hogs}!=set(range(len(hogs))):raise ValueError('Nonconsecutive HOG IDs')
        after=path.stat()
        if (before.st_size,before.st_mtime_ns,before.st_ino)!=(after.st_size,after.st_mtime_ns,after.st_ino):raise ValueError('HOG file changed during audit')
        records.append(dict(node=node,descendant_taxa=len(desc[node]),hogs=len(hogs),source_families=len(families),gene_assignments=members,bytes=before.st_size,sha256=digest.hexdigest()))
        print(node,len(hogs),members,flush=True)
    for path,digest in pins.items():
        if sha(path)!=digest:raise ValueError('Source changed during audit')
    table=output/'node_summary.tsv'
    with table.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(records[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(records)
    receipt=dict(status='passed_all_hog_identity_and_species_clade_checks',source_proteins=count,source_families=family_count,
                 node_tables=len(records),hog_rows=sum(r['hogs'] for r in records),gene_assignments_across_levels=sum(r['gene_assignments'] for r in records),
                 input_hashes=pins,script_sha256=sha(__file__),cluster_parser_sha256=sha(Path(__file__).with_name('assess_small_family_output_exposure.py')),
                 artifacts={'node_summary.tsv':sha(table)},scope='Every emitted HOG protein checked against source taxon/protein/family identities, species-tree clade containment and nonduplication within each hierarchical level. Expected internal-node table set checked. Repeated genes across levels are expected. Does not prove completeness of gene assignment, ancestral HOG membership, parent gene-tree clades, nested orthology, rooting or duplication events.')
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text())
    for path,digest in plan['pins'].items():
        if sha(path)!=digest:raise ValueError('Changed pinned dependency')
    stage=json.loads(Path(plan['native_stage']).read_text())['stage']
    if stage['status']!='native_execution_complete_pending_full_output_readback' or stage['result']!=plan['result']:
        raise ValueError('Native stage is incomplete or different')
    audit(Path(plan['source']),Path(plan['result']),Path(plan['output']))


if __name__=='__main__':main()
