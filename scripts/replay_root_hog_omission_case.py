#!/usr/bin/env python3
"""Replay the smallest unflagged omission family from its resolved production tree."""
import csv
import hashlib
import json
from pathlib import Path
from orthofinder.gene_tree_inference import trees2ologs_of as native


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    plan=json.loads(Path('metadata/profile_root_hog_disposition_plan.json').read_text())
    result=Path(plan['result']);source=Path(plan['source'])
    missing=Path(plan['output'])/'missing_genes.tsv'
    with missing.open() as f:
        candidates=[r for r in csv.DictReader(f,delimiter='\t') if r['native_flagged_misplaced']=='False']
    family=min(candidates,key=lambda r:(int(r['family_size']),r['family']))['family']
    trees=result/'Resolved_Gene_Trees/Resolved_Gene_Trees.txt'
    with trees.open() as f:
        newick=next(line.split(': ',1)[1] for line in f if line.startswith(family+': '))
    species={}
    for line in (source/'SpeciesIDs.txt').read_text().splitlines():
        sid,name=line.split(': ',1);species[sid]=name.rsplit('.',1)[0]
    reverse={v:k for k,v in species.items()}
    species_path=result/'Species_Tree/SpeciesTree_rooted_node_labels.txt'
    st=native.tree_lib.Tree(str(species_path),format=1)
    for leaf in st:leaf.name=reverse[leaf.name]
    tree=native.tree_lib.Tree(newick,format=1);labels={}
    for i,leaf in enumerate(tree):
        taxon,protein=leaf.name.split('_',1);key=reverse[taxon]+'_'+str(i)
        labels[key]=leaf.name;leaf.name=key
    _,tree,suspect,dups=native.GetOrthologues_from_tree(int(family[2:]),tree,st,native.GeneToSpecies_dash,native.GetSpeciesNeighbours(st),q_get_dups=True,qNoRecon=True)
    writer=native.HogWriter(st,[n.name for n in st.traverse() if not n.is_leaf()],labels,species,list(map(int,species)),False,write_output=False)
    writer.mark_dups_below(tree)
    replay={}
    for node in tree.traverse('preorder'):
        for level,row in writer.write_clade_v2(node,family):
            if level=='N0':replay[row[1]]={g for cell in row[2:] for g in cell.split(', ') if g}
    hogs=result/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv';csv.field_size_limit(16*1024*1024);expected={}
    with hogs.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if row['OG']!=family:continue
            expected[row['Gene Tree Parent Clade']]={t+'_'+p for t in species.values() for p in row[t].split(', ') if p}
    assert replay==expected and replay
    omitted=set(labels.values())-set().union(*replay.values())
    expected_missing={r['taxon_id']+'_'+r['protein_id'] for r in candidates if r['family']==family}
    assert omitted==expected_missing
    paths=[missing,trees,species_path,source/'SpeciesIDs.txt',hogs,Path(native.__file__),Path(__file__)]
    receipt={'status':'passed_native_resolved_tree_case_replay','family':family,'genes':len(labels),'omitted':sorted(omitted),'native_suspect_genes':sorted(labels[g] for g in suspect),'root_parent_clades':sorted(replay),'replayed_duplications':[{'species_node':d[0],'gene_node':d[1],'support':d[2]} for d in dups], 'pins':{str(p):sha(p) for p in paths},'scope':'Smallest affected family selected deterministically; native reconciliation classification rerun on its already resolved tree, followed by exact root parent/membership comparison. Does not rerun initial rooting/tree resolution, independently validate biological events or establish all-family causality.'}
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
