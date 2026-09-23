#!/usr/bin/env python3
"""Account for root HOG genes using source membership and native flagged-gene lists."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from assess_small_family_output_exposure import groups,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan_sha=sha(args.plan)
    def verify():
        if sha(args.plan)!=plan_sha:raise ValueError('Plan changed')
        for path,digest in plan['pins'].items():
            if sha(path)!=digest:raise ValueError('Input changed: '+path)
    verify();source=Path(plan['source']);result=Path(plan['result']);out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    species={}
    for line in (source/'SpeciesIDs.txt').read_text().splitlines():
        native,name=line.split(': ',1);species[native]=name.rsplit('.',1)[0]
    flagged=set();flag_files=list((result/'Phylogenetically_Misplaced_Genes').glob('*.txt'))
    for path in flag_files:
        if path.stem not in species.values():raise ValueError('Unknown flagged taxon')
        for name in path.read_text().splitlines():
            key=(path.stem,name.strip())
            if not key[1] or key in flagged:raise ValueError('Empty/duplicate native flagged gene')
            flagged.add(key)
    gene_family={};sizes={}
    for family,members in groups(source/'clusters_OrthoFinder.txt_id_pairs.txt'):
        sizes[family]=len(members)
        for gene in members:
            if gene in gene_family:raise ValueError('Repeated source family member')
            gene_family[gene]=family
    source_count=len(gene_family);lookup={t:{} for t in species.values()};singleton_count=0;known_flags=set();flag_rows=[]
    with (source/'SequenceIDs.txt').open() as handle:
        for line in handle:
            native,protein=line.rstrip('\n').split(': ',1);family=gene_family.pop(native);taxon=species[native.split('_')[0]]
            key=(taxon,protein)
            if key in flagged:
                known_flags.add(key);flag_rows.append(dict(taxon_id=taxon,protein_id=protein,native_gene_id=native,family=family,family_size=sizes[family]))
            if sizes[family]==1:singleton_count+=1;continue
            if protein in lookup[taxon]:raise ValueError('Duplicate source protein label')
            lookup[taxon][protein]=(native,family)
    if gene_family or known_flags!=flagged:raise ValueError('Incomplete source identities or unknown flagged genes')
    del gene_family
    csv.field_size_limit(16*1024*1024);observed=0;flagged_present=set();hogs=set()
    with (result/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv').open() as handle:
        table=csv.reader(handle,delimiter='\t');header=next(table)
        if header[:3]!=['HOG','OG','Gene Tree Parent Clade'] or set(header[3:])!=set(lookup):raise ValueError('Wrong root HOG header')
        for row in table:
            if len(row)!=len(header) or row[0] in hogs:raise ValueError('Malformed or duplicate root HOG')
            hogs.add(row[0])
            for taxon,cell in zip(header[3:],row[3:]):
                if not cell:continue
                for protein in cell.split(', '):
                    entry=lookup[taxon].pop(protein,None)
                    if entry is None or entry[1]!=row[1]:raise ValueError('Unknown, singleton, duplicate or wrong-family root assignment')
                    observed+=1
                    if (taxon,protein) in flagged:flagged_present.add((taxon,protein))
    missing={(taxon,protein) for taxon,proteins in lookup.items() for protein in proteins}
    out.mkdir(parents=True);missing_rows=[];counts=Counter()
    for taxon,protein in sorted(missing):
        native,family=lookup[taxon][protein];is_flagged=(taxon,protein) in flagged
        missing_rows.append(dict(taxon_id=taxon,protein_id=protein,native_gene_id=native,family=family,family_size=sizes[family],native_flagged_misplaced=is_flagged));counts[taxon]+=1
    taxon_rows=[dict(taxon_id=t,root_missing_nonsingleton_genes=counts[t]) for t in sorted(lookup)]
    artifacts={}
    for name,rows,fields in [('missing_genes.tsv',missing_rows,['taxon_id','protein_id','native_gene_id','family','family_size','native_flagged_misplaced']),('taxon_disposition.tsv',taxon_rows,['taxon_id','root_missing_nonsingleton_genes'])]:
        with (out/name).open('w') as handle:
            w=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
        artifacts[name]=sha(out/name)
    if observed+len(missing)+singleton_count!=source_count:raise ValueError('Disposition totals inconsistent')
    verify()
    receipt=dict(status='complete_root_hog_gene_disposition',source_genes=source_count,source_families=len(sizes),singleton_genes=singleton_count,nonsingleton_genes=source_count-singleton_count,root_hogs=len(hogs),root_assigned_genes=observed,missing_nonsingleton_genes=len(missing),native_flagged_genes=len(flagged),native_flag_files=len(flag_files),missing_exactly_matches_native_flagged=(missing==flagged),unexplained_missing=len(missing-flagged),native_flagged_but_not_missing=len(flagged-missing),flagged_present_in_root=len(flagged_present),affected_families=len({r['family'] for r in missing_rows}),affected_taxa=len(counts),plan_sha256=plan_sha,artifacts=artifacts,scope='Exact source/root/native-flag disposition. Native misplaced classification is an algorithmic output, not independently established contamination, horizontal transfer or biological error. No validated ancestral membership, duplication event or loss inference.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)


if __name__=='__main__':main()
