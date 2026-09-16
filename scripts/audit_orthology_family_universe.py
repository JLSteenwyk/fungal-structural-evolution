#!/usr/bin/env python3
"""Check every family FASTA identity set against the saved clustering output."""
import argparse,hashlib,json,re,time
from collections import Counter
from pathlib import Path


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def groups(path):
    begun=False;number=None;genes=[]
    with path.open() as f:
        for line in f:
            if not begun:
                if 'begin' in line:begun=True
                continue
            if ')' in line:break
            parts=line.strip().split()
            if not parts:continue
            if not line[0].isspace():
                if number is not None:yield number,genes
                number=int(parts.pop(0));genes=[]
            genes.extend(x for x in parts if x!='$')
        if number is not None:yield number,genes


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--working',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    source=a.working/'clusters_OrthoFinder.txt_id_pairs.txt';seen=set();families=set();counts=Counter();taxa=set();started=time.time()
    for number,genes in groups(source):
        family=f'OG{number:07d}'
        if family in families or not genes or len(genes)!=len(set(genes)):raise ValueError('Duplicate/empty cluster')
        if any(not re.fullmatch(r'\d+_\d+',g) for g in genes):raise ValueError('Unexpected gene identity or profile token')
        ids=set(genes)
        if seen.intersection(ids):raise ValueError('Gene assigned to multiple clusters')
        seen.update(ids);families.add(family);taxa.update(g.split('_')[0] for g in genes)
        fa=a.working/'Sequences_ids'/(family+'.fa');headers=[]
        with fa.open() as f:
            for line in f:
                if line.startswith('>'):headers.append(line[1:].strip())
        if len(headers)!=len(set(headers)) or set(headers)!=ids:raise ValueError('FASTA/cluster mismatch: '+family)
        counts['genes']+=len(ids);counts['families']+=1
        counts['singleton_families']+=len(ids)==1;counts['two_sequence_families']+=len(ids)==2
        counts['alignment_eligible_at_least_two']+=len(ids)>=2;counts['tree_eligible_at_least_three']+=len(ids)>=3
        if counts['families']%10000==0:print('Verified',counts['families'],'family identity sets',flush=True)
    if families!={p.stem for p in (a.working/'Sequences_ids').glob('*.fa')}:raise ValueError('Unaccounted family FASTAs')
    species_file=a.working/'SpeciesIDs.txt';species={line.partition(':')[0].strip() for line in species_file.read_text().splitlines() if line.strip() and not line.startswith('#')}
    if not taxa<=species:raise ValueError('Unknown species ID in clustering')
    result={'status':'passed_complete_cluster_to_family_identity_readback','counts':dict(counts),'species_in_clusters':len(taxa),'species_in_mapping':len(species),'species_without_clustered_genes':sorted(species-taxa),'elapsed_seconds':time.time()-started,'cluster_sha256':sha(source),'species_ids_sha256':sha(species_file),'script_sha256':sha(Path(__file__)),'interpretation':'Exact gene identity universe and unique partition, not orthology/homology validation. Native >=2 alignment and >=3 tree eligibility assumes no restart-index exclusions. Source sequence contents are audited separately.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
