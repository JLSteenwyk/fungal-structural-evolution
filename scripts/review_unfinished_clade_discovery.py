#!/usr/bin/env python3
"""Reproduce discovery-clade dimensions and singleton-overlap checks."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from orthofinder.orthogroups import accelerate
from orthofinder.tools import mcl


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    review=json.loads(args.review.read_text())
    for name,expected in review['pins'].items():
        if sha(Path(name))!=expected: raise ValueError('Changed review input '+name)
    root=Path('results/orthology/full-protein-family-coverage-v1')
    producer=json.loads((root/'receipt.json').read_text())
    if sha(root/'taxon_coverage.tsv')!=producer['artifacts']['taxon_coverage.tsv']: raise ValueError('Changed taxon table')
    with (root/'taxon_coverage.tsv').open() as handle:
        rows={int(r['species_id']):r for r in csv.DictReader(handle,delimiter='\t')}
    core_file=Path('results/orthology/core-v1/Results_Sep12/WorkingDirectory/SpeciesIDs.txt')
    core={int(line.split(':',1)[0]) for line in core_file.read_text().splitlines() if line.strip()}
    for expected in review['guides']:
        tree=Path('results/orthology/reconciliation-guide-inputs-v1')/expected['guide']/'species_tree_ids.nwk'
        if sha(tree)!=expected['guide_sha256']: raise ValueError('Changed guide')
        clades=accelerate.get_new_species_clades(str(tree),core)
        if clades!=[c['species'] for c in expected['chunks']]: raise ValueError('Native clades differ')
        seen=set(); pairs=0; proteins=0
        for clade,record in zip(clades,expected['chunks']):
            if seen.intersection(clade): raise ValueError('Overlapping clades')
            seen.update(clade)
            sizes=[int(rows[s]['intermediate_unassigned_proteins']) for s in clade]
            ordered=sum(n*m for n in sizes for m in sizes)
            cross=sum(n*m for i,n in enumerate(sizes) for m in sizes[i+1:])
            if ordered!=record['ordered_candidate_pairs_including_self'] or cross!=record['unordered_cross_species_candidate_pairs'] or sum(sizes)!=record['proteins']:
                raise ValueError('Clade workload differs')
            pairs+=ordered; proteins+=sum(sizes)
        if pairs!=expected['ordered_candidate_pairs_including_self'] or proteins!=expected['total_intermediate_proteins'] or set(rows)-seen!=set(expected['uncovered_taxa']):
            raise ValueError('Clade summary differs')
    wd=Path('results/orthology/core-v1/Results_full526v2/WorkingDirectory')
    groups=mcl.GetPredictedOGs(str(wd/'clusters_OrthoFinder.txt_id_pairs.txt'))
    singleton={next(iter(group)) for group in groups if len(group)==1}
    retained=set.union(*groups); intermediate=set()
    hashes=json.loads((root/'intermediate_source_checksums.json').read_text())
    for name,h in hashes.items():
        path=Path(name)
        if sha(path)!=h: raise ValueError('Intermediate changed')
        with path.open() as handle:
            intermediate.update(line[1:].strip() for line in handle if line.startswith('>'))
    if intermediate.intersection(retained)!=singleton or len(singleton)!=review['retained_intermediate_overlap'] or len(intermediate)!=review['intermediate_proteins']:
        raise ValueError('Singleton-overlap identity differs')
    result=dict(status='passed_native_clade_and_singleton_overlap_readback',guides=len(review['guides']),
                clades_per_guide=[r['clades'] for r in review['guides']],intermediate_proteins=len(intermediate),
                retained_singleton_overlap=len(singleton),review_sha256=sha(args.review),script_sha256=sha(Path(__file__)),
                scope='Re-executes installed native clade selection, independently sums all clade species-pair products, and verifies exact intermediate/retained overlap. Does not validate scientific correctness of native clade restriction or perform discovery.')
    args.output.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))


if __name__=='__main__': main()
