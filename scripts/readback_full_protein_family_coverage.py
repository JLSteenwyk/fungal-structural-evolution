#!/usr/bin/env python3
"""Read back the full protein partition using OrthoFinder's native cluster reader."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from orthofinder.tools import mcl


def sha(p):
    with p.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    plan=json.loads(args.plan.read_text()); root=Path(plan['output']); wd=Path(plan['working_directory'])
    receipt=json.loads((root/'receipt.json').read_text())
    if receipt['plan_sha256'] != sha(args.plan): raise ValueError('Plan mismatch')
    for name,h in receipt['artifacts'].items():
        if sha(root/name)!=h: raise ValueError('Artifact changed')
    with (root/'taxon_coverage.tsv').open() as handle:
        rows={int(r['species_id']):r for r in csv.DictReader(handle,delimiter='\t')}
    native=mcl.GetPredictedOGs(str(wd/'clusters_OrthoFinder.txt_id_pairs.txt'))
    masks={sid:bytearray(int(r['proteins'])) for sid,r in rows.items()}
    counts={sid:Counter() for sid in rows}; sizes=Counter()
    for family in native:
        category=min(len(family),3); sizes[category]+=1
        for gene in family:
            sid,index=map(int,gene.split('_'))
            if index<0 or masks[sid][index]: raise ValueError('Duplicate native membership')
            masks[sid][index]=category; counts[sid][category]+=1
    del native
    seen={sid:bytearray(len(mask)) for sid,mask in masks.items()}
    missing_count=0
    with (root/'outside_retained_partition.tsv').open() as handle, (wd/'SequenceIDs.txt').open() as source:
        missing=iter(csv.DictReader(handle,delimiter='\t'))
        for line in source:
            gene,label=line.rstrip('\n').split(': ',1); sid,index=map(int,gene.split('_'))
            if masks[sid][index]: continue
            row=next(missing,None)
            if row != {'native_id':gene,'taxon_id':rows[sid]['taxon_id'],'protein_id':label}:
                raise ValueError('Outside-partition identity export differs from native complement')
            if seen[sid][index]: raise ValueError('Repeated outside ID')
            seen[sid][index]=1; missing_count+=1
        if next(missing,None) is not None: raise ValueError('Extra outside IDs')
    labels={1:'singleton_family_proteins',2:'two_member_family_proteins',3:'tree_eligible_family_proteins'}
    for sid,row in rows.items():
        if any(counts[sid][k]!=int(row[label]) for k,label in labels.items()): raise ValueError('Taxon category mismatch')
        if sum(seen[sid])!=int(row['outside_retained_partition']) or sum(counts[sid].values())+sum(seen[sid])!=len(masks[sid]):
            raise ValueError('Taxon partition does not close')
        if abs(float(row['outside_fraction'])-sum(seen[sid])/len(masks[sid]))>1e-14: raise ValueError('Fraction mismatch')
    if dict(sizes)!={int(k):v for k,v in receipt['family_counts_by_size_class'].items()} or missing_count!=receipt['totals']['outside_retained_partition']:
        raise ValueError('Receipt totals mismatch')
    result=dict(status='passed_native_reader_full_partition_readback',taxa=len(rows),families=sum(sizes.values()),
                proteins=sum(map(len,masks.values())),outside_proteins=missing_count,
                receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),native_reader_sha256=sha(Path(mcl.__file__)),
                scope='Native cluster reader independently reconstructs all size classes and the exact complement. Every exported outside ID, source protein label and taxon matches the sequence mapping; all taxa close. Historical intermediate-unassigned comparisons and raw FASTA headers are checked by the producer, not independently reparsed here.')
    args.output.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))


if __name__=='__main__':main()
