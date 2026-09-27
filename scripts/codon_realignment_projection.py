#!/usr/bin/env python3
"""Project source codons into a local protein alignment and track correspondences."""
from Bio.Data import CodonTable
AA=set('ACDEFGHIKLMNPQRSTVWY')
def project(aligned,dna,code):
    assert aligned and set(aligned)==set(dna) and len({len(s) for s in aligned.values()})==1
    table=CodonTable.unambiguous_dna_by_id[code];positions={};codons={}
    for taxon,protein in aligned.items():
        source=dna[taxon];assert len(source)%3==0;pos=0;indices=[];triplets=[]
        for amino in protein:
            if amino=='-':indices.append(0);triplets.append('???');continue
            codon=source[3*pos:3*pos+3];assert len(codon)==3;pos+=1
            canonical=set(codon)<=set('ACGT')
            if canonical:assert codon in table.forward_table and table.forward_table[codon]==amino
            called=canonical and amino in AA
            indices.append(pos if called else 0);triplets.append(codon if called else '???')
        assert pos*3==len(source);positions[taxon]=indices;codons[taxon]=triplets
    n=len(aligned);width=len(next(iter(aligned.values())))
    retained=[i for i in range(width) if 5*sum(positions[t][i]>0 for t in aligned)>=4*n]
    return positions,codons,retained

def correspondence(left,right):return {(a,b) for a,b in zip(left,right) if a and b}
def compare(old,raw,retained):
    shared=old&retained;raw_shared=old&raw;union=old|retained
    assert retained<=raw
    return dict(original_pairs=len(old),raw_local_pairs=len(raw),retained_local_pairs=len(retained),original_preserved_raw=len(raw_shared),original_preserved_retained=len(shared),original_absent_from_raw=len(old-raw),original_preserved_but_filtered=len(raw_shared-retained),local_retained_not_original=len(retained-old),original_fraction_retained=len(shared)/len(old) if old else '',retained_pair_jaccard=len(shared)/len(union) if union else '')
