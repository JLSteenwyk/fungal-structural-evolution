#!/usr/bin/env python3
"""Measure uncorrected codon/AA differences in coverage-screened groups."""
import argparse,csv,itertools,json,statistics
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.Data import CodonTable
from audit_busco_gene_copies import ROOT,sha,read_table


def encode(sequence,code):
    if len(sequence)%3:raise ValueError('Non-triplet aligned sequence')
    b=np.frombuffer(sequence.encode('ascii'),dtype=np.uint8).reshape(-1,3)
    lookup=np.full(256,-1,dtype=np.int16)
    for i,x in enumerate(b'ACGT'):lookup[x]=i
    v=lookup[b];valid=(v>=0).all(axis=1);ids=(v[:,0]*16+v[:,1]*4+v[:,2]).clip(0,63)
    table=CodonTable.unambiguous_dna_by_id[code]
    aa=np.array([ord(table.forward_table.get(''.join(c),'*')) for c in itertools.product('ACGT',repeat=3)],dtype=np.uint8)[ids]
    if np.any(valid & (aa==ord('*'))):raise ValueError('Unexpected stop in validated codon alignment')
    return b,valid,ids,aa


def differences(left,right):
    a,av,ac,aa=left;b,bv,bc,ba=right
    if a.shape!=b.shape:raise ValueError('Different alignment lengths')
    shared=av&bv;n=int(shared.sum());codon=int(np.count_nonzero(ac[shared]!=bc[shared]));amino=int(np.count_nonzero(aa[shared]!=ba[shared]))
    nt=(a[shared]!=b[shared]).sum(axis=0)
    return {'shared_called_codons':n,'different_codons':codon,'different_amino_acids':amino,'different_codons_same_amino_acid':codon-amino,'position1_differences':int(nt[0]),'position2_differences':int(nt[1]),'position3_differences':int(nt[2]),'codon_difference_fraction':codon/n if n else '', 'amino_acid_difference_fraction':amino/n if n else '', 'third_position_difference_fraction':int(nt[2])/n if n else '', 'nucleotide_difference_fraction':int(nt.sum())/(3*n) if n else '', 'passes_pair_overlap_screen':n>=100}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable pair assessment')
    group=ROOT/'results/cds/genus-code-coverage-v1';gr=json.loads((group/'receipt.json').read_text())
    if gr['status']!='complete_full_marker_genus_code_coverage_screen':raise ValueError('Completed group coverage required')
    for name,h in gr['artifacts'].items():
        if sha(group/name)!=h:raise ValueError('Changed group artifact')
    base=ROOT/'results/cds/marker-codon-alignments-v1';cr=json.loads((base/'receipt.json').read_text())
    if sha(base/'receipt.json')!=gr['source_codon_receipt_sha256']:raise ValueError('Changed codon provenance')
    source_rows=read_table(base/'sequence_audit.tsv');code_of={(x['marker'],x['taxon_id']):int(x['translation_table']) for x in source_rows if x['status']=='translation_verified_codon_alignment'}
    groups=read_table(group/'group_codon_coverage.tsv');selected=[x for x in groups if x['passes_coverage_screen']=='True'];by_marker=defaultdict(list)
    for x in selected:by_marker[x['marker']].append(x)
    rows=[];summaries=[];unique_pairs=set()
    for marker,local in by_marker.items():
        path=base/(marker+'.fna')
        if sha(path)!=cr['artifacts'][path.name]:raise ValueError('Changed codon alignment')
        seqs={x.id:str(x.seq) for x in SeqIO.parse(path,'fasta')};encoded={};cache={}
        for g in local:
            taxa=sorted(json.loads(g['retained_taxa_json']));code=int(g['translation_table']);localrows=[]
            for t in taxa:
                if code_of[marker,t]!=code:raise ValueError('Group genetic code mismatch')
                if t not in encoded:encoded[t]=encode(seqs[t],code)
            for t,u in itertools.combinations(taxa,2):
                key=(t,u)
                if key not in cache:cache[key]=differences(encoded[t],encoded[u])
                row={'marker':marker,'genus_label':g['genus_label'],'translation_table':code,'policy':g['policy'],'taxon_a':t,'taxon_b':u,**cache[key]}
                rows.append(row);localrows.append(row);unique_pairs.add((marker,code,t,u))
            passing=[x for x in localrows if x['passes_pair_overlap_screen']]
            summary={'marker':marker,'genus_label':g['genus_label'],'translation_table':code,'policy':g['policy'],'taxa':len(taxa),'pairs':len(localrows),'pairs_at_least_100_codons':len(passing),'minimum_shared_codons':min(x['shared_called_codons'] for x in localrows)}
            for key in ['codon_difference_fraction','amino_acid_difference_fraction','third_position_difference_fraction']:
                values=[x[key] for x in passing]
                summary['median_'+key]=statistics.median(values) if values else ''
                summary['maximum_'+key]=max(values) if values else ''
            summaries.append(summary)
        print(marker,len(rows),flush=True)
    a.output.mkdir(parents=True)
    for name,data in [('pair_differences.tsv',rows),('group_pair_summary.tsv',summaries)]:
        with (a.output/name).open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    result={'status':'complete_coverage_group_observed_codon_divergence','group_policy_rows':len(summaries),'pair_policy_rows':len(rows),'unique_marker_code_pairs':len(unique_pairs),'pair_overlap_pass_by_policy':dict(Counter(x['policy'] for x in rows if x['passes_pair_overlap_screen'])),'group_receipt_sha256':sha(group/'receipt.json'),'codon_receipt_sha256':sha(base/'receipt.json'),'script_sha256':sha(Path(__file__)),'interpretation':'Uncorrected observed differences on pairwise shared unambiguous codons in existing masks. Different codons encoding the same amino acid are counted as observations, not reconstructed synonymous substitutions. No dN/dS, dS, saturation test, branch rate or independent-sample inference. Pairs share taxa/sites/ancestry and repeated policy rows are not replicates. Same genus label is not evidence of monophyly. Summaries of divergence fractions use pairs with at least 100 shared called codons; all lower-overlap pairs remain explicit.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))

if __name__=='__main__':main()
