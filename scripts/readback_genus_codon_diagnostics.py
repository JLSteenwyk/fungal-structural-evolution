#!/usr/bin/env python3
"""Independently verify the complete exported genus/code diagnostic alignments."""
import argparse
import json
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['assessment', 'codons', 'coverage', 'hybrids', 'output']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new readback receipt')
    r=checked_receipt(a.assessment);checked_receipt(a.codons);checked_receipt(a.coverage)
    for name in ['codons','coverage']:
        if r['source_receipts'][name]!=sha(getattr(a,name)/'receipt.json'):raise ValueError('Changed source')
    if r['hybrid_table_sha256']!=sha(a.hybrids):raise ValueError('Changed hybrid policy')
    hybrids={x['taxon_id'] for x in read_table(a.hybrids)}
    source_groups={(x['genus_label'],x['marker'],x['translation_table']):x for x in read_table(a.coverage/'group_codon_coverage.tsv') if x['policy']=='exclude_recorded_annotation_gene_label_flags'}
    cases=read_table(a.assessment/'case_summary.tsv');seen=set();ready=0;observations=0;sites=0;cache={}
    source_columns=json.loads((a.codons/'protein_source_columns.json').read_text())
    for row in cases:
        marker=row['marker'];key=row['genus_label'],marker,row['translation_table']
        if key in seen:raise ValueError('Duplicate case')
        seen.add(key);group=source_groups[key];taxa=set(json.loads(group['retained_taxa_json']))-hybrids
        if len(taxa)!=int(row['retained_taxa']):raise ValueError('Retained taxon count differs')
        if marker not in cache:cache[marker]={x.id:str(x.seq) for x in SeqIO.parse(a.codons/(marker+'.fna'),'fasta')}
        original={t:cache[marker][t] for t in taxa};width=len(source_columns[marker])
        kept=[]
        for i in range(width):
            called=sum(all(c in 'ACGT' for c in s[3*i:3*i+3]) for s in original.values())
            if taxa and 5*called>=4*len(taxa):kept.append(i)
        called_by_taxon={t:sum(all(c in 'ACGT' for c in original[t][3*i:3*i+3]) for i in kept) for t in taxa}
        eligible=len(taxa)>=4 and len(kept)>=100 and all(called_by_taxon.values())
        if eligible!=(row['status']=='ready_for_tree_and_divergence_diagnostics') or len(kept)!=int(row['retained_codon_columns']):raise ValueError('Eligibility/mask differs')
        folder=a.assessment/row['case_id']
        if not eligible:
            if folder.exists():raise ValueError('Unexpected excluded-case alignment')
            continue
        dna_records=list(SeqIO.parse(folder/'codons.fna','fasta'));aa_records=list(SeqIO.parse(folder/'amino_acids.faa','fasta'))
        dna={x.id:str(x.seq) for x in dna_records};aa={x.id:str(x.seq) for x in aa_records}
        if len(dna)!=len(dna_records) or len(aa)!=len(aa_records) or set(dna)!=taxa or set(aa)!=taxa:raise ValueError('Exported taxa differ')
        columns=read_table(folder/'columns.tsv')
        if len(columns)!=len(kept):raise ValueError('Exported column count differs')
        for j,i in enumerate(kept):
            if columns[j]!={'diagnostic_codon_column_1based':str(j+1),'source_codon_column_1based':str(i+1),'source_mafft_protein_column_1based':str(source_columns[marker][i])}:raise ValueError('Exported source mapping differs')
        for t in taxa:
            if len(dna[t])!=3*len(kept) or len(aa[t])!=len(kept):raise ValueError('Exported dimensions differ')
            for j,i in enumerate(kept):
                old=original[t][3*i:3*i+3];codon=dna[t][3*j:3*j+3]
                expected=old if all(c in 'ACGT' for c in old) else '???'
                if codon!=expected:raise ValueError('Exported codon differs')
                if codon=='???':
                    if aa[t][j]!='?':raise ValueError('Missingness differs')
                else:
                    if str(Seq(codon).translate(table=int(row['translation_table'])))!=aa[t][j] or aa[t][j]=='*':raise ValueError('Translation differs')
                    observations+=1
        ready+=1;sites+=len(kept)
    if seen!=set(source_groups) or ready!=r['ready_cases'] or observations!=r['observed_taxon_codons_in_ready_cases']:raise ValueError('Incomplete full case/codon grid')
    result={'status':'passed_full_genus_codon_diagnostic_readback','cases':len(cases),'ready_cases':ready,'retained_columns_across_ready_cases':sites,'observed_taxon_codons':observations,'source_receipt_sha256':sha(a.assessment/'receipt.json'),'script_sha256':sha(Path(__file__)),'scope':'All case identities, updated taxon membership, integer 80% occupancy masks, emitted DNA/AA, translation codes and source coordinates independently checked. No phylogeny, saturation, orthology or selection validity claim.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
