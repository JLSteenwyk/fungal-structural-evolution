#!/usr/bin/env python3
"""Read back every codon alignment against its original masked amino-acid columns."""
import argparse,csv,json
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from audit_busco_gene_copies import ROOT,sha,read_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--alignments',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable readback directory')
    r=json.loads((a.alignments/'receipt.json').read_text())
    if r['status']!='complete_full_marker_codon_alignment_projection':raise ValueError('Completed full codon projection required')
    for name,h in r['artifacts'].items():
        if sha(a.alignments/name)!=h:raise ValueError('Changed codon artifact')
    rows=read_table(a.alignments/'sequence_audit.tsv');keys={(x['marker'],x['taxon_id']) for x in rows}
    original=read_table(ROOT/'results/cds/marker-source-index-v1/marker_cds_source_index.tsv')
    if len(keys)!=len(rows) or keys!={(x['marker'],x['taxon_id']) for x in original}:raise ValueError('Incomplete original identity coverage')
    accepted={(x['marker'],x['taxon_id']):x for x in rows if x['status']=='translation_verified_codon_alignment'}
    masks=json.loads((a.alignments/'protein_source_columns.json').read_text())
    audit_path=ROOT/'results/phylogeny/mafft-audit-full-v1/receipt.json'
    if sha(audit_path)!=r['protein_alignment_audit_receipt_sha256']:raise ValueError('Changed protein audit')
    protein_hashes={x['marker']:x['alignment_sha256'] for x in json.loads(audit_path.read_text())['alignment_receipts']}
    summaries=[];seen=set();translated_codons=0
    for m in r['markers']:
        marker=m['marker'];pp=ROOT/'results/phylogeny/alignments-full-v1'/(marker+'.faa')
        if sha(pp)!=protein_hashes[marker]:raise ValueError('Changed protein alignment')
        aas={x.id:str(x.seq).upper() for x in SeqIO.parse(pp,'fasta')};columns=masks[marker];codes=Counter();n=0
        for rec in SeqIO.parse(a.alignments/(marker+'.fna'),'fasta'):
            key=(marker,rec.id)
            if key in seen or key not in accepted:raise ValueError('Unexpected or duplicate codon row')
            seen.add(key);row=accepted[key];code=int(row['translation_table']);dna=str(rec.seq)
            if len(dna)!=3*len(columns) or len(dna)!=m['nucleotide_columns']:raise ValueError('Wrong codon alignment length')
            codons=[dna[i:i+3] for i in range(0,len(dna),3)]
            if any('-' in c and c!='---' for c in codons):raise ValueError('Partial codon gap')
            translated=''.join('-' if c=='---' else str(Seq(c).translate(table=code)) for c in codons)
            expected=''.join(aas[rec.id][i-1] for i in columns)
            if translated!=expected:raise ValueError('Codons differ from masked protein columns')
            if '*' in translated:raise ValueError('Internal aligned stop')
            codes[str(code)]+=1;n+=1;translated_codons+=sum(c!='---' for c in codons)
        if n!=m['accepted_taxa'] or dict(codes)!=m['translation_codes']:raise ValueError('Marker summary mismatch')
        summaries.append({'marker':marker,'input_taxa':m['input_taxa'],'accepted_taxa':n,'excluded_taxa':m['input_taxa']-n,'codon_columns':len(columns),'code_count':len(codes),'translation_codes_json':json.dumps(dict(sorted(codes.items())))})
    if seen!=set(accepted):raise ValueError('Missing codon rows')
    a.output.mkdir(parents=True);table=a.output/'marker_codon_summary.tsv'
    with table.open('w') as f:
        w=csv.DictWriter(f,list(summaries[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summaries)
    result={'status':'complete_full_marker_codon_alignment_readback','source_receipt_sha256':sha(a.alignments/'receipt.json'),'marker_links':len(rows),'verified_aligned_sequences':len(seen),'excluded_sequences':len(rows)-len(seen),'verified_non_gap_codons':translated_codons,'markers':len(summaries),'mixed_code_markers':sum(x['code_count']>1 for x in summaries),'artifacts':{table.name:sha(table)},'script_sha256':sha(Path(__file__)),'interpretation':'Every output codon translated under its recorded code and matched the exact original masked amino-acid column; all source identity/status and output artifact checks passed. Not evidence of saturation suitability, independent ecological replication, gene reconciliation or selection.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
