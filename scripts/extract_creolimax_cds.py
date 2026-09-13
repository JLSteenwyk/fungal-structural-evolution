#!/usr/bin/env python3
"""Extract Creolimax CDS with exact identifier, strand, phase and translation checks."""
import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from audit_busco_gene_copies import ROOT, sha


def assemble(parts, genome):
    if len({(p['sequence'],p['strand']) for p in parts}) != 1:
        raise ValueError('mixed_sequence_or_strand')
    strand=parts[0]['strand'];name=parts[0]['sequence']
    if strand not in ['+','-'] or name not in genome:raise ValueError('unknown_strand_or_sequence')
    parts=sorted(parts,key=lambda p:p['start'],reverse=strand=='-')
    intervals=sorted((p['start'],p['end']) for p in parts)
    if any(a[1]>=b[0] for a,b in zip(intervals,intervals[1:])):raise ValueError('overlapping_or_duplicate_cds_intervals')
    if parts[0]['phase'] != '0':raise ValueError('nonzero_or_missing_initial_phase')
    pieces=[];length=0
    for part in parts:
        if not 1<=part['start']<=part['end']<=len(genome[name]):raise ValueError('coordinates_out_of_bounds')
        if part['phase'] not in ['0','1','2'] or int(part['phase']) != (3-length%3)%3:
            raise ValueError('inconsistent_internal_phase')
        piece=genome[name][part['start']-1:part['end']]
        if strand=='-':piece=str(Seq(piece).reverse_complement())
        pieces.append(piece);length+=len(piece)
    cds=''.join(pieces).upper()
    if not set(cds)<=set('ACGTRYSWKMBDHVN'):raise ValueError('non_iupac_dna')
    if len(cds)%3:raise ValueError('non_triplet_cds_length')
    return cds,parts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a new immutable extraction audit')
    source_path=ROOT/'metadata/external_genome_receipts.json'
    sources={r['name']:r for r in json.loads(source_path.read_text()) if r.get('article_id')==1403592}
    needed=['Creolimax_fragrantissima.gtf.gz','Creolimax_fragrantissima.genome.fasta.gz','Creolimax_fragrantissima.pep.fasta.gz']
    for name in needed:
        if sha(ROOT/sources[name]['path'])!=sources[name]['sha256']:raise ValueError('Changed deposited source')
    inputs_path=ROOT/'metadata/qc_input_receipts.json'
    input_row=next(r for r in json.loads(inputs_path.read_text()) if r['taxon_id']=='OFS1403592')
    protein_path=ROOT/input_row['input_path']
    if sha(protein_path)!=input_row['sha256'] or input_row['source_sha256']!=sources[needed[2]]['sha256']:raise ValueError('Normalized protein provenance differs')
    with gzip.open(ROOT/sources[needed[1]]['path'],'rt') as handle:
        records=list(SeqIO.parse(handle,'fasta'))
    genome={r.id:str(r.seq) for r in records}
    if len(genome)!=len(records):raise ValueError('Duplicate assembly sequence IDs')
    proteins=list(SeqIO.parse(protein_path,'fasta'))
    if len({r.id for r in proteins})!=len(proteins):raise ValueError('Duplicate protein IDs')
    transcripts=defaultdict(list);genes=defaultdict(set)
    with gzip.open(ROOT/sources[needed[0]]['path'],'rt') as handle:
        for line_number,line in enumerate(handle,1):
            if not line.strip() or line.startswith('#'):continue
            f=line.rstrip('\n').split('\t')
            if len(f)!=9:raise ValueError('Invalid GTF schema')
            if f[2]!='CDS':continue
            attrs=dict(re.findall(r'(\w+) "([^\"]*)"',f[8]))
            tid,gid=attrs['transcript_id'],attrs['gene_id']
            transcripts[tid].append({'sequence':f[0],'start':int(f[3]),'end':int(f[4]),'strand':f[6],'phase':f[7],'gene_id':gid,'gtf_line':line_number})
            genes[gid].add(tid)
    rows=[];verified={};used=set()
    for protein in proteins:
        pid=protein.id
        candidates=[pid] if pid in transcripts else sorted(genes.get(pid,[]))
        row={'protein_id':pid,'protein_length':len(protein.seq),'mapping_basis':'exact_transcript_id' if pid in transcripts else 'exact_gene_id',
            'candidate_transcripts_json':json.dumps(candidates),'status':'','transcript_id':'','cds_length':'','terminal_stop_in_cds':'','ordered_segments_json':''}
        if len(candidates)!=1:
            row['status']='missing_transcript' if not candidates else 'ambiguous_gene_transcripts'
        else:
            tid=candidates[0];used.add(tid);row['transcript_id']=tid
            try:
                if len({p['gene_id'] for p in transcripts[tid]})!=1:raise ValueError('multiple_gene_ids_for_transcript')
                cds,ordered=assemble(transcripts[tid],genome)
                row['cds_length']=len(cds);row['ordered_segments_json']=json.dumps(ordered)
                aa=str(Seq(cds).translate(table=1))
                row['terminal_stop_in_cds']=aa.endswith('*')
                if aa.endswith('*'):aa=aa[:-1]
                if aa!=str(protein.seq):raise ValueError('translation_mismatch')
                row['status']='exact_translation';verified[pid]=cds
            except ValueError as error:row['status']=str(error)
        rows.append(row)
    args.output.mkdir(parents=True)
    fasta=args.output/'verified_cds.fna'
    with fasta.open('w') as handle:
        for pid,sequence in sorted(verified.items()):handle.write('>'+pid+'\n'+sequence+'\n')
    table=args.output/'protein_cds_audit.tsv'
    with table.open('w') as handle:
        writer=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    unused=args.output/'annotation_transcripts_without_selected_protein.json'
    unused.write_text(json.dumps(sorted(set(transcripts)-used),indent=2)+'\n')
    result={'status':'complete_creolimax_cds_extraction_audit','taxon_id':'OFS1403592','proteins_screened':len(proteins),
        'annotation_cds_transcripts':len(transcripts),'status_counts':dict(Counter(r['status'] for r in rows)),
        'mapping_basis_counts':dict(Counter(r['mapping_basis'] for r in rows)),
        'unused_annotation_transcripts':len(set(transcripts)-used),'verified_cds_records':len(verified),
        'translation_table_tested':1,'source_receipt_sha256':sha(source_path),'qc_input_receipts_sha256':sha(inputs_path),
        'normalized_protein_sha256':sha(protein_path),'sources':{name:sources[name] for name in needed},'script_sha256':sha(Path(__file__)),
        'interpretation':'CDS blocks ordered in transcription direction, reverse-complemented on minus strand, bounds/non-overlap and phase continuity checked. Internal phase bases are retained across splice junctions, not individually trimmed. Requires initial phase zero and a multiple-of-three CDS. Standard-code translation must match the complete normalized protein after removing at most one terminal stop from the translation; DNA terminal stops are retained and flagged. Exceptions remain excluded from verified CDS, not repaired. Exact gene-ID fallback only when one transcript is explicitly linked. This does not establish codon alignment or selection-test suitability.',
        'artifacts':{p.name:sha(p) for p in args.output.iterdir()}}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='sources'},indent=2))


if __name__=='__main__':main()
