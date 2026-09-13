#!/usr/bin/env python3
"""Validate genome-derived complete-codon spans with explicit partial boundaries."""
import argparse
import csv
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from audit_busco_gene_copies import ROOT,sha,read_table
from assess_pae_sensitivity import checked_receipt
from map_proteins_to_genes import attributes


def project(parts,genome,published):
    if not parts:raise ValueError('no_cds_annotation')
    if len({(p['sequence'],p['strand']) for p in parts})!=1:raise ValueError('mixed_sequence_or_strand')
    name=parts[0]['sequence'];strand=parts[0]['strand']
    if name not in genome or strand not in ['+','-']:raise ValueError('unknown_sequence_or_strand')
    ordered=sorted(parts,key=lambda p:p['start'],reverse=strand=='-')
    intervals=sorted((p['start'],p['end']) for p in parts)
    if any(a[1]>=b[0] for a,b in zip(intervals,intervals[1:])):raise ValueError('overlapping_cds_intervals')
    if ordered[0]['phase'] not in ['0','1','2']:raise ValueError('missing_initial_phase')
    initial=int(ordered[0]['phase']);length=0;pieces=[]
    for index,part in enumerate(ordered):
        if not 1<=part['start']<=part['end']<=len(genome[name]):raise ValueError('out_of_bounds')
        expected=initial if index==0 else (3-(length-initial)%3)%3
        if part['phase']!=str(expected):raise ValueError('inconsistent_internal_phase')
        segment=genome[name][part['start']-1:part['end']]
        if strand=='-':segment=str(Seq(segment).reverse_complement())
        pieces.append(segment);length+=len(segment)
    raw=''.join(pieces).upper()
    if raw!=published:raise ValueError('genome_annotation_differs_from_published_cds')
    if len(raw)<=initial:raise ValueError('no_complete_codon')
    tail=(len(raw)-initial)%3;end=len(raw)-tail
    codons=raw[initial:end]
    if not codons:raise ValueError('no_complete_codon')
    return codons,initial,tail,ordered


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--strict',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a new immutable codon projection')
    strict=checked_receipt(args.strict)
    if strict['status']!='complete_published_outgroup_cds_translation_audit':raise ValueError('Completed strict baseline required')
    source_path=ROOT/'metadata/external_genome_receipts.json';sources=json.loads(source_path.read_text())
    if sha(source_path)!=strict['external_source_receipts_sha256']:raise ValueError('Changed source inventory')
    manifest_path=ROOT/'metadata/analysis_manifest.tsv';manifest={r['taxon_id']:r for r in read_table(manifest_path)}
    inputs_path=ROOT/'metadata/qc_input_receipts.json';inputs={r['taxon_id']:r for r in json.loads(inputs_path.read_text())}
    if sha(inputs_path)!=strict['qc_input_receipts_sha256']:raise ValueError('Changed protein inventory')
    baseline={(r['taxon_id'],r['protein_id']):r for r in read_table(args.strict/'protein_cds_audit.tsv')}
    rows=[];summary=[];source_records=[]
    args.output.mkdir(parents=True)
    for taxon in strict['taxa']:
        t=taxon['taxon_id'];m=manifest[t]
        selected={}
        for kind,predicate in [('genome',lambda r:r['url']==m['genome_url']),('cds',lambda r:r['url']==m['cds_url']),
                               ('gff',lambda r:str(r['article_id'])==t[3:] and r['name'].endswith('.gff'))]:
            matches=[r for r in sources if predicate(r)]
            if len(matches)!=1:raise ValueError('Ambiguous deposited source')
            r=matches[0];path=ROOT/r['path']
            if sha(path)!=r['sha256']:raise ValueError('Changed deposited source')
            selected[kind]=path;source_records.append(dict(taxon_id=t,kind=kind,**r))
        dictionaries={}
        for kind in ['genome','cds']:
            records=list(SeqIO.parse(selected[kind],'fasta'));d={r.id:str(r.seq).upper() for r in records}
            if len(d)!=len(records):raise ValueError('Repeated source sequence identity')
            dictionaries[kind]=d
        protein_path=ROOT/inputs[t]['input_path']
        if sha(protein_path)!=taxon['normalized_protein_sha256']:raise ValueError('Changed protein source')
        proteins={r.id:str(r.seq) for r in SeqIO.parse(protein_path,'fasta')}
        parts=defaultdict(list)
        with selected['gff'].open() as handle:
            for number,line in enumerate(handle,1):
                if line.startswith('##FASTA'):break
                if not line.strip() or line.startswith('#'):continue
                f=line.rstrip('\n').split('\t')
                if len(f)!=9:raise ValueError('Invalid GFF row')
                if f[2]!='CDS':continue
                parents=attributes(f[8]).get('Parent',[])
                for parent in parents:parts[parent].append({'sequence':f[0],'start':int(f[3]),'end':int(f[4]),'strand':f[6],'phase':f[7],'gff_line':number})
        local=[];verified={}
        for pid,protein in proteins.items():
            row={'taxon_id':t,'protein_id':pid,'strict_status':baseline[t,pid]['status'],'projection_status':'',
                'initial_partial_bases':'','terminal_partial_bases':'','complete_codon_span_length':'','terminal_stop_in_span':'',
                'omitted_initial_bases':'','omitted_terminal_bases':'','ordered_segments_json':''}
            try:
                published=dictionaries['cds'][pid]
                codons,initial,tail,ordered=project(parts[pid],dictionaries['genome'],published)
                row.update(initial_partial_bases=initial,terminal_partial_bases=tail,complete_codon_span_length=len(codons),
                    omitted_initial_bases=published[:initial],omitted_terminal_bases=published[len(published)-tail:] if tail else '',ordered_segments_json=json.dumps(ordered))
                aa=str(Seq(codons).translate(table=1));stop=aa.endswith('*');row['terminal_stop_in_span']=stop
                if stop:aa=aa[:-1]
                if aa!=protein:raise ValueError('projected_translation_mismatch')
                row['projection_status']='exact_genome_linked_codon_translation';verified[pid]=codons
            except (ValueError,KeyError) as error:row['projection_status']=str(error)
            rows.append(row);local.append(row)
        output=args.output/(t+'.verified_codon_spans.fna')
        with output.open('w') as handle:
            for pid,dna in sorted(verified.items()):handle.write('>'+pid+'\n'+dna+'\n')
        summary.append({'taxon_id':t,'proteins_screened':len(proteins),'verified_codon_spans':len(verified),
            'status_counts':dict(Counter(r['projection_status'] for r in local)),
            'verified_with_partial_boundaries':sum(r['projection_status']=='exact_genome_linked_codon_translation' and (r['initial_partial_bases']>0 or r['terminal_partial_bases']>0) for r in local),
            'strict_exceptions_with_verified_projection':sum(r['projection_status']=='exact_genome_linked_codon_translation' and r['strict_status']!='exact_translation' for r in local)})
    table=args.output/'codon_projection_audit.tsv'
    with table.open('w') as handle:
        writer=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    result={'status':'complete_published_outgroup_genome_codon_projection_audit','taxa':summary,'proteins_screened':len(rows),
        'status_counts':dict(Counter(r['projection_status'] for r in rows)), 'strict_receipt_sha256':sha(args.strict/'receipt.json'),
        'manifest_sha256':sha(manifest_path),'sources':source_records,'script_sha256':sha(Path(__file__)),
        'interpretation':'Genome-derived ordered CDS blocks must reproduce the complete published CDS string. Initial offset comes only from annotated phase; internal phases must agree. Only the annotated initial partial codon and a terminal 1–2-base remainder are omitted, with bases recorded. Every complete translated protein residue must match; no frame or code search. Projected spans are not claims of complete genes, corrected annotations or selection eligibility. Strict unmodified-CDS results remain preserved; boundary-aware spans require separate downstream inclusion policy. Terminal stop codons remain flagged in DNA.',
        'artifacts':{p.name:sha(p) for p in args.output.iterdir()}}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='sources'},indent=2))


if __name__=='__main__':main()
