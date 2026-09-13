#!/usr/bin/env python3
"""Validate complete published CDS/protein pairs for four external outgroups."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from audit_busco_gene_copies import ROOT, sha, read_table


def translate_match(cds, protein):
    if len(cds)%3:return 'non_triplet_cds_length',False
    if not cds or not set(cds)<=set('ACGTRYSWKMBDHVN'):return 'invalid_dna',False
    aa=str(Seq(cds).translate(table=1));stop=aa.endswith('*')
    if stop:aa=aa[:-1]
    return ('exact_translation' if aa==protein else 'translation_mismatch'),stop


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a new immutable translation audit')
    inventory_path=ROOT/'metadata/external_cds_source_inventory.json'
    inventory=json.loads(inventory_path.read_text())
    taxa=[r for r in inventory['taxa'] if r['status']=='available_published_cds_translation_pending']
    expected={'OFS5426470','OFS5426494','OFS5426506','OFS5426458'}
    if {r['taxon_id'] for r in taxa}!=expected or len(taxa)!=4:raise ValueError('Expected four published outgroup CDS sources')
    input_path=ROOT/'metadata/qc_input_receipts.json'
    inputs={r['taxon_id']:r for r in json.loads(input_path.read_text())}
    source_path=ROOT/'metadata/external_genome_receipts.json'
    sources=json.loads(source_path.read_text())
    marker_path=ROOT/'results/phylogeny/markers-full-v1/protein_mapping.tsv'
    markers=read_table(marker_path)
    args.output.mkdir(parents=True)
    all_rows=[];summaries=[]
    for taxon in taxa:
        name=taxon['taxon_id'];source=ROOT/taxon['path'];protein_path=ROOT/inputs[name]['input_path']
        if sha(source)!=taxon['sha256'] or sha(protein_path)!=inputs[name]['sha256']:raise ValueError('Changed CDS or protein source')
        original=[r for r in sources if r['sha256']==inputs[name]['source_sha256']]
        if len(original)!=1 or sha(ROOT/original[0]['path'])!=inputs[name]['source_sha256']:raise ValueError('Deposited protein provenance differs')
        cds_records=list(SeqIO.parse(source,'fasta'));protein_records=list(SeqIO.parse(protein_path,'fasta'))
        cds={r.id:str(r.seq).upper() for r in cds_records};proteins={r.id:str(r.seq) for r in protein_records}
        if len(cds)!=len(cds_records) or len(proteins)!=len(protein_records):raise ValueError('Repeated record ID')
        rows=[];verified={}
        for pid,protein in proteins.items():
            dna=cds.get(pid)
            status,stop=translate_match(dna,protein) if dna is not None else ('no_exact_cds_identifier',False)
            row={'taxon_id':name,'protein_id':pid,'protein_length':len(protein),'cds_length':len(dna) if dna is not None else '',
                'status':status,'terminal_stop_in_cds':stop,'ambiguous_dna':not set(dna)<=set('ACGT') if dna is not None else '',
                'cds_sequence_sha256':hashlib.sha256(dna.encode()).hexdigest() if dna is not None else '',
                'protein_sequence_sha256':hashlib.sha256(protein.encode()).hexdigest()}
            rows.append(row)
            if status=='exact_translation':verified[pid]=dna
        target=args.output/(name+'.verified_cds.fna')
        with target.open('w') as handle:
            for pid,dna in sorted(verified.items()):handle.write('>'+pid+'\n'+dna+'\n')
        extra=args.output/(name+'.cds_without_protein.json')
        extra.write_text(json.dumps(sorted(set(cds)-set(proteins)),indent=2)+'\n')
        local_markers=[r for r in markers if r['taxon_id']==name]
        summaries.append({'taxon_id':name,'species_name':taxon['species_name'],'proteins_screened':len(proteins),
            'source_cds_records':len(cds),'verified_cds_records':len(verified),'status_counts':dict(Counter(r['status'] for r in rows)),
            'cds_without_protein':len(set(cds)-set(proteins)),'recovered_single_copy_markers':len(local_markers),
            'markers_with_verified_cds':sum(r['protein_id'] in verified for r in local_markers),
            'source_cds_sha256':sha(source),'normalized_protein_sha256':sha(protein_path),
            'verified_cds_path':str(target),'verified_cds_sha256':sha(target)})
        all_rows.extend(rows)
    table=args.output/'protein_cds_audit.tsv'
    with table.open('w') as handle:
        writer=csv.DictWriter(handle,list(all_rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(all_rows)
    result={'status':'complete_published_outgroup_cds_translation_audit','taxa':summaries,
        'translation_table_tested':1,'proteins_screened':len(all_rows),'status_counts':dict(Counter(r['status'] for r in all_rows)),
        'external_inventory_sha256':sha(inventory_path),'qc_input_receipts_sha256':sha(input_path),
        'external_source_receipts_sha256':sha(source_path),'marker_mapping_sha256':sha(marker_path),'script_sha256':sha(Path(__file__)),
        'interpretation':'Exact CDS/protein identifiers; complete table-1 translation compared to normalized proteins, removing at most one terminal translated stop. DNA stop codons retained and flagged. No identifier heuristics, frame shifts, sequence trimming or alternative-code search. Non-triplet/mismatching entries remain explicit and excluded from verified CDS. Published CDS/genome coordinate agreement, codon alignment and selection suitability are separate checks.',
        'artifacts':{p.name:sha(p) for p in args.output.iterdir()}}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
