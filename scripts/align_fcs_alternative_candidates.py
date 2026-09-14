#!/usr/bin/env python3
"""Align alternative proteins to original HMM match states and frozen phylogenetic sites."""
import argparse,csv,json,hashlib,shutil,subprocess
from collections import defaultdict
from pathlib import Path
from Bio import AlignIO,SeqIO
AA=set('ACDEFGHIKLMNPQRSTVWY')


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['candidates','matrix','profiles','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable output')
    cr=json.loads((a.candidates/'receipt.json').read_text());mr=json.loads((a.matrix/'receipt.json').read_text())
    for root,r in [(a.candidates,cr),(a.matrix,mr)]:
        for name,digest in r['artifacts'].items():
            if sha(root/name)!=digest:raise ValueError('Changed candidate/matrix artifact')
    records={x.id:str(x.seq) for x in SeqIO.parse(a.candidates/'candidate_proteins.faa','fasta')};review=table(a.candidates/'candidate_review.tsv');by_marker=defaultdict(list)
    for row in review:by_marker[row['marker']].append(row)
    sites=defaultdict(list)
    for row in table(a.matrix/'site_mapping.tsv'):sites[row['marker']].append(row)
    baseline={x.id:str(x.seq) for x in SeqIO.parse(a.matrix/'matrix.faa','fasta')}
    exe=Path(shutil.which('hmmalign'));version=subprocess.run([str(exe),'-h'],capture_output=True,text=True,check=True).stdout.splitlines()[1]
    a.output.mkdir(parents=True);summaries=[];proofs=[];residue_rows=[]
    for marker,local in sorted(by_marker.items()):
        pr=json.loads((a.profiles/(marker+'.receipt.json')).read_text());model=Path('data/busco_downloads/lineages/eukaryota_odb12.2/hmms')/(marker+'.hmm')
        if sha(model)!=pr['profile_sha256'] or version!=pr['version']:raise ValueError('Different profile model/version')
        folder=a.output/marker;folder.mkdir();source=folder/'candidate_unaligned.faa'
        with source.open('w') as f:
            for row in local:f.write('>'+row['review_sequence_id']+'\n'+records[row['review_sequence_id']]+'\n')
        sto=folder/'candidates.sto';command=[str(exe),'--amino','-o',str(sto),str(model),str(source)]
        with (folder/'alignment.log').open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
        alignment=AlignIO.read(sto,'stockholm');rf=alignment.column_annotations['reference_annotation'];match=[i for i,c in enumerate(rf) if c not in '.- ']
        length=int(next(line.split()[1] for line in model.read_text().splitlines() if line.startswith('LENG ')))
        if len(match)!=length or {x.id for x in alignment}!={x['review_sequence_id'] for x in local}:raise ValueError('Profile coordinate/identity grid differs')
        model_to_sto={i+1:j for i,j in enumerate(match)};columns=sites[marker];selected=[model_to_sto[int(x['alignment_column_1based'])] for x in columns];projected={};original_positions={}
        for rec in alignment:
            sequence=str(rec.seq).upper().replace('.','-')
            if sequence.replace('-','')!=records[rec.id]:raise ValueError('Stockholm changed original residues')
            pos=[];count=0
            for char in sequence:
                if char!='-':count+=1;pos.append(count)
                else:pos.append(None)
            projected[rec.id]=''.join(sequence[i] if sequence[i] in AA else '?' for i in selected)
            original_positions[rec.id]=pos
            for j,(slot,col) in enumerate(zip(selected,columns),1):
                aa=sequence[slot];position=pos[slot]
                if position is not None and records[rec.id][position-1]!=aa:raise ValueError('Residue mapping differs')
                residue_rows.append({'marker':marker,'review_sequence_id':rec.id,'marker_matrix_column_1based':j,'matrix_column_1based':col['matrix_column_1based'],'profile_match_state_1based':col['alignment_column_1based'],'protein_residue_1based':position if position is not None else '', 'amino_acid':aa,'analysis_state':aa if aa in AA else '?'})
        ids={row['review_sequence_id']:row for row in local}
        with (folder/'candidate_projected.faa').open('w') as f:
            for ident,seq in projected.items():f.write('>'+ident+'\n'+seq+'\n')
        matrix_indices=[int(x['matrix_column_1based'])-1 for x in columns];full={taxon:''.join(seq[i] for i in matrix_indices) for taxon,seq in baseline.items()}
        with (folder/'baseline_plus_candidates.faa').open('w') as f:
            for ident,seq in full.items():f.write('>'+ident+'\n'+seq+'\n')
            for ident,seq in projected.items():
                if ident in full:raise ValueError('Candidate collides with baseline taxon identifier')
                f.write('>'+ident+'\n'+seq+'\n')
        reread={x.id:str(x.seq) for x in SeqIO.parse(folder/'baseline_plus_candidates.faa','fasta')}
        if reread!={**full,**projected}:raise ValueError('Augmented alignment readback differs')
        for ident,seq in projected.items():
            row=ids[ident];summaries.append({'marker':marker,'taxon_id':row['taxon_id'],'protein_id':row['protein_id'],'review_sequence_id':ident,'matrix_marker_columns':len(seq),'observed_canonical_residues':sum(c!='?' for c in seq),'observed_fraction':sum(c!='?' for c in seq)/len(seq),'cds_audit_statuses':row['cds_audit_statuses'],'baseline_taxa_preserved':len(full),'review_status':'candidate_placement_and_copy_review_pending'})
        proofs.append({'marker':marker,'profile_model_sha256':sha(model),'baseline_profile_receipt_sha256':sha(a.profiles/(marker+'.receipt.json')),'command':command,'artifacts':{p.name:sha(p) for p in folder.iterdir()}})
    for name,rows in [('candidate_alignment_review.tsv',summaries),('candidate_residue_mapping.tsv',residue_rows)]:
        with (a.output/name).open('w',newline='') as f:
            w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    r={'status':'complete_candidate_profile_alignment_and_frozen_site_projection','candidates':len(summaries),'markers':len(by_marker),'projected_candidate_site_rows':len(residue_rows),'observed_canonical_candidate_residues':sum(x['observed_canonical_residues'] for x in summaries),'candidate_source_receipt_sha256':sha(a.candidates/'receipt.json'),'baseline_matrix_receipt_sha256':sha(a.matrix/'receipt.json'),'script_sha256':sha(Path(__file__)),'hmmalign_executable_sha256':sha(exe),'hmmalign_version':version,'marker_proofs':proofs,'artifacts':{p.name:sha(p) for p in a.output.iterdir() if p.is_file()},'interpretation':'Candidate proteins aligned against the same marker HMMs and projected to the unchanged baseline matrix columns. Full sequence and residue-coordinate readback passed; original 526 baseline entries retained in augmented review files, including fully missing per-marker rows. No phylogeny inferred, species count increased, orthology accepted, candidate replaced, CDS frame repaired or structure predicted.'}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ['marker_proofs','artifacts']},indent=2))


if __name__=='__main__':main()
