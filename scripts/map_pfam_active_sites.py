#!/usr/bin/env python3
"""Project Pfam active-site HMM positions to marker proteins, retaining ambiguity."""
import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from Bio import AlignIO,SeqIO
from assess_pae_sensitivity import checked_receipt
from prepare_pfam import ROOT,digest
from prepare_paired_phylogenetic_inputs import write_table


def parse_sites(path):
    result={};family=None
    for line in path.read_text().splitlines():
        if line.startswith('ID '):
            if family is not None:raise ValueError('Unterminated active-site family')
            family=line.split()[1]
            if family in result:raise ValueError('Repeated active-site family')
            result[family]=[]
        elif line=='//':family=None
        elif line.strip():
            fields=line.split()
            if family is None or len(fields)<2 or any(not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYU][1-9][0-9]*',x) for x in fields[1:]):raise ValueError('Invalid active-site pattern')
            sites=[(int(x[1:]),x[0]) for x in fields[1:]]
            if len({x[0] for x in sites})!=len(sites):raise ValueError('Repeated profile position within pattern')
            result[family].append((fields[0],sites))
    if family is not None:raise ValueError('Truncated active-site file')
    return result


def project(aligned,reference,start,original):
    sequence=aligned.upper().replace('.','-')
    if sequence.replace('-','')!=original:raise ValueError('Profile alignment altered protein residues')
    if len(sequence)!=len(reference):raise ValueError('Profile/reference lengths differ')
    match=0;residue=start-1;out={}
    for aa,rf in zip(sequence,reference):
        if aa!='-':residue+=1
        if rf not in '.- ':
            match+=1;out[match]=(residue if aa!='-' else None,aa)
    return out


def pattern_status(mapping,sites):
    if any(pos not in mapping for pos,aa in sites):raise ValueError('Active-site position outside profile')
    if any(mapping[pos][0] is None for pos,aa in sites):return 'pattern_has_unaligned_residue'
    if any(mapping[pos][1]!=aa for pos,aa in sites):return 'pattern_residues_not_all_conserved'
    return 'complete_pattern_conserved'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output=a.output.resolve()
    if a.output.exists():raise FileExistsError('Use new immutable active-site projection')
    annotations=ROOT/'results/domains/marker-annotations-v1';inputs=ROOT/'data/domains/marker-inputs-v1'
    checked_receipt(annotations);checked_receipt(inputs)
    release=ROOT/'metadata/pfam_release_receipt.json';pr=json.loads(release.read_text());base=ROOT/'data/pfam/38.2'
    for file in pr['files']:
        if file['uncompressed_file'] in ['active_site.dat','Pfam-A.hmm'] and digest(base/file['uncompressed_file'])!=file['uncompressed_sha256']:raise ValueError('Changed Pfam source')
    patterns=parse_sites(base/'active_site.dat')
    with (annotations/'raw_annotated_hits.tsv').open() as f:hits=[r for r in csv.DictReader(f,delimiter='\t') if r['pfam_name'] in patterns]
    with (annotations/'overlapping_hits.tsv').open() as f:overlaps={r[k] for r in csv.DictReader(f,delimiter='\t') for k in ['hit_a','hit_b']}
    sequences={r.id:str(r.seq) for r in SeqIO.parse(inputs/'sequences.faa','fasta')};groups=defaultdict(list)
    for r in hits:
        seq=sequences[r['sequence_id']]
        if 'S'+hashlib.sha256(seq.encode()).hexdigest()!=r['sequence_id'] or len(seq)!=int(r['protein_length']):raise ValueError('Changed protein identity')
        groups[r['pfam_accession']].append(r)
    a.output.mkdir(parents=True);profile_dir=a.output/'profiles';profile_dir.mkdir();block=[];found=set()
    with (base/'Pfam-A.hmm').open() as f:
        for line in f:
            block.append(line)
            if line.strip()=='//':
                acc=next(x.split()[1] for x in block if x.startswith('ACC '))
                if acc in groups:
                    if acc in found:raise ValueError('Duplicate profile')
                    found.add(acc);(profile_dir/(acc+'.hmm')).write_text(''.join(block))
                block=[]
    if block or found!=set(groups):raise ValueError('Incomplete selected profiles')
    align=shutil.which('hmmalign');version=subprocess.run([align,'-h'],capture_output=True,text=True,check=True).stdout.splitlines()[1]
    if 'HMMER 3.4' not in version:raise ValueError('Unexpected HMMER build')
    def run(item):
        acc,rows=item;folder=a.output/acc;folder.mkdir();profile=profile_dir/(acc+'.hmm');length=int(next(x.split()[1] for x in profile.read_text().splitlines() if x.startswith('LENG ')))
        source=folder/'segments.faa';sto=folder/'alignment.sto';segments={};by_id={r['hit_id']:r for r in rows}
        if len(by_id)!=len(rows):raise ValueError('Repeated domain hit')
        with source.open('w') as f:
            for r in rows:
                segment=sequences[r['sequence_id']][int(r['alignment_start'])-1:int(r['alignment_end'])]
                segments[r['hit_id']]=segment;f.write('>'+r['hit_id']+'\n'+segment+'\n')
        cmd=[align,'--amino','-o',str(sto),str(profile),str(source)]
        with (folder/'alignment.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=f,check=True)
        alignment=AlignIO.read(sto,'stockholm');ref=alignment.column_annotations['reference_annotation']
        if len(alignment)!=len(rows) or {r.id for r in alignment}!=set(by_id) or sum(x not in '.- ' for x in ref)!=length:raise ValueError('Profile alignment universe differs')
        output=[]
        for rec in alignment:
            r=by_id[rec.id];mapping=project(str(rec.seq),ref,int(r['alignment_start']),segments[rec.id])
            for number,(reference_id,sites) in enumerate(patterns[r['pfam_name']],1):
                status=pattern_status(mapping,sites);unambiguous=rec.id not in overlaps;coverage=float(r['hmm_coverage'])>=.5
                for pos,expected in sites:
                    protein_pos,observed=mapping[pos]
                    if protein_pos is not None and sequences[r['sequence_id']][protein_pos-1]!=observed:raise ValueError('Full protein residue readback differs')
                    output.append({'sequence_id':r['sequence_id'],'hit_id':rec.id,'pfam_accession':acc,'pfam_name':r['pfam_name'],'pfam_type':r['pfam_type'],'pattern_index':number,'reference_uniprot_accession':reference_id,'profile_match_position':pos,'expected_amino_acid':expected,'protein_residue_1based':protein_pos if protein_pos is not None else '', 'observed_amino_acid':observed,'pattern_status':status,'overlaps_another_ga_hit':not unambiguous,'hmm_coverage':r['hmm_coverage'],'candidate_active_site':status=='complete_pattern_conserved' and unambiguous and coverage})
        write_table(folder/'site_projection.tsv',output)
        return {'pfam_accession':acc,'hits':len(rows),'site_rows':len(output),'pattern_status_counts':dict(Counter(x['pattern_status'] for x in output)),'candidate_site_rows':sum(x['candidate_active_site'] for x in output),'command':cmd,'profile_sha256':digest(profile),'artifacts':{str(p.relative_to(a.output)):digest(p) for p in folder.iterdir() if p.is_file()}}
    with ThreadPoolExecutor(max_workers=2) as pool:runs=list(pool.map(run,sorted(groups.items())))
    combined=a.output/'site_projection.tsv'
    with combined.open('w') as target:
        for i,run in enumerate(runs):
            with (a.output/run['pfam_accession']/'site_projection.tsv').open() as source:
                header=next(source)
                if i==0:target.write(header)
                for line in source:target.write(line)
    result={'status':'complete_marker_pfam_active_site_projection','annotation_receipt_sha256':digest(annotations/'receipt.json'),'input_receipt_sha256':digest(inputs/'receipt.json'),'pfam_release_receipt_sha256':digest(release),'active_site_sha256':digest(base/'active_site.dat'),'families':len(runs),'domain_hits':len(hits),'unique_proteins':len({r['sequence_id'] for r in hits}),'site_rows':sum(r['site_rows'] for r in runs),'candidate_site_rows':sum(r['candidate_site_rows'] for r in runs),'aligner_sha256':digest(Path(align)),'aligner_version':version,'script_sha256':digest(Path(__file__)),'runs':runs,'interpretation':'Candidate functional correspondence, not experimental validation. HMM positions from Pfam active_site.dat are realigned to each significant hit segment. Candidate flag requires all residues of one reference pattern conserved, no overlap with any other retained GA hit and >=50% original profile coverage. Substitution/gap and alternative overlap rows remain explicit; nonconservation does not establish inactivity. Multiple reference patterns can map to the same site and are not independent observations. Structural confidence, gene history and site-specific evolutionary tests remain separate.','artifacts':{'site_projection.tsv':digest(combined)}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['runs','artifacts']},indent=2))


if __name__=='__main__':main()
