#!/usr/bin/env python3
"""Align ordered TFIIB domain pairs with explicit full-protein residue maps."""
import argparse,csv,json,shutil,subprocess
from collections import Counter,defaultdict
from pathlib import Path
from Bio import AlignIO,SeqIO
from audit_busco_gene_copies import ROOT,sha,read_table


def choose_pair(hits,minimum=.7):
    selected=[x for x in hits if x['pfam_accession']=='PF00382.25']
    if len(selected)!=2:return None,'requires_exactly_two_TFIIB_hits'
    selected=sorted(selected,key=lambda x:int(x['alignment_start']))
    if int(selected[0]['alignment_end'])>=int(selected[1]['alignment_start']):return None,'overlapping_TFIIB_hits'
    if any((int(x['hmm_end'])-int(x['hmm_start'])+1)/int(x['hmm_length'])<minimum for x in selected):return None,'TFIIB_profile_coverage_below_70pct'
    return selected,'eligible_ordered_pair'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable domain alignment')
    base=ROOT/'results/domains/tfiib-candidates-v1';r=json.loads((base/'receipt.json').read_text())
    if r['status']!='complete_full_tfiib_domain_candidate_collection':raise ValueError('Complete candidate universe required')
    for name,h in r['artifacts'].items():
        if sha(base/name)!=h:raise ValueError('Changed candidate source')
    focus=ROOT/'results/domains/tfiib-family-search-v1';fr=json.loads((focus/'receipt.json').read_text())
    if sha(focus/'receipt.json')!=r['focus_receipt_sha256'] or sha(focus/'selected_profiles.hmm')!=fr['artifacts']['selected_profiles.hmm']:raise ValueError('Changed source profiles')
    rows=read_table(base/'candidate_protein_mapping.tsv');hits=defaultdict(list)
    for x in read_table(base/'domain_hits.tsv'):hits[x['sequence_id']].append(x)
    seqs={x.id:str(x.seq) for x in SeqIO.parse(base/'candidate_proteins.faa','fasta')};audit=[];selected={};segments=[{},{}]
    for row in rows:
        entry=row['entry_id'];pair,status=choose_pair(hits[row['sequence_id']]);out=dict(row,pair_status=status)
        if pair:
            sequence=seqs[entry]
            parts=[sequence[int(x['alignment_start'])-1:int(x['alignment_end'])] for x in pair]
            if any(not set(s)<=set('ACDEFGHIKLMNPQRSTVWYX') for s in parts):out['pair_status']='nonstandard_domain_residues'
            else:
                selected[entry]=pair
                for i in range(2):segments[i][entry]=parts[i]
        audit.append(out)
    if not selected:raise ValueError('No eligible domain pairs')
    a.output.mkdir(parents=True);profile=a.output/'TFIIB.hmm';fetch=shutil.which('hmmfetch');align=shutil.which('hmmalign')
    with profile.open('w') as f:subprocess.run([fetch,str(focus/'selected_profiles.hmm'),'PF00382.25'],stdout=f,check=True)
    length=int(next(x.split()[1] for x in profile.read_text().splitlines() if x.startswith('LENG ')))
    version=subprocess.run([align,'-h'],capture_output=True,text=True,check=True).stdout.splitlines()[1]
    if 'HMMER 3.4' not in version:raise ValueError('Expected HMMER 3.4')
    projected=[{},{}];positions=[];runs=[]
    for repeat in range(2):
        source=a.output/('repeat'+str(repeat+1)+'.faa');sto=source.with_suffix('.sto')
        with source.open('w') as f:
            for entry,seq in sorted(segments[repeat].items()):f.write('>'+entry+'\n'+seq+'\n')
        command=[align,'--amino','-o',str(sto),str(profile),str(source)]
        with (a.output/('repeat'+str(repeat+1)+'.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=log,check=True)
        aln=AlignIO.read(sto,'stockholm');ref=aln.column_annotations['reference_annotation'];columns=[i for i,x in enumerate(ref) if x not in '.- ']
        if len(columns)!=length or len(aln)!=len(selected) or {x.id for x in aln}!=set(selected):raise ValueError('Profile column or identity mismatch')
        for rec in aln:
            entry=rec.id;seq=str(rec.seq).upper().replace('.','-')
            if seq.replace('-','')!=segments[repeat][entry]:raise ValueError('Profile alignment changed source domain residues')
            projected[repeat][entry]=''.join(seq[i] for i in columns);index=0;residue_positions={}
            for col,aa in enumerate(seq):
                if aa!='-':index+=1;residue_positions[col]=int(selected[entry][repeat]['alignment_start'])+index-1
            for match,col in enumerate(columns,1):positions.append({'entry_id':entry,'repeat_order':repeat+1,'profile_match_position':match,'concatenated_column':repeat*length+match,'full_protein_position':residue_positions.get(col,''),'amino_acid':seq[col]})
        runs.append({'repeat_order':repeat+1,'command':command,'input_sha256':sha(source),'stockholm_sha256':sha(sto),'match_columns_1based':[i+1 for i in columns]})
    target=a.output/'paired_domains.faa'
    with target.open('w') as f:
        for entry in sorted(selected):f.write('>'+entry+'\n'+projected[0][entry]+projected[1][entry]+'\n')
    for name,data in [('candidate_pair_audit.tsv',audit),('protein_residue_mapping.tsv',positions)]:
        with (a.output/name).open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    result={'status':'complete_ordered_tfiib_domain_pair_alignment','candidate_proteins':len(rows),'aligned_proteins':len(selected),'aligned_taxa':len({x['taxon_id'] for x in rows if x['entry_id'] in selected}),'profile_columns_per_repeat':length,'concatenated_columns':2*length,'pair_status_counts':dict(Counter(x['pair_status'] for x in audit)),'candidate_receipt_sha256':sha(base/'receipt.json'),'profile_sha256':sha(profile),'aligner_sha256':sha(Path(align)),'aligner_version':version,'runs':runs,'script_sha256':sha(Path(__file__)),'interpretation':'Exactly two nonoverlapping GA TFIIB hits, each covering at least 70 percent of the Pfam HMM, are aligned separately in N-to-C order and concatenated. All candidates/exclusions and protein-residue mappings remain explicit. Positional repeat homology is a working correspondence hypothesis; domain order, gene conversion and repeat-specific phylogenetic sensitivity require assessment. No orthogroup assignment, rooted duplication history, positive selection or resolved species phylogeny is implied.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['runs','artifacts']},indent=2))

if __name__=='__main__':main()
