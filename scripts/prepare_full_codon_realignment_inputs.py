#!/usr/bin/env python3
"""Recover exact full source sequences for every prepared codon diagnostic case."""
import argparse,csv,json,hashlib
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fasta(p):return {r.id:str(r.seq).upper() for r in SeqIO.parse(p,'fasta')}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text());ph=sha(a.plan)
 for path,h in p['pins'].items():assert sha(path)==h,path
 root=Path(p['inputs']);receipt=json.loads((root/'receipt.json').read_text());proof=json.loads(Path(p['readback']).read_text())
 assert proof['status']=='passed_full_genus_codon_diagnostic_readback' and proof['source_receipt_sha256']==sha(root/'receipt.json')
 for name,h in receipt['artifacts'].items():assert sha(root/name)==h
 index=Path(p['source_index']);ir=json.loads((index/'receipt.json').read_text());dna=fasta(index/'marker_source_cds.fna');assert sha(index/'marker_source_cds.fna')==ir['artifacts']['marker_source_cds.fna']
 ar=json.loads(Path(p['alignment_audit']).read_text());ah={r['marker']:r['alignment_sha256'] for r in ar['alignment_receipts']}
 with open(p['case_summary']) as f:cases=[r for r in csv.DictReader(f,delimiter='\t') if r['status']=='ready_for_tree_and_divergence_diagnostics']
 assert len(cases)==receipt['ready_cases']==1712
 out=Path(p['output']);out.mkdir(parents=True,exist_ok=False);cache={};summary=[];artifacts={};ntotal=0
 for case in cases:
  cid=case['case_id'];marker=case['marker'];code=int(case['translation_table']);taxa=sorted(fasta(root/cid/'codons.fna'))
  if marker not in cache:
   path=Path(p['global_alignments'])/(marker+'.faa');assert sha(path)==ah[marker];cache[marker]=fasta(path)
  global_aa=cache[marker];columns=list(csv.DictReader(open(root/cid/'columns.tsv'),delimiter='\t'));old_dna=fasta(root/cid/'codons.fna');old_aa=fasta(root/cid/'amino_acids.faa');sequences={};positions={}
  for taxon in taxa:
   aligned=global_aa[taxon];protein=aligned.replace('-','');coding=dna[marker+'|'+taxon];assert len(coding)%3==0
   translated=str(Seq(coding).translate(table=code))
   if translated.endswith('*'):translated=translated[:-1];coding=coding[:-3]
   assert translated==protein and '*' not in protein,(cid,taxon,'full translation')
   assert len(coding)==3*len(protein)
   pos=[];residue=0
   for x in aligned:
    if x!='-':residue+=1
    pos.append(residue if x!='-' else 0)
   original=[]
   for j,row in enumerate(columns):
    g=int(row['source_mafft_protein_column_1based'])-1;i=pos[g];codon=coding[3*(i-1):3*i] if i else '???';valid=i and set(codon)<=set('ACGT') and aligned[g] in set('ACDEFGHIKLMNPQRSTVWY');expected=codon if valid else '???'
    assert old_dna[taxon][3*j:3*j+3]==expected,(cid,taxon,j)
    assert old_aa[taxon][j]==(aligned[g] if valid else '?')
    original.append(i if valid else 0)
   sequences[taxon]=(protein,coding);positions[taxon]=original
  folder=out/cid;folder.mkdir()
  for name,k in [('full_amino_acids.faa',0),('full_codons.fna',1)]:
   with (folder/name).open('w') as f:
    for t in taxa:f.write('>'+t+'\n'+sequences[t][k]+'\n')
  (folder/'original_residue_positions.json').write_text(json.dumps(positions,separators=(',',':'))+'\n')
  for path in folder.iterdir():artifacts[str(path.relative_to(out))]=sha(path)
  summary.append(dict(case_id=cid,marker=marker,genus_label=case['genus_label'],translation_table=code,taxa=len(taxa),original_codon_columns=len(columns),full_amino_acids=sum(len(v[0]) for v in sequences.values()),marker_copy_caveat=case['marker_copy_caveat'],status='ready_for_full_group_realignment',selection_eligibility='not_established'));ntotal+=len(taxa)
 with (out/'cases.tsv').open('w') as f:
  w=csv.DictWriter(f,list(summary[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
 artifacts['cases.tsv']=sha(out/'cases.tsv')
 for path,h in p['pins'].items():assert sha(path)==h,path
 r=dict(status='complete_full_codon_realignment_inputs_pending_readback',plan_sha256=ph,cases=len(summary),case_taxon_sequences=ntotal,artifacts=artifacts,scope='All 1712 existing prepared cases, exact retained taxon membership and translation codes. Full source CDS translates exactly to original unaligned source protein, one terminal stop removed explicitly. Every original emitted codon and AA is reconstructed through source-MAFFT columns and full-sequence residue positions. No new alignment, selection test, copy validation or saturation claim.')
 (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='artifacts'},indent=2))
if __name__=='__main__':main()
