#!/usr/bin/env python3
"""Independently verify full source exports and original residue coordinate maps."""
import argparse,csv,json,hashlib
from pathlib import Path
import numpy as np
from Bio.Data import CodonTable

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fasta(p):
 records={};key=None
 for line in Path(p).read_text().splitlines():
  if line.startswith('>'):
   key=line[1:].split()[0];assert key not in records;records[key]=''
  else:records[key]+=line.strip().upper()
 return records
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
assert r['plan_sha256']==sha(a.plan)
for path,h in plan['pins'].items():assert sha(path)==h
for name,h in r['artifacts'].items():assert sha(root/name)==h
source=fasta(Path(plan['source_index'])/'marker_source_cds.fna');original=Path(plan['inputs']);alignment_cache={};count=observed=0
with (root/'cases.tsv').open() as f:cases=list(csv.DictReader(f,delimiter='\t'))
with open(plan['case_summary']) as f:expected={x['case_id']:x for x in csv.DictReader(f,delimiter='\t') if x['status']=='ready_for_tree_and_divergence_diagnostics'}
assert len(cases)==len(expected)==1712 and {x['case_id'] for x in cases}==set(expected)
for case in cases:
 cid=case['case_id'];folder=root/cid;aa=fasta(folder/'full_amino_acids.faa');dna=fasta(folder/'full_codons.fna');old_dna=fasta(original/cid/'codons.fna');old_aa=fasta(original/cid/'amino_acids.faa');mapping=json.loads((folder/'original_residue_positions.json').read_text());table=CodonTable.unambiguous_dna_by_id[int(case['translation_table'])]
 assert set(aa)==set(dna)==set(old_dna)==set(old_aa)==set(mapping) and len(aa)==int(case['taxa'])
 for field in ['marker','genus_label','translation_table','marker_copy_caveat']:assert case[field]==expected[cid][field]
 marker=case['marker']
 if marker not in alignment_cache:alignment_cache[marker]=fasta(Path(plan['global_alignments'])/(marker+'.faa'))
 with (original/cid/'columns.tsv').open() as f:columns=[int(x['source_mafft_protein_column_1based'])-1 for x in csv.DictReader(f,delimiter='\t')]
 for taxon in aa:
  raw=source[marker+'|'+taxon];trimmed=raw[:-3] if raw[-3:] in table.stop_codons else raw;assert dna[taxon]==trimmed
  global_aa=alignment_cache[marker][taxon];assert aa[taxon]==global_aa.replace('-','') and len(dna[taxon])==3*len(aa[taxon])
  for i,amino in enumerate(aa[taxon]):
   codon=dna[taxon][3*i:3*i+3]
   if set(codon)<=set('ACGT'):assert table.forward_table[codon]==amino
  positions=np.cumsum(np.array(list(global_aa))!='-')
  derived=[]
  for j,col in enumerate(columns):
   pos=int(positions[col]);codon=dna[taxon][3*(pos-1):3*pos] if global_aa[col]!='-' else '???';called=global_aa[col] in 'ACDEFGHIKLMNPQRSTVWY' and set(codon)<=set('ACGT')
   derived.append(pos if called else 0);assert old_dna[taxon][3*j:3*j+3]==(codon if called else '???');assert old_aa[taxon][j]==(global_aa[col] if called else '?');observed+=called
  assert derived==mapping[taxon];count+=1
 assert int(case['original_codon_columns'])==len(columns) and int(case['full_amino_acids'])==sum(map(len,aa.values()))
assert count==r['case_taxon_sequences']
for name,h in r['artifacts'].items():assert sha(root/name)==h
result=dict(status='passed_full_codon_realignment_source_and_position_readback',producer_receipt_sha256=sha(rp),cases=len(cases),case_taxon_sequences=count,observed_original_codons=int(observed),script_sha256=sha(__file__),scope='Independent FASTA parser and cumulative-index reconstruction; exact full source DNA/AA exports, canonical codon-table translations, memberships, original codons and all residue coordinate maps checked. Ambiguous full-source codons retain exact source identity; their translation identity inherits producer and prior source audit. No alignment correctness or selection eligibility claim.')
a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
