#!/usr/bin/env python3
"""Choose one provisional annotated assembly per species, not the final sample."""
import csv,json
from pathlib import Path
from collections import defaultdict
ROOT=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((ROOT/'metadata/assembly_candidates_taxonomy.tsv').open(),delimiter='\t'))
level={'Complete Genome':4,'Chromosome':3,'Scaffold':2,'Contig':1}
def score(r):
 return (r['refseq_category']=='reference genome',r['catalog_source']=='refseq_fungi',level.get(r['assembly_level'],0),r['seq_rel_date'],r['assembly_accession'])
groups=defaultdict(list)
for r in rows:
 if r['annotation_name'] not in ('','na') and r['version_status']=='latest' and r['ftp_path'] not in ('','na'):
  groups[r['species_taxid']].append(r)
selected=[]
for sid,rs in sorted(groups.items(),key=lambda x:int(x[0])):
 r=max(rs,key=score).copy()
 base=r['ftp_path'].replace('ftp://','https://').rstrip('/')
 stem=base.rsplit('/',1)[-1]
 r.update(provisional_proteome_url=f'{base}/{stem}_protein.faa.gz',selection_status='candidate_only_availability_and_quality_unverified',assembly_choice_reason='Annotated latest assembly; prefer reference then RefSeq then assembly level then release date. Quality reassessment can replace this choice.')
 selected.append(r)
with (ROOT/'metadata/species_assembly_candidates.tsv').open('w') as f:
 w=csv.DictWriter(f,list(selected[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(selected)
assert len({r['species_taxid'] for r in selected})==len(selected)
print(len(selected),'unique annotated species candidates; NOT final 500')
