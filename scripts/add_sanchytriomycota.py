#!/usr/bin/env python3
"""Add published genome-derived annotations for two otherwise missing fungal taxa."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
base=list(csv.DictReader((ROOT/'metadata/fungal_sampling_draft.tsv').open(),delimiter='\t'))
receipts=json.loads((ROOT/'metadata/external_genome_receipts.json').read_text())
catalog=list(csv.DictReader((ROOT/'metadata/assembly_candidates_taxonomy.tsv').open(),delimiter='\t'))
for genus in ['Amoeboradix','Sanchytrium']:
 row=next(r for r in catalog if r['taxon_name'].startswith(genus+' '))
 pep=next(r for r in receipts if r['article_id']==16411629 and r['name'].startswith(genus+'_'))
 genome=next(r for r in receipts if r['article_id']==14816085 and r['name'].startswith(genus+'_'))
 assert pep['fasta_status']=='validated'
 m={k:'' for k in base[0]}
 m.update(taxon_id='F'+row['species_taxid'],species_name=row['taxon_name'],species_taxid=row['species_taxid'],study_role='ingroup',lineage='Sanchytriomycota',assembly_accession=row['assembly_accession'],annotation_source_version='Figshare 16411629 version '+str(pep['article_version']),proteome_url=pep['url'],genome_url=genome['url'],assembly_level=row['assembly_level'],contamination_status='not_assessed',ecology='unknown',structure_coverage='not_assessed',inclusion_reason='Fill missing Sanchytriomycota using genome-derived annotations from doi:10.1038/s41467-021-25308-w; strain/assembly sequence correspondence requires further audit',status='provisional_selected_for_QC')
 base.append(m)
assert len(base)==502 and len({r['species_taxid'] for r in base})==502
with (ROOT/'metadata/fungal_sampling_expanded_draft.tsv').open('w') as f:
 w=csv.DictWriter(f,list(base[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(base)
print('502 distinct fungal candidates; Sanchytriomycota proteomes verified, cross-source assembly correspondence pending')
