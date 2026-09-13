#!/usr/bin/env python3
"""Record protozoan catalogs and explicit candidate sources for close outgroups."""
import csv,hashlib,json
from pathlib import Path
from catalog_inventory import parse_catalog
ROOT=Path(__file__).resolve().parents[1]
receipts_path=ROOT/'metadata/source_receipts.json'
receipts=json.loads(receipts_path.read_text())
records=[]
terms=['nuclear','fonticula','parvular','capsaspora','chromosphaera','monosiga','salpingoeca','sphaeroforma','creolimax','thecamonas','acanthamoeba','dictyostelium','protostelium','polysphondylium','ichthyophonus','abeoforma','pirum','corallochytrium','ministeria']
for source in ['genbank','refseq']:
 path=ROOT/f'data/raw/{source}_protozoa_assembly_summary.txt'
 data=path.read_bytes()
 receipts[f'{source}_protozoa']={'url':f'https://ftp.ncbi.nlm.nih.gov/genomes/{source}/protozoa/assembly_summary.txt','path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)}
 for row in parse_catalog(data.decode()):
  if any(t in row['organism_name'].lower() for t in terms):
   row['catalog_source']=source
   row['status']='candidate_taxonomy_and_annotation_review_required'
   records.append(row)
with (ROOT/'metadata/outgroup_assembly_candidates.tsv').open('w') as f:
 w=csv.DictWriter(f,list(records[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(records)
receipts_path.write_text(json.dumps(receipts,indent=2)+'\n')
print(len(records),'assembly records from targeted close-outgroup name search; not 25 selected taxa')
