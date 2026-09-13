#!/usr/bin/env python3
"""Import published sampling evidence without assigning project taxon roles."""
import csv,hashlib,io,json,zipfile
from pathlib import Path
import openpyxl
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'data/raw/timetree_data.zip'
z=zipfile.ZipFile(p)
w=openpyxl.load_workbook(io.BytesIO(z.read('Data/supplementaryTables.xlsx')),read_only=True,data_only=True)
records=[]
for sheet in w.sheetnames[:2]:
 rows=list(w[sheet].values)[2:]
 for row in rows:
  if not row[0]:continue
  if sheet.startswith('T1'):
   name,abbr,t1,t2,t3,source=row[:6]
   lineage='; '.join(str(x) for x in (t1,t2,t3) if x)
  else:
   name,abbr,lineage,source=row[:4]
  records.append(dict(species_name=name,abbreviation=abbr,published_lineage=lineage,source=source,sheet=sheet,citation='https://doi.org/10.1038/s41559-025-02851-z',status='published_candidate_not_selected'))
with (ROOT/'metadata/timetree_published_taxa.tsv').open('w') as f:
 writer=csv.DictWriter(f,list(records[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(records)
d=json.loads((ROOT/'metadata/source_receipts.json').read_text())
article=json.loads((ROOT/'metadata/timetree_source.json').read_text())
entry=article['files'][0]
assert hashlib.md5(p.read_bytes()).hexdigest()==entry['computed_md5'],'Figshare MD5 mismatch'
d['timetree_supplement']={'url':entry['download_url'],'doi':'10.6084/m9.figshare.28046594','figshare_version':article['version'],'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'publisher_md5_verified':True,'bytes':p.stat().st_size,'path':str(p.relative_to(ROOT))}
(ROOT/'metadata/source_receipts.json').write_text(json.dumps(d,indent=2)+'\n')
print('Imported',len(records),'published taxa; Figshare checksum verified')
for r in records:
 if not r['published_lineage'].startswith(('Dikarya','Fungi')):print(r['species_name'],r['published_lineage'])
