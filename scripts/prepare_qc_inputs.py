#!/usr/bin/env python3
"""Create raw-proteome QC inputs, stripping terminal markers with an audit trail."""
import csv,gzip,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def convert(source,target):
 opener=gzip.open if source.suffix=='.gz' else open
 records=[];pid=None;seq=''
 def add():
  if pid is not None:records.append((pid,seq))
 with opener(source,'rt') as f:
  for line in f:
   line=line.strip()
   if line.startswith('>'):
    add();pid=line[1:].split()[0];seq=''
   else:seq+=line.upper()
 add()
 audit={'proteins':0,'terminal_markers_removed':0,'internal_stop_records':0,'noncanonical_records':0}
 temp=target.with_suffix('.partial')
 with temp.open('w') as out:
  for pid,seq in records:
   clean=seq.rstrip('*.');audit['terminal_markers_removed']+=len(seq)-len(clean)
   if '*' in clean:audit['internal_stop_records']+=1;continue
   if not clean or set(clean)-set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):audit['noncanonical_records']+=1;continue
   out.write('>'+pid+'\n'+clean+'\n');audit['proteins']+=1
 temp.replace(target)
 audit['sha256']=hashlib.sha256(target.read_bytes()).hexdigest()
 return audit

def main():
 wanted={r['proteome_url']:r for r in csv.DictReader((ROOT/'metadata/fungal_sampling_expanded_draft.tsv').open(),delimiter='\t')}
 downloads={}
 for line in (ROOT/'data/raw/proteome_downloads.jsonl').read_text().splitlines():
  try:r=json.loads(line)
  except json.JSONDecodeError:continue
  if r['status']=='validated':downloads[r['url']]=r
 for r in json.loads((ROOT/'metadata/external_genome_receipts.json').read_text()):
  if 'proteins' in r:downloads[r['url']]=r
 output=[];folder=ROOT/'data/qc_proteomes';folder.mkdir(parents=True,exist_ok=True)
 for url,row in wanted.items():
  if url not in downloads:continue
  source=ROOT/downloads[url]['path'];target=folder/(row['taxon_id']+'.faa')
  assert hashlib.sha256(source.read_bytes()).hexdigest()==downloads[url]['sha256'],'Raw checksum changed'
  audit=convert(source,target)
  output.append(dict(taxon_id=row['taxon_id'],species_name=row['species_name'],input_path=str(target.relative_to(ROOT)),source_sha256=downloads[url]['sha256'],**audit))
 (ROOT/'metadata/qc_input_receipts.json').write_text(json.dumps(output,indent=2)+'\n')
 print('Prepared',len(output),'taxa for raw-proteome completeness; isoform selection pending')
if __name__=='__main__':main()
