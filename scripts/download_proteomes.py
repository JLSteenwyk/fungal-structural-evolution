#!/usr/bin/env python3
"""Download full-scale QC proteomes; verify publisher MD5 and parse every sequence."""
import argparse,csv,gzip,hashlib,json,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import urlopen
ROOT=Path(__file__).resolve().parents[1]
def digest(path,algorithm):
 h=hashlib.new(algorithm)
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def download(row):
 result={'taxon_id':row['taxon_id'],'assembly_accession':row['assembly_accession'],'url':row['proteome_url'],'checked_at_utc':datetime.now(timezone.utc).isoformat()}
 path=ROOT/'data/proteomes'/ (row['assembly_accession']+'_protein.faa.gz');path.parent.mkdir(parents=True,exist_ok=True)
 temp=path.with_suffix('.partial')
 try:
  checksum_url=row['proteome_url'].rsplit('/',1)[0]+'/md5checksums.txt'
  with urlopen(checksum_url,timeout=30) as r:checksums=r.read().decode()
  filename=row['proteome_url'].rsplit('/',1)[-1]
  matches=[line.split()[0] for line in checksums.splitlines() if len(line.split())==2 and line.split()[1].lstrip('./')==filename]
  if len(matches)!=1:raise ValueError('Missing or ambiguous publisher checksum')
  expected=matches[0]
  if not path.exists():
   with urlopen(row['proteome_url'],timeout=60) as response,temp.open('wb') as out:
    while True:
     chunk=response.read(1024*1024)
     if not chunk:break
     out.write(chunk)
   if digest(temp,'md5')!=expected:raise ValueError('Downloaded file MD5 mismatch')
   temp.replace(path)
  if digest(path,'md5')!=expected:raise ValueError('Cached file MD5 mismatch')
  lengths=[];ids=set();length=None;unknown=0
  with gzip.open(path,'rt') as f:
   for line in f:
    line=line.strip()
    if not line:continue
    if line.startswith('>'):
     if length is not None:
      if length==0:raise ValueError('Empty sequence')
      lengths.append(length)
     name=line[1:].split()[0]
     if name in ids:raise ValueError('Duplicate protein ID')
     ids.add(name);length=0
    else:
     if length is None:raise ValueError('FASTA sequence before header')
     if set(line.upper())-set('ACDEFGHIKLMNPQRSTVWYBXZJUO*'):raise ValueError('Invalid protein alphabet')
     length+=len(line);unknown+=line.upper().count('X')
  if length is None or length==0:raise ValueError('No proteins or empty last sequence')
  lengths.append(length)
  result.update(status='validated',path=str(path.relative_to(ROOT)),publisher_md5=expected,sha256=digest(path,'sha256'),compressed_bytes=path.stat().st_size,proteins=len(lengths),residues=sum(lengths),max_length=max(lengths),unknown_residues=unknown)
 except Exception as e:result.update(status='error',error=str(e))
 finally:temp.unlink(missing_ok=True)
 return result

def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',default='metadata/fungal_sampling_draft.tsv');a=p.parse_args()
 rows=list(csv.DictReader((ROOT/a.manifest).open(),delimiter='\t'))
 cache=ROOT/'data/raw/proteome_downloads.jsonl';done={}
 if cache.exists():
  for line in cache.read_text().splitlines():
   r=json.loads(line)
   if r['status']=='validated':done[r['url']]=r
 todo=[r for r in rows if r['proteome_url'] not in done or not (ROOT/done[r['proteome_url']]['path']).exists()]
 print('Pending downloads:',len(todo),flush=True)
 with cache.open('a') as out,ThreadPoolExecutor(max_workers=2) as pool:
  for n,future in enumerate(as_completed([pool.submit(download,r) for r in todo]),1):
   r=future.result();out.write(json.dumps(r)+'\n');out.flush();done[r['url']]=r
   if n%10==0:print('Processed',n,'of',len(todo),flush=True)
 selected=[done[r['proteome_url']] for r in rows]
 (ROOT/('metadata/proteome_download_receipts.json' if a.manifest == 'metadata/fungal_sampling_draft.tsv' else 'metadata/'+Path(a.manifest).stem+'_download_receipts.json')).write_text(json.dumps(selected,indent=2)+'\n')
 print('Validated:',sum(r['status']=='validated' for r in selected),'of',len(rows),flush=True)
if __name__=='__main__':main()
