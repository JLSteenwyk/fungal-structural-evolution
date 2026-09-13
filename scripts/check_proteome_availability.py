#!/usr/bin/env python3
"""Resume modest-concurrency HEAD checks; availability is not sequence validation."""
import argparse,csv,json,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError
ROOT=Path(__file__).resolve().parents[1]
def check(row):
 url=row['provisional_proteome_url']
 result={'assembly_accession':row['assembly_accession'],'species_taxid':row['species_taxid'],'url':url,'checked_at_utc':datetime.now(timezone.utc).isoformat()}
 for attempt in range(3):
  try:
   with urlopen(Request(url,method='HEAD'),timeout=20) as response:
    result.update(http_status=response.status,content_length=response.headers.get('Content-Length',''),content_type=response.headers.get('Content-Type',''),status='head_available' if response.status==200 else 'unexpected_status')
   break
  except HTTPError as e:
   result.update(http_status=e.code,status='not_found' if e.code==404 else 'http_error')
   if e.code not in (429,500,502,503,504):break
  except Exception as e:result.update(status='request_error',error=str(e))
  time.sleep(2**attempt)
 return result

def main():
 p=argparse.ArgumentParser();p.add_argument('--retry-errors',action='store_true');a=p.parse_args()
 rows=list(csv.DictReader((ROOT/'metadata/species_assembly_candidates.tsv').open(),delimiter='\t'))
 cache=ROOT/'data/raw/proteome_availability.jsonl'
 existing={}
 if cache.exists():
  for line in cache.read_text().splitlines():
   r=json.loads(line);existing[r['url']]=r
 todo=[r for r in rows if r['provisional_proteome_url'] not in existing or (a.retry_errors and existing[r['provisional_proteome_url']]['status'] in ('http_error','request_error'))]
 print('Pending:',len(todo),flush=True)
 with cache.open('a') as out,ThreadPoolExecutor(max_workers=2) as pool:
  for n,f in enumerate(as_completed([pool.submit(check,r) for r in todo]),1):
   r=f.result();out.write(json.dumps(r)+'\n');out.flush();existing[r['url']]=r
   if n%100==0:print('Checked',n,'of',len(todo),flush=True)
 fields=['assembly_accession','species_taxid','url','checked_at_utc','http_status','content_length','content_type','status','error']
 with (ROOT/'metadata/proteome_availability.tsv').open('w') as out:
  w=csv.DictWriter(out,fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for row in rows:w.writerow(existing[row['provisional_proteome_url']])
 from collections import Counter
 print(dict(Counter(existing[r['provisional_proteome_url']]['status'] for r in rows)),flush=True)
if __name__=='__main__':main()
