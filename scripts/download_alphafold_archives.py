#!/usr/bin/env python3
"""Retrieve version-pinned AFDB archives with range resume and publisher MD5 checks."""
import base64,csv,fcntl,hashlib,json,os,shutil,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import Request,urlopen
ROOT=Path(__file__).resolve().parents[1]
def hashes(path):
 md5=hashlib.md5();sha=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):md5.update(b);sha.update(b)
 return base64.b64encode(md5.digest()).decode(),sha.hexdigest()
def fetch(obj):
 folder=ROOT/'data/alphafold_archives';path=folder/Path(obj['name']).name
 partial=path.with_suffix('.partial');receipt=path.with_suffix('.receipt.json')
 expected=int(obj['size'])
 try:
  if path.exists():
   md5,sha=hashes(path)
   if path.stat().st_size!=expected or md5!=obj['md5Hash']:raise ValueError('Existing archive checksum mismatch')
  else:
   for attempt in range(4):
    offset=partial.stat().st_size if partial.exists() else 0
    if offset>expected:raise ValueError('Partial larger than expected object')
    if offset==expected:break
    try:
     if shutil.disk_usage(folder).free < expected-offset+50*1024**3:raise RuntimeError('Less than 50 GiB safety margin')
     headers={'Range':f'bytes={offset}-'} if offset else {}
     with urlopen(Request(obj['mediaLink'],headers=headers),timeout=90) as response:
      if offset:
       if response.status!=206 or not response.headers.get('Content-Range','').startswith(f'bytes {offset}-'):raise ValueError('Server did not honor resume range')
      with partial.open('ab' if offset else 'wb') as out:
       while True:
        b=response.read(4*1024*1024)
        if not b:break
        out.write(b)
     if partial.stat().st_size!=expected:raise ValueError('Incomplete transfer')
     break
    except Exception:
     if attempt==3:raise
     time.sleep(2**attempt)
   md5,sha=hashes(partial)
   if partial.stat().st_size!=expected or md5!=obj['md5Hash']:raise ValueError('Downloaded archive checksum mismatch')
   os.replace(partial,path)
  r={'name':obj['name'],'generation':obj['generation'],'url':obj['mediaLink'],'bytes':expected,'md5_base64':md5,'sha256':sha,'path':str(path.relative_to(ROOT)),'verified_at_utc':datetime.now(timezone.utc).isoformat(),'status':'verified'}
  receipt.write_text(json.dumps(r,indent=2)+'\n');return r
 except Exception as e:return {'name':obj['name'],'status':'error','error':str(e)}
def main():
 folder=ROOT/'data/alphafold_archives';folder.mkdir(parents=True,exist_ok=True)
 lock=(folder/'.download.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 selected=json.loads((ROOT/'metadata/alphafold_selected_archives.json').read_text())['selected_objects']
 taxa={r['taxon_id'] for r in csv.DictReader((ROOT/'metadata/analysis_manifest.tsv').open(),delimiter='\t')}
 ids={tid for r in csv.DictReader((ROOT/'metadata/alphafold_bulk_coverage.tsv').open(),delimiter='\t') if r['taxon_id'] in taxa for tid in r['queried_taxids'].split(';') if tid}
 selected=[o for o in selected if o['name'].split('tax_id-')[1].split('-')[0] in ids]
 # Small archives first to enable sequence reconciliation during acquisition.
 selected.sort(key=lambda o:int(o['size']))
 print('Archives:',len(selected),'bytes:',sum(int(o['size']) for o in selected),flush=True)
 with (ROOT/'data/raw/alphafold_downloads.jsonl').open('a') as log,ThreadPoolExecutor(max_workers=2) as pool:
  for n,f in enumerate(as_completed([pool.submit(fetch,o) for o in selected]),1):
   r=f.result();log.write(json.dumps(r)+'\n');log.flush();print(n,r['name'],r['status'],flush=True)
if __name__=='__main__':main()
