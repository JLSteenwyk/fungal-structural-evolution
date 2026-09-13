#!/usr/bin/env python3
"""Restore missing source snapshots only if published project checksums match."""
import hashlib,json,os,shutil
from pathlib import Path
from urllib.request import urlopen
ROOT=Path(__file__).resolve().parents[1]
for name,entry in json.loads((ROOT/'metadata/source_receipts.json').read_text()).items():
 path=ROOT/entry['path']
 if not path.exists():
  path.parent.mkdir(parents=True,exist_ok=True)
  temp=path.with_suffix(path.suffix+'.partial')
  try:
   with urlopen(entry['url'],timeout=180) as response,temp.open('wb') as out:
    shutil.copyfileobj(response,out)
   with temp.open('rb') as inp:
    h=hashlib.sha256()
    for block in iter(lambda:inp.read(1024*1024),b''):h.update(block)
   if h.hexdigest()!=entry['sha256']:
    raise RuntimeError(f'{name}: remote source changed; obtain the archived snapshot or explicitly revise provenance')
   os.replace(temp,path)
  finally:
   temp.unlink(missing_ok=True)
 with path.open('rb') as inp:
  h=hashlib.sha256()
  for block in iter(lambda:inp.read(1024*1024),b''):h.update(block)
 if h.hexdigest()!=entry['sha256']:raise RuntimeError(f'{name}: checksum mismatch')
 print(name,'verified')
