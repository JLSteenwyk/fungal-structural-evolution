#!/usr/bin/env python3
"""Download publisher-versioned genome bundles and verify each file checksum."""
import json,hashlib,os,gzip
from pathlib import Path
from urllib.request import urlopen
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
def hashfile(p,algorithm):
 h=hashlib.new(algorithm)
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 receipts=[]
 for spec in json.loads((ROOT/'config/external_genomes.json').read_text()):
  aid=spec['article_id'];folder=ROOT/f'data/external/{aid}';folder.mkdir(parents=True,exist_ok=True)
  meta=ROOT/f'metadata/figshare_{aid}.json'
  if not meta.exists():
   with urlopen(f'https://api.figshare.com/v2/articles/{aid}',timeout=30) as r:d=json.load(r)
   meta.write_text(json.dumps({k:d[k] for k in ('id','title','version','doi','files','license','description')},indent=2)+'\n')
  d=json.loads(meta.read_text())
  for entry in d['files']:
   name=entry['name']
   if Path(name).name!=name:raise ValueError('Unsafe publisher filename')
   path=folder/name
   if not path.exists():
    temp=path.with_suffix(path.suffix+'.partial')
    with urlopen(entry['download_url'],timeout=120) as response,temp.open('wb') as out:
     while True:
      b=response.read(1024*1024)
      if not b:break
      out.write(b)
    if hashfile(temp,'md5')!=entry['computed_md5']:raise ValueError('Publisher MD5 mismatch: '+name)
    os.replace(temp,path)
   if hashfile(path,'md5')!=entry['computed_md5']:raise ValueError('Cached MD5 mismatch: '+name)
   receipt=dict(spec,article_version=d['version'],name=name,url=entry['download_url'],path=str(path.relative_to(ROOT)),sha256=hashfile(path,'sha256'),publisher_md5=entry['computed_md5'],bytes=path.stat().st_size,verified_at_utc=datetime.now(timezone.utc).isoformat())
   if 'pep.fasta' in name or 'proteome.fasta' in name:
    opener=gzip.open if name.endswith('.gz') else open
    lengths=[];ids=set();current=None;terminal_dots=0;internal_dots=0
    with opener(path,'rt') as f:
     for line in f:
      line=line.strip()
      if not line:continue
      if line.startswith('>'):
       if current is not None:
        if not current:raise ValueError('Empty protein')
        lengths.append(len(current.rstrip('.')));terminal_dots+=current.endswith('.');internal_dots+='.' in current.rstrip('.')
       pid=line[1:].split()[0]
       if pid in ids:raise ValueError('Duplicate protein ID')
       ids.add(pid);current=''
      else:
       if current is None or set(line.upper())-set('ACDEFGHIKLMNPQRSTVWYBXZJUO*.'):raise ValueError('Invalid protein FASTA')
       current+=line.upper()
    if not current:raise ValueError('Missing protein sequence')
    lengths.append(len(current.rstrip('.')));terminal_dots+=current.endswith('.');internal_dots+='.' in current.rstrip('.')
    receipt.update(proteins=len(lengths),residues=sum(lengths),max_length=max(lengths),fasta_status='validated' if not (terminal_dots or internal_dots) else 'noncanonical_dot_markers',terminal_dot_records=terminal_dots,internal_dot_records=internal_dots)
   receipts.append(receipt)
   (ROOT/'metadata/external_genome_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
   print(aid,name,'verified',receipt.get('proteins',''),flush=True)
if __name__=='__main__':main()
