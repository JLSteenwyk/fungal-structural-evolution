#!/usr/bin/env python3
"""Independently verify complete full-atlas alias identity-map coverage."""
import argparse,csv,gzip,hashlib,json
from collections import Counter
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def lookup(p):
 s=set()
 for n,l in enumerate(open(p),1):
  x=l.rstrip('\n').split('\t')
  if len(x)!=3 or not x[1] or x[1] in s:raise ValueError(f'bad lookup {n}')
  s.add(x[1])
 return s
def main():
 a=argparse.ArgumentParser();a.add_argument('--map',type=Path,required=True);a.add_argument('--producer',type=Path,required=True);a.add_argument('--lookup',type=Path,required=True);a.add_argument('--receipt',type=Path,required=True);x=a.parse_args()
 if x.receipt.exists():raise FileExistsError(x.receipt)
 p=json.load(open(x.producer));
 if p.get('status')!='completed_full_atlas_cluster_identity_map' or p.get('output_sha256')!=sha(x.map):raise ValueError('producer map receipt mismatch')
 expected=lookup(x.lookup); observed=set();sources=Counter();taxa=set();rows=0
 with gzip.open(x.map,'rt',newline='') as f:
  r=csv.DictReader(f,delimiter='\t'); required={'foldseek_alias','taxon_id','protein_id','source','model_id','model_version','sequence_sha256','length','availability','paired_model_alias'}
  if set(r.fieldnames or [])!=required:raise ValueError('identity-map header mismatch')
  for row in r:
   alias=row['foldseek_alias']
   if alias not in expected or not row['taxon_id'] or not row['protein_id'] or row['source'] not in ('AFDB','ESMFold'):raise ValueError('invalid identity relation')
   observed.add(alias);sources[row['source']]+=1;taxa.add(row['taxon_id']);rows+=1
 if observed!=expected:raise ValueError(f'incomplete identity map {len(observed)} of {len(expected)}')
 result={'status':'passed_independent_full_atlas_cluster_identity_map_readback','models':len(observed),'identity_rows':rows,'taxa':len(taxa),'relations_by_source':dict(sources),'map_sha256':sha(x.map),'producer_receipt_sha256':sha(x.producer),'lookup_sha256':sha(x.lookup),'scientific_eligibility':False,'scope':'Independent complete alias-coverage and schema replay of the identity map; no homology or evolutionary inference.'}
 x.receipt.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
