#!/usr/bin/env python3
"""Find exact sequence matches to UniProt proteins carrying AlphaFoldDB crossrefs."""
import csv,gzip,hashlib,json,io,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
ROOT=Path(__file__).resolve().parents[1]
def match(row,inp):
 folder=ROOT/'data/uniprot';folder.mkdir(parents=True,exist_ok=True)
 path=folder/(row['taxon_id']+'.tsv.gz');receipt=folder/(row['taxon_id']+'.receipt.json')
 result={'taxon_id':row['taxon_id'],'species_taxid':row['species_taxid'],'status':'pending'}
 try:
  if not row['species_taxid']:return dict(result,status='taxid_unresolved')
  params={'query':'taxonomy_id:'+row['species_taxid'],'format':'tsv','fields':'accession,organism_id,sequence,xref_alphafolddb'}
  url='https://rest.uniprot.org/uniprotkb/stream?'+urlencode(params)
  if not path.exists():
   for attempt in range(3):
    try:
     with urlopen(url,timeout=180) as response,gzip.open(path.with_suffix('.partial'),'wb') as out:
      release=response.headers.get('X-UniProt-Release','unknown')
      while True:
       b=response.read(1024*1024)
       if not b:break
       out.write(b)
     path.with_suffix('.partial').replace(path)
     receipt.write_text(json.dumps({'url':url,'uniprot_release':release,'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()},indent=2)+'\n')
     break
    except Exception:
     if attempt==2:raise
     time.sleep(2**attempt)
  source=json.loads(receipt.read_text())
  assert hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256']
  target=ROOT/inp['input_path'];assert hashlib.sha256(target.read_bytes()).hexdigest()==inp['sha256']
  lookup={};pid=None;seq=''
  def add():
   if pid is not None:lookup.setdefault(seq,[]).append(pid)
  for line in target.read_text().splitlines():
   if line.startswith('>'):add();pid=line[1:].split()[0];seq=''
   else:seq+=line.strip()
  add();matched=set();count=0;references=0
  with gzip.open(path,'rt') as ref,(folder/(row['taxon_id']+'.matches.tsv')).open('w') as out:
   reader=csv.DictReader(ref,delimiter='\t')
   if reader.fieldnames!=['Entry','Organism (ID)','Sequence','AlphaFoldDB']:raise ValueError('Unexpected UniProt schema '+str(reader.fieldnames))
   writer=csv.writer(out,delimiter='\t',lineterminator='\n');writer.writerow(['taxon_id','protein_id','uniprot_accession','uniprot_organism_id','sequence_sha256'])
   for r in reader:
    count+=1
    if not r['AlphaFoldDB']:continue
    references+=1
    sequence=r['Sequence']
    for protein in lookup.get(sequence,[]):
     matched.add(protein);writer.writerow([row['taxon_id'],protein,r['Entry'],r['Organism (ID)'],hashlib.sha256(sequence.encode()).hexdigest()])
  return dict(result,status='matched',uniprot_records=count,records_with_afdb_crossref=references,exact_matched_proteins=len(matched),query_proteins=inp['proteins'],source=source,match_path=str((folder/(row['taxon_id']+'.matches.tsv')).relative_to(ROOT)),note='Crossrefs nominate models; current AFDB sequence and coordinates still require verification')
 except Exception as e:return dict(result,status='error',error=str(e))
def main():
 rows=list(csv.DictReader((ROOT/'metadata/analysis_manifest.tsv').open(),delimiter='\t'))
 inputs={r['taxon_id']:r for r in json.loads((ROOT/'metadata/qc_input_receipts.json').read_text())}
 cache=ROOT/'data/raw/uniprot_matches.jsonl';done={}
 if cache.exists():
  for line in cache.read_text().splitlines():
   r=json.loads(line)
   if r['status']=='matched':done[r['taxon_id']]=r
 with cache.open('a') as out,ThreadPoolExecutor(max_workers=2) as pool:
  jobs=[pool.submit(match,r,inputs[r['taxon_id']]) for r in rows if r['taxon_id'] not in done]
  for n,f in enumerate(as_completed(jobs),1):
   r=f.result();out.write(json.dumps(r)+'\n');out.flush();print(n,r['taxon_id'],r['status'],r.get('exact_matched_proteins',r.get('error','')),flush=True)
if __name__=='__main__':main()
