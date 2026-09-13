#!/usr/bin/env python3
"""Find exact sequence matches to UniProt proteins carrying AlphaFoldDB crossrefs."""
import csv,gzip,hashlib,json,io,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
ROOT=Path(__file__).resolve().parents[1]
def apply_taxonomy_overrides(rows):
 path=ROOT/'config/taxonomy_overrides.json'
 if not path.exists():return rows
 overrides=json.loads(path.read_text());by_id={r['taxon_id']:r for r in rows}
 for taxon,override in overrides.items():
  if taxon not in by_id:raise ValueError('Unknown override taxon')
  row=by_id[taxon]
  if row['species_name']!=override['species_name'] or row['species_taxid'] not in ['',override['species_taxid']]:raise ValueError('Conflicting taxonomy override')
  evidence_path=ROOT/override['evidence_path']
  if hashlib.sha256(evidence_path.read_bytes()).hexdigest()!=override['evidence_sha256']:raise ValueError('Changed taxonomy evidence')
  evidence=json.loads(evidence_path.read_text())
  if evidence['status']!='verified_deposition_taxonomy_link' or evidence['taxon_id']!=taxon or evidence['species_taxid']!=override['species_taxid']:raise ValueError('Taxonomy evidence identity mismatch')
  for source in evidence['sources'].values():
   if hashlib.sha256((ROOT/source['path']).read_bytes()).hexdigest()!=source['sha256']:raise ValueError('Changed taxonomy source')
  row['species_taxid']=override['species_taxid'];row['taxonomy_evidence_path']=override['evidence_path'];row['taxonomy_evidence_sha256']=override['evidence_sha256']
 return rows

def match(row,inp):
 folder=ROOT/'data/uniprot';folder.mkdir(parents=True,exist_ok=True)
 path=folder/(row['taxon_id']+'.tsv.gz');receipt=folder/(row['taxon_id']+'.receipt.json')
 result={'taxon_id':row['taxon_id'],'species_taxid':row['species_taxid'],'status':'pending'}
 if row.get('taxonomy_evidence_path'):result.update(taxonomy_evidence_path=row['taxonomy_evidence_path'],taxonomy_evidence_sha256=row['taxonomy_evidence_sha256'])
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
  if source['url']!=url:raise ValueError('Cached UniProt query differs from current taxonomy; preserve cache and use a reviewed new snapshot')
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
 rows=apply_taxonomy_overrides(list(csv.DictReader((ROOT/'metadata/analysis_manifest.tsv').open(),delimiter='\t')))
 inputs={r['taxon_id']:r for r in json.loads((ROOT/'metadata/qc_input_receipts.json').read_text())}
 cache=ROOT/'data/raw/uniprot_matches.jsonl';done={}
 if cache.exists():
  for line in cache.read_text().splitlines():
   r=json.loads(line)
   if r['status']=='matched':done[r['taxon_id']]=r
 with cache.open('a') as out,ThreadPoolExecutor(max_workers=2) as pool:
  jobs=[pool.submit(match,r,inputs[r['taxon_id']]) for r in rows if r['taxon_id'] not in done or done[r['taxon_id']]['species_taxid']!=r['species_taxid']]
  for n,f in enumerate(as_completed(jobs),1):
   r=f.result();out.write(json.dumps(r)+'\n');out.flush();print(n,r['taxon_id'],r['status'],r.get('exact_matched_proteins',r.get('error','')),flush=True)
if __name__=='__main__':main()
