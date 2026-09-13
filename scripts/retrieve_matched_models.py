#!/usr/bin/env python3
"""Retrieve current full-length AFDB monomers with exact sequence and coordinate checks."""
import csv,fcntl,hashlib,io,json,re,time,signal
from collections import defaultdict,deque
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
from pathlib import Path
from urllib.request import urlopen
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
ROOT=Path(__file__).resolve().parents[1]
def prioritize(links,marker_keys,done):
 groups=defaultdict(set);remaining=set()
 for row in links:
  accession=row['uniprot_accession']
  if accession in done:continue
  remaining.add(accession)
  if (row['taxon_id'],row['protein_id'],row['sequence_sha256']) in marker_keys:groups[row['taxon_id']].add(accession)
 queues=[deque(sorted(groups[taxon])) for taxon in sorted(groups)];ordered=[];seen=set()
 while any(queues):
  for queue in queues:
   while queue and queue[0] in seen:queue.popleft()
   if queue:
    accession=queue.popleft();ordered.append(accession);seen.add(accession)
 priority_count=len(ordered)
 ordered.extend(sorted(remaining-seen))
 return ordered,priority_count,len(groups)
def fetch(accession,expected):
 folder=ROOT/'data/structures/afdb';meta=folder/(accession+'.api.json')
 try:
  if not meta.exists():
   with urlopen('https://alphafold.ebi.ac.uk/api/prediction/'+accession,timeout=45) as r:data=json.load(r)
   temp=meta.with_suffix('.partial');temp.write_text(json.dumps(data));temp.replace(meta)
  records=json.loads(meta.read_text());accepted=[]
  for model in records:
   if model.get('isComplex',False):continue
   sequence=model.get('sequence','');sha=hashlib.sha256(sequence.encode()).hexdigest()
   if sha not in expected or model.get('sequenceStart')!=1 or model.get('sequenceEnd')!=len(sequence):continue
   ident=model['modelEntityId'];version=model['latestVersion']
   if not re.fullmatch(r'[A-Za-z0-9_-]+',ident):raise ValueError('Invalid model ID')
   path=folder/f'{ident}-v{version}.cif';url=model['cifUrl']
   if not url.startswith('https://alphafold.ebi.ac.uk/files/'):raise ValueError('Unexpected CIF provider')
   if not path.exists():
    with urlopen(url,timeout=60) as r:content=r.read()
   else:content=path.read_bytes()
   cif=MMCIF2Dict(io.StringIO(content.decode()))
   seqs=[''.join(s.split()) for s in cif.get('_entity_poly.pdbx_seq_one_letter_code_can',[])]
   if seqs!=[sequence]:raise ValueError('CIF polymer sequence does not match API and input')
   atoms=cif['_atom_site.label_atom_id'];positions=cif['_atom_site.label_seq_id'];bfactors=cif['_atom_site.B_iso_or_equiv']
   ca={int(pos):float(b) for atom,pos,b in zip(atoms,positions,bfactors) if atom=='CA'}
   if set(ca)!=set(range(1,len(sequence)+1)):raise ValueError('Incomplete CA coverage')
   if not path.exists():
    temp=path.with_suffix('.partial');temp.write_bytes(content);temp.replace(path)
   accepted.append({'model_id':ident,'version':version,'sequence_sha256':sha,'length':len(sequence),'mean_ca_plddt':sum(ca.values())/len(ca),'fraction_ca_plddt_below50':sum(x<50 for x in ca.values())/len(ca),'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(content).hexdigest(),'url':url,'model_created_date':model.get('modelCreatedDate'),'sequence_version_date':model.get('sequenceVersionDate'),'provider':model.get('providerId'),'tool':model.get('toolUsed'),'pae_url':model.get('paeDocUrl'),'pae_downloaded':False})
  return {'uniprot_accession':accession,'status':'verified' if accepted else 'no_exact_full_length_model','models':accepted,'api_path':str(meta.relative_to(ROOT)),'api_sha256':hashlib.sha256(meta.read_bytes()).hexdigest()}
 except Exception as e:return {'uniprot_accession':accession,'status':'error','error':str(e)}
def main():
 folder=ROOT/'data/structures/afdb';folder.mkdir(parents=True,exist_ok=True)
 lock=(folder/'.retrieval.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 matches={};links=[]
 completed={}
 for line in (ROOT/'data/raw/uniprot_matches.jsonl').read_text().splitlines():
  try:r=json.loads(line)
  except json.JSONDecodeError:continue
  if r['status']=='matched':completed[r['taxon_id']]=r
 for r in completed.values():
  with (ROOT/r['match_path']).open() as handle:
   for row in csv.DictReader(handle,delimiter='\t'):
    a=row['uniprot_accession']
    if not re.fullmatch(r'[A-Za-z0-9]+',a):raise ValueError('Invalid accession')
    matches.setdefault(a,set()).add(row['sequence_sha256']);links.append(row)
 # Retain every original protein mapping even when structures are downloaded once.
 with (folder/'input_links.tsv').open('w') as out:
  if links:
   w=csv.DictWriter(out,list(links[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(links)
 cache=ROOT/'data/raw/afdb_models.jsonl';done={}
 if cache.exists():
  for line in cache.read_text().splitlines():
   r=json.loads(line)
   if r['status']=='verified' and all((ROOT/m['path']).exists() and hashlib.sha256((ROOT/m['path']).read_bytes()).hexdigest()==m['sha256'] for m in r['models']):done[r['uniprot_accession']]=r
 with (ROOT/'results/phylogeny/markers-full-v1/protein_mapping.tsv').open() as handle:
  marker_keys={(r['taxon_id'],r['protein_id'],r['sequence_sha256']) for r in csv.DictReader(handle,delimiter='\t')}
 pending,priority_count,priority_taxa=prioritize(links,marker_keys,done)
 plan={'queued_accessions':len(pending),'priority_marker_accessions':priority_count,'priority_taxa':priority_taxa,'verified_cached_accessions':len(done),'completed_matching_taxa':len(completed),'policy':'Round-robin marker accessions across taxa, then all other exact-sequence candidates; no atlas candidates discarded.'}
 (ROOT/'metadata/structure_retrieval_queue_snapshot.json').write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps(plan),flush=True)
 stopping=[False]
 def stop(signum,frame):stopping[0]=True
 signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
 with cache.open('a') as out,ThreadPoolExecutor(max_workers=2) as pool:
  iterator=iter(pending);active=set();n=0
  def submit():
   a=next(iterator,None)
   if a is not None:active.add(pool.submit(fetch,a,matches[a]))
  for _ in range(4):submit()
  while active:
   ready,_=wait(active,return_when=FIRST_COMPLETED)
   for future in ready:
    active.remove(future);r=future.result();out.write(json.dumps(r)+'\n');out.flush();n+=1
    if n%100==0 or r['status']=='error':print(n,r['uniprot_accession'],r['status'],r.get('error',''),flush=True)
    if not stopping[0]:submit()
  print('Stopped cleanly' if stopping[0] else 'Snapshot queue completed',n,flush=True)
if __name__=='__main__':main()
