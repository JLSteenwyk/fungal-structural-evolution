#!/usr/bin/env python3
"""Inventory public AFDB bulk objects by species and selected assembly taxids."""
import csv,json,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
ROOT=Path(__file__).resolve().parents[1]
BUCKET='public-datasets-deepmind-alphafold-v4'
def query(taxid):
 result={'taxid':taxid,'bucket':BUCKET,'checked_at_utc':datetime.now(timezone.utc).isoformat(),'objects':[]}
 try:
  token=None
  while True:
   params={'prefix':f'proteomes/proteome-tax_id-{taxid}-','maxResults':1000}
   if token:params['pageToken']=token
   with urlopen('https://storage.googleapis.com/storage/v1/b/'+BUCKET+'/o?'+urlencode(params),timeout=30) as response:d=json.load(response)
   result['objects'].extend({k:o.get(k,'') for k in ('name','size','generation','md5Hash','crc32c','updated','mediaLink')} for o in d.get('items',[]))
   token=d.get('nextPageToken')
   if not token:break
  result['status']='listed'
 except Exception as e:result.update(status='error',error=str(e))
 return result

def main():
 taxa=list(csv.DictReader((ROOT/'metadata/sampling_manifest.tsv').open(),delimiter='\t'))
 assembly_taxids={}
 for filename in ('assembly_candidates.tsv','outgroup_assembly_candidates.tsv'):
  for r in csv.DictReader((ROOT/'metadata'/filename).open(),delimiter='\t'):assembly_taxids[r['assembly_accession']]=r['taxid']
 # Include selected animal assembly taxids from their complete source catalogs.
 from catalog_inventory import parse_catalog
 for name in ('refseq_invertebrate','genbank_invertebrate'):
  for r in parse_catalog((ROOT/f'data/raw/{name}_assembly_summary.txt').read_text()):assembly_taxids[r['assembly_accession']]=r['taxid']
 requested={r['taxon_id']:sorted({x for x in (r['species_taxid'],assembly_taxids.get(r['assembly_accession'],'')) if x}) for r in taxa}
 ids=sorted({x for values in requested.values() for x in values})
 cache=ROOT/'data/raw/alphafold_bulk_inventory.jsonl';done={}
 if cache.exists():
  for line in cache.read_text().splitlines():
   r=json.loads(line)
   if r['status']=='listed':done[r['taxid']]=r
 pending=[t for t in ids if t not in done];print('Taxonomy IDs to query:',len(pending),flush=True)
 with cache.open('a') as out,ThreadPoolExecutor(max_workers=2) as pool:
  for n,f in enumerate(as_completed([pool.submit(query,t) for t in pending]),1):
   r=f.result();done[r['taxid']]=r;out.write(json.dumps(r)+'\n');out.flush()
   if n%100==0:print('Listed',n,flush=True)
 summary=[];objects={}
 for taxon in taxa:
  hits={};errors=[]
  for tid in requested[taxon['taxon_id']]:
   r=done[tid]
   if r['status']=='error':errors.append(tid)
   for obj in r['objects']:hits[obj['name']]=obj;objects[obj['name']]=obj
  summary.append({'taxon_id':taxon['taxon_id'],'species_name':taxon['species_name'],'queried_taxids':';'.join(requested[taxon['taxon_id']]),'object_count':len(hits),'archive_bytes':sum(int(x['size']) for x in hits.values()),'query_errors':';'.join(errors),'status':'candidate_archives_found' if hits else ('query_error' if errors else 'no_archive_at_queried_ids' if requested[taxon['taxon_id']] else 'taxid_unresolved')})
 with (ROOT/'metadata/alphafold_bulk_coverage.tsv').open('w') as out:
  w=csv.DictWriter(out,list(summary[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
 (ROOT/'metadata/alphafold_bulk_objects.json').write_text(json.dumps({'bucket':BUCKET,'objects':list(objects.values()),'note':'Archive file suffix versions may differ from bucket name. Sequence identity to selected annotations remains unverified. Absence at queried IDs is not absence from current AFDB.'},indent=2)+'\n')
 print('Taxa with archives:',sum(r['object_count']>0 for r in summary),'unique archive bytes:',sum(int(x['size']) for x in objects.values()),flush=True)
if __name__=='__main__':main()
