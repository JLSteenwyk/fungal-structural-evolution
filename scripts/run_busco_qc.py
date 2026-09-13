#!/usr/bin/env python3
"""Run bounded full-dataset BUSCO protein QC; resumable by verified success receipts."""
import json,subprocess,hashlib,fcntl
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def run(row):
 name=row['taxon_id'];receipt=ROOT/f'results/busco/{name}.receipt.json'
 if receipt.exists():
  old=json.loads(receipt.read_text())
  if old['returncode']==0 and old['input_sha256']==row['sha256'] and list((ROOT/f'results/busco/{name}').glob('short_summary.specific.*.json')):return old
 assert hashlib.sha256((ROOT/row['input_path']).read_bytes()).hexdigest()==row['sha256']
 cmd=['conda','run','--prefix',str(ROOT/'.cache/envs/busco'),'busco','-i',str(ROOT/row['input_path']),'-m','proteins','-l',str(ROOT/'data/busco_downloads/lineages/eukaryota_odb12.2'),'--offline','-c','4','-o',name,'--out_path',str(ROOT/'results/busco')]
 # Partial output requires manual review; never force overwrite or silently reuse it.
 if (ROOT/f'results/busco/{name}').exists():return {'taxon_id':name,'returncode':None,'status':'partial_output_requires_review'}
 with (ROOT/f'logs/busco_{name}.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT)
 result={'taxon_id':name,'input_sha256':row['sha256'],'returncode':r.returncode,'command':cmd}
 receipt.write_text(json.dumps(result,indent=2)+'\n');return result

def main():
 (ROOT/'results/busco').mkdir(parents=True,exist_ok=True)
 lock=(ROOT/'results/busco/.batch.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 rows=json.loads((ROOT/'metadata/qc_input_receipts.json').read_text());print('QC taxa:',len(rows),flush=True)
 with ThreadPoolExecutor(max_workers=4) as pool:
  for n,f in enumerate(as_completed([pool.submit(run,r) for r in rows]),1):
   r=f.result();print(n,r['taxon_id'],r['returncode'],flush=True)
if __name__=='__main__':main()
