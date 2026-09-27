#!/usr/bin/env python3
"""Wait for the pinned producer, then independently verify all bootstrap edges."""
import argparse,csv,json,time,fcntl
from pathlib import Path
from collections import defaultdict,Counter
from concurrent.futures import ProcessPoolExecutor
import psutil
from ecology_bootstrap_flow_readback import sha,audit_tree

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text());ph=sha(a.plan)
 def verify():
  assert sha(a.plan)==ph
  for path,h in p['pins'].items():assert sha(path)==h,path
 verify();producer=p['producer']
 while True:
  try:
   proc=psutil.Process(producer['pid'])
   if abs(proc.create_time()-producer['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
   assert proc.cmdline()==producer['cmdline']
  except psutil.NoSuchProcess:break
  time.sleep(30)
 verify();sourceplan=json.loads(Path(p['source_plan']).read_text());source=Path(sourceplan['output']);rp=source/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
 assert r['status']=='complete_ecology_bootstrap_edge_mapping_pending_readback' and r['plan_sha256']==sha(p['source_plan']) and r['trees']==2000
 for name,h in r['artifacts'].items():assert sha(source/name)==h
 with open(sourceplan['evidence']) as f:ev=list(csv.DictReader(f,delimiter='\t'))
 states={x['taxon_id']:1 if x['state']=='ectomycorrhizal' else 0 for x in ev if x['state'] in ['ectomycorrhizal','saprotrophic','asymbiotic']}
 with open(sourceplan['manifest']) as f:taxa=[x['taxon_id'] for x in csv.DictReader(f,delimiter='\t')]
 jobs=[]
 for spec in sourceplan['trees']:
  lines=[x for x in Path(spec['bootstrap']).read_text().splitlines() if x.strip()];assert len(lines)==1000
  for i,line in enumerate(lines,1):
   name=f"{spec['label']}-{i:04d}.json.gz";jobs.append((spec['label'],i,line,states,taxa,str(source),r['shards'][name]))
 assert len(jobs)==len(r['shards'])==2000
 out=Path(p['output']);out.mkdir(parents=True,exist_ok=True);lock=(out/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 # Full reruns intentionally recompute network flows; incomplete checkpoints are not proofs.
 aggregate=defaultdict(Counter);sides={};summaries=[];costs=edges=0;proofs={}
 with ProcessPoolExecutor(max_workers=p['workers']) as pool:
  for i,proof in enumerate(pool.map(audit_tree,jobs,chunksize=1),1):
   label=proof['tree_source'];summaries.extend(proof['summaries']);costs+=proof['constrained_costs'];edges+=len(proof['contributions'])
   for coding,key,status,side in proof['contributions']:
    aggregate[label,coding,key]['present']+=1;aggregate[label,coding,key][status]+=1
    if key in sides:assert sides[key]==side
    else:sides[key]=side
   checkpoint=out/f"{label}-{proof['tree_index']:04d}.json";compact={k:v for k,v in proof.items() if k!='contributions'};compact['plan_sha256']=ph;checkpoint.write_text(json.dumps(compact,indent=2)+'\n');proofs[checkpoint.name]=sha(checkpoint)
   if i%10==0:
    state=dict(checked_trees=i,total_trees=2000,checked_costs=costs);(out/'state.json').write_text(json.dumps(state)+'\n');print(json.dumps(state),flush=True)
 with (source/'tree_summaries.tsv').open() as f:observed=list(csv.DictReader(f,delimiter='\t'))
 assert len(observed)==len(summaries)==r['coding_tree_combinations']
 for actual,expected in zip(observed,summaries):assert actual=={k:str(v) for k,v in expected.items()}
 with (source/'split_frequencies.tsv').open() as f:observed=list(csv.DictReader(f,delimiter='\t'))
 assert len(observed)==len(aggregate)==r['split_summary_rows'];seen=set()
 for row in observed:
  key=row['tree_source'],row['coding'],row['split_id'];assert key not in seen;seen.add(key);c=aggregate[key]
  assert row['canonical_side']==sides[key[2]] and int(row['ensemble_trees'])==1000
  for k in ['present','required_change','optional_change','no_change_in_any_optimum']:assert int(row[k])==c[k]
 assert edges==r['edge_rows'] and costs==4*edges
 verify();assert sha(rp)==rh
 for name,h in r['artifacts'].items():assert sha(source/name)==h
 result=dict(status='passed_full_bootstrap_ecology_edge_network_flow_readback',plan_sha256=ph,producer_receipt_sha256=rh,trees=len(jobs),edge_rows=edges,constrained_costs=costs,split_summary_rows=len(aggregate),proofs=proofs,scope='Every endpoint cost on every bootstrap edge reconstructed by independent network maximum flow. Every split identity, per-tree summary, split-presence denominator and required/optional/no-change count verified. Conditional topology diagnostic only, not biological origins or ecological effects.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
