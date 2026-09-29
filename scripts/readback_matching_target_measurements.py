#!/usr/bin/env python3
"""Check every target/fit field and endpoint against original nodes and alignment queue."""
import csv,gzip,hashlib,json,subprocess,time
from collections import Counter
from pathlib import Path
import psutil
from screen_duplication_alignment_reuse import sha

def main():
 bindings={}
 def load(path):path=Path(path);bindings[str(path)]=sha(path);return json.loads(path.read_text())
 launch=load('metadata/matching_target_measurements_launch_20260928.json')
 while psutil.pid_exists(launch['pid']):
  try:
   p=psutil.Process(launch['pid'])
   if abs(p.create_time()-launch['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
   assert p.cmdline()==launch['cmdline']
  except psutil.NoSuchProcess:break
  time.sleep(30)
 state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 root=Path('results/structural_comparisons/matching-target-measurements-20260928-v1');r=load(root/'receipt.json');assert r['status']=='complete_matching_target_measurement_links_pending_independent_readback'
 for path,h in r['source_hashes'].items():assert sha(path)==h;bindings[path]=h
 for name,h in r['artifacts'].items():assert sha(root/name)==h;bindings[str(root/name)]=h
 native=load('results/structural_comparisons/duplication-alignments-20260926-v1/receipt.json');planpath=Path('metadata/duplication_alignment_plan_20260926.json');plan=load(planpath)
 assert native['plan_sha256']==sha(planpath)
 primary=load('results/structural_comparisons/primary-usable-orders-20260927-v1/receipt.json');assert primary['native_receipt_sha256']==sha('results/structural_comparisons/duplication-alignments-20260926-v1/receipt.json')
 queue=Path(plan['queue']);qr=load(queue/'receipt.json')
 for name in ['receipt.json','model_pairs.tsv','models.jsonl']:
  path=queue/name;assert sha(path)==plan['pins'][str(path)];bindings[str(path)]=sha(path)
 pairs={}
 for row in csv.DictReader((queue/'model_pairs.tsv').open(),delimiter='\t'):
  assert row['pair_key'] not in pairs;pairs[row['pair_key']]=[(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))]
 models={}
 for line in (queue/'models.jsonl').open():
  n=json.loads(line);key=n['model_id'],n['version'];assert key not in models;models[key]=n
 nodes={}
 for line in Path('results/orthology/background-match-graph-20260927-v1/target_nodes.jsonl').open():
  n=json.loads(line);assert n['node_id'] not in nodes;nodes[n['node_id']]=n
 fits={}
 for row in csv.DictReader(Path('results/structural_comparisons/primary-usable-orders-20260927-v1/pair_mask_order_summary.tsv').open(),delimiter='\t'):
  key=row['pair_key'],row['mask'];assert key not in fits;fits[key]=row
 fields=list(next(iter(fits.values())));seen=set();counts=Counter();orientation=Counter();used=set();sourcevalues=0
 with gzip.open(root/'target_mask_measurements.tsv.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   key=row['node_id'],row['mask'];assert key not in seen and key[1] in ['full','plddt70'];seen.add(key);n=nodes[key[0]]
   assert all(row[k]==str(v) for k,v in n.items());ends=[(n['model_id_'+x],n['version_'+x]) for x in ['a','b']]
   for side,end in zip(['a','b'],ends):
    model=models[end]
    for field in ['sha256','sequence_sha256','length','mean_ca_plddt','fraction_ca_plddt_below50']:assert n[field+'_'+side]==model[field]
   same=ends[0]==ends[1];assert bool(n['same_model'])==same
   canon=ends if same else pairs[n['pair_key']];assert canon==sorted(ends)
   assert hashlib.sha256(json.dumps(canon,separators=(',',':')).encode()).hexdigest()==n['pair_key']
   expected_order='same_as_canonical' if ends==canon else 'reversed_from_canonical';assert row['target_endpoint_order']==expected_order;orientation[expected_order]+=1
   for side,end in zip(['a','b'],canon):assert row['canonical_model_'+side]==end[0] and int(row['canonical_version_'+side])==end[1]
   assert row['measurement_disposition']==('identical_model_no_alignment' if same else 'measured_pair')
   expected={k:'' for k in fields} if same else fits[n['pair_key'],key[1]]
   assert set(row)==set(n)|{'mask','measurement_disposition','target_endpoint_order','canonical_model_a','canonical_version_a','canonical_model_b','canonical_version_b'}|{'fit_'+k for k in fields}
   for k,v in expected.items():assert row['fit_'+k]==v;sourcevalues+=1
   if not same:used.add(n['pair_key'])
   counts[n['guide']+':'+key[1]+':'+('identical_model_no_alignment' if same else expected['order_summary_status'])]+=1
   if len(seen)%100000==0:print('Verified target measurement rows',len(seen),'/436946',flush=True)
 assert seen=={(k,m) for k in nodes for m in ['full','plddt70']} and len(seen)==r['target_mask_rows']==436946 and len(nodes)==r['targets']==218473
 assert used==set(pairs) and len(used)==r['distinct_measured_pairs'] and dict(counts)==r['counts']
 for path,h in bindings.items():assert sha(path)==h,path
 result=dict(status='passed_full_matching_target_measurement_readback',source_receipt_sha256=sha(root/'receipt.json'),checker_sha256=sha(__file__),producer_terminal_state=state,targets=len(nodes),target_mask_rows=len(seen),fit_fields_checked=sourcevalues,distinct_measured_pairs=len(used),endpoint_orientation_counts=dict(orientation),queue_receipt_sha256=sha(queue/'receipt.json'),scope='All node and fit fields, full row membership, blanks, dispositions, original queue endpoint order and model length/sequence/coordinate/confidence bindings checked. Neither prediction independence, outcome calibration nor evolutionary effects established.')
 with Path('metadata/matching_target_measurements_completed_readback_20260928.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result),flush=True)
if __name__=='__main__':main()
