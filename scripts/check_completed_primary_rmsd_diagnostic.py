import csv,json,hashlib,math
import argparse,subprocess,time
import psutil
from pathlib import Path
from collections import Counter

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
parser=argparse.ArgumentParser()
parser.add_argument('--wait-for',type=Path,required=True)
args=parser.parse_args()
launch_hash=sha(args.wait_for)
launch=json.loads(args.wait_for.read_text())
while True:
 try:
  process=psutil.Process(launch['pid'])
  if process.create_time()!=launch['created'] or process.status()==psutil.STATUS_ZOMBIE:break
  assert process.cmdline()==launch['cmdline']
 except psutil.NoSuchProcess:break
 time.sleep(30)
state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),state
assert sha(args.wait_for)==launch_hash
planpath=Path('metadata/primary_alignment_rmsd_diagnostic_plan_20260927.json')
plan=json.loads(planpath.read_text());root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
assert r['status']=='complete_primary_alignment_rmsd_diagnostic_not_scientific_acceptance' and r['scientific_eligibility'] is False
assert r['plan_sha256']==sha(planpath)==launch['plan_sha256']
for p,h in plan['pins'].items():assert sha(p)==h,p
source=json.loads(Path(plan['source_plan']).read_text());producer=Path(source['output']);pr=json.loads((producer/'receipt.json').read_text())
assert sha(producer/'receipt.json')==r['producer_receipt_sha256']
for name,h in r['artifacts'].items():assert sha(root/name)==h
assert sha(producer/'checkpoint_manifest.tsv')==pr['artifacts']['checkpoint_manifest.tsv']
expected=set();skips=0;disposition_counts=Counter()
for x in csv.DictReader((producer/'checkpoint_manifest.tsv').open(),delimiter='\t'):
 key=Path(x['path']).stem
 mask=key.rsplit('-',2)[1]
 disposition_counts[mask+':'+x['status']]+=1
 if x['status']=='aligned':
  assert key not in expected
  expected.add(key)
 else:assert x['status']=='input_unavailable';skips+=1
assert dict(disposition_counts)==pr['counts']==r['counts']
claimed={(x['pair_key'],x['mask'],x['order']):x for x in r['rmsd_discrepancies']}
assert len(claimed)==len(r['rmsd_discrepancies'])
seen=set();counts=Counter();short=Counter();bad=[];maximum=0
for row in csv.DictReader((root/'numeric_readback.tsv').open(),delimiter='\t'):
 key=f"{row['pair_key']}-{row['mask']}-{row['order']}"
 assert key not in seen;seen.add(key)
 error=float(row['rmsd_rounding_error']);assert math.isfinite(error)
 assert abs(abs(float(row['rmsd_recomputed'])-float(row['rmsd_native']))-error)<1e-12
 status='outside_printed_rounding' if error>.00501 else 'within_printed_rounding'
 assert row['rmsd_status']==status
 maximum=max(maximum,error);counts[status]+=1
 n=int(row['aligned_length']);assert n>=1
 short[str(n) if n<3 else 'at_least_3']+=1
 if status=='outside_printed_rounding':
  identity=(row['pair_key'],row['mask'],int(row['order']))
  bad.append(identity)
  assert identity in claimed and set(claimed[identity])==set(row)
  for field,value in claimed[identity].items():
   assert type(value)(row[field])==value,(identity,field)
assert seen==expected
assert len(seen)==r['numerically_checked_alignments']
assert len(seen)+skips==r['directed_dispositions']==pr['directed_dispositions']
assert len(bad)==len(r['rmsd_discrepancies'])
assert set(bad)=={(x['pair_key'],x['mask'],x['order']) for x in r['rmsd_discrepancies']}
assert maximum==r['maximum_rmsd_rounding_error']
assert dict(counts)==r['rmsd_status_counts']
result=dict(script_sha256=sha(__file__),terminal_state=state,status='verified_completed_primary_rmsd_diagnostic_not_scientific_acceptance',receipt=str(rp),receipt_sha256=sha(rp),plan_sha256=sha(planpath),numeric_rows=len(seen),explicit_unavailable=skips,rmsd_status_counts=dict(counts),aligned_length_counts=dict(short),maximum_rmsd_discrepancy=maximum,discrepancies=r['rmsd_discrepancies'],verification='All plan pins, producer binding, artifact hashes, exact successful checkpoint key set, table counts, every discrepancy field and error classification checked. No independent recomputation of all coordinates in this completion check; full diagnostic performs that reconstruction. Downstream eligibility remains false.')
out=Path('metadata/primary_alignment_rmsd_diagnostic_completed_20260927.json');out.open('x').write(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
