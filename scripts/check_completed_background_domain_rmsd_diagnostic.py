import csv,json,hashlib,math
from pathlib import Path
from collections import Counter

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
planpath=Path('metadata/background_domain_rmsd_diagnostic_plan_20260927.json')
plan=json.loads(planpath.read_text());root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
assert r['status']=='complete_background_domain_rmsd_diagnostic_not_scientific_acceptance' and r['scientific_eligibility'] is False
assert r['plan_sha256']==sha(planpath)
for p,h in plan['pins'].items():assert sha(p)==h,p
source=json.loads(Path(plan['source_plan']).read_text());producer=Path(source['output']);pr=json.loads((producer/'receipt.json').read_text())
assert sha(producer/'receipt.json')==r['producer_receipt_sha256']
for name,h in r['artifacts'].items():assert sha(root/name)==h
assert sha(producer/'checkpoint_manifest.tsv')==pr['artifacts']['checkpoint_manifest.tsv']
expected=set();skips=0
for x in csv.DictReader((producer/'checkpoint_manifest.tsv').open(),delimiter='\t'):
 if x['status']=='aligned':expected.add(Path(x['path']).stem)
 else:assert x['status']=='input_unavailable';skips+=1
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
 if status=='outside_printed_rounding':bad.append((row['pair_key'],row['mask'],int(row['order'])))
assert seen==expected
assert len(seen)==r['numerically_checked_alignments']
assert len(seen)+skips==r['directed_dispositions']==pr['directed_dispositions']
assert len(bad)==len(r['rmsd_discrepancies'])
assert set(bad)=={(x['pair_key'],x['mask'],x['order']) for x in r['rmsd_discrepancies']}
assert maximum==r['maximum_rmsd_rounding_error']
assert dict(counts)==r['rmsd_status_counts']
result=dict(status='verified_completed_background_domain_rmsd_diagnostic_not_scientific_acceptance',receipt=str(rp),receipt_sha256=sha(rp),plan_sha256=sha(planpath),numeric_rows=len(seen),explicit_unavailable=skips,rmsd_status_counts=dict(counts),aligned_length_counts=dict(short),maximum_rmsd_discrepancy=maximum,discrepancies=r['rmsd_discrepancies'],verification='All plan pins, producer binding, artifact hashes, exact successful checkpoint key set, table counts and error classification checked. No independent recomputation of all coordinates in this completion check; full diagnostic performs that reconstruction. Downstream eligibility remains false.')
out=Path('metadata/background_domain_rmsd_diagnostic_completed_20260927.json');out.open('x').write(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
