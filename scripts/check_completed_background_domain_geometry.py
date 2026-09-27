#!/usr/bin/env python3
"""Archive the completed full background geometry audit and exact short-fit set."""
import csv,json,hashlib,subprocess
from pathlib import Path
from collections import Counter
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
pp=Path('metadata/background_domain_geometry_readback_plan_20260927.json');plan=json.loads(pp.read_text())
for p,h in plan['pins'].items():assert sha(p)==h
proof=Path(plan['output']);audit=json.loads(proof.read_text());assert audit['status']=='passed_full_background_domain_geometry_readback' and audit['plan_sha256']==sha(pp)
gp=json.loads(Path(plan['source_plan']).read_text());root=Path(gp['output']);r=json.loads((root/'receipt.json').read_text());assert audit['producer_receipt_sha256']==sha(root/'receipt.json')
for name,h in r['artifacts'].items():assert sha(root/name)==h
rows=list(csv.DictReader((root/'alignment_geometry.tsv').open(),delimiter='\t'));key=lambda x:(x['pair_key'],x['mask'],int(x['order']))
assert len(rows)==len({key(x) for x in rows})==audit['alignments_checked']==r['alignments']==267246
assert dict(Counter(x['mask']+':'+x['geometry_status'] for x in rows))==audit['counts']==r['counts']
degenerate={key(x) for x in rows if x['geometry_status']=='degenerate_at_numeric_tolerance'}
short=json.loads(Path('metadata/background_domain_short_geometry_readback_20260927.json').read_text());assert degenerate=={key(x) for x in short['rows']}=={key(x) for x in audit['degenerate_alignments']}
states={}
for unit in ['fungal-background-domain-geometry-20260927.service','fungal-background-domain-geometry-readback-20260927.service','fungal-background-domain-usable-orders-20260927.service','fungal-domain-coverage-geometry-readback-20260927.service']:
 state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines());assert state=={'ActiveState':'inactive','ExecMainStatus':'0'};states[unit]=state
result=dict(status='verified_completed_full_background_domain_geometry',alignments=len(rows),degenerate_alignments=len(degenerate),counts=audit['counts'],maximum_scaled_quaternion_curvature_error=audit['maximum_scaled_quaternion_curvature_error'],readback_sha256=sha(proof),geometry_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),terminal_states=states,scope='Complete key counts, spectra-category totals, source pins and exact agreement with all independently solved short fits checked. Full coordinate reconstruction performed by completed geometry verifier. No biological qualification.')
with Path('metadata/background_domain_geometry_completed_20260927.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
