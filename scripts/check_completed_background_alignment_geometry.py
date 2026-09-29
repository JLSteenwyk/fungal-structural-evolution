#!/usr/bin/env python3
"""Record complete background geometry validation without clearing RMSD flags."""
import csv,json,subprocess
from collections import Counter
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha

def main():
 states={}
 for unit in ['fungal-background-alignment-geometry-20260928.service','fungal-background-geometry-readback-20260928.service']:
  state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
  if state!=dict(ActiveState='inactive',Result='success',ExecMainStatus='0'):raise ValueError('Geometry checks not complete: '+unit+' '+str(state))
  states[unit]=state
 pp=Path('metadata/background_alignment_geometry_readback_plan_20260928.json');plan=json.loads(pp.read_text())
 for p,h in plan['pins'].items():
  if sha(p)!=h:raise ValueError('Changed pin '+p)
 proof=Path(plan['output']);audit=json.loads(proof.read_text())
 assert audit['status']=='passed_full_background_geometry_readback' and audit['plan_sha256']==sha(pp)
 gp=json.loads(Path(plan['source_plan']).read_text());root=Path(gp['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
 assert audit['producer_receipt_sha256']==sha(rp) and r['plan_sha256']==sha(plan['source_plan'])
 for name,h in r['artifacts'].items():assert sha(root/name)==h
 diag=Path(gp['diagnostic']);dr=json.loads((diag/'receipt.json').read_text())
 assert sha(diag/'receipt.json')==r['diagnostic_receipt_sha256']
 assert dr['scientific_eligibility'] is False and dr['plan_sha256']==sha(gp['diagnostic_plan'])
 assert sha(diag/'numeric_readback.tsv')==dr['artifacts']['numeric_readback.tsv']
 key=lambda x:(x['pair_key'],x['mask'],int(x['order']))
 expected={}
 for row in csv.DictReader((diag/'numeric_readback.tsv').open(),delimiter='\t'):
  assert key(row) not in expected;expected[key(row)]=row
 assert len(expected)==dr['numerically_checked_alignments']
 seen=set();counts=Counter();cross=Counter();degenerate={};rmsd=Counter()
 for row in csv.DictReader((root/'alignment_geometry.tsv').open(),delimiter='\t'):
  k=key(row);assert k not in seen and k in expected;seen.add(k)
  assert row['rmsd_status']==expected[k]['rmsd_status'] and row['aligned_length']==expected[k]['aligned_length']
  counts[row['mask']+':'+row['geometry_status']]+=1;rmsd[row['rmsd_status']]+=1
  cross[row['mask']+':'+row['geometry_status']+':'+row['rmsd_status']]+=1
  if row['geometry_status']=='degenerate_at_numeric_tolerance':degenerate[k]=row
 assert seen==set(expected) and len(seen)==r['alignments']==audit['alignments_checked']
 assert dict(counts)==r['counts']==audit['counts'] and dict(rmsd)==dr['rmsd_status_counts']
 assert len(audit['degenerate_alignments'])==len(degenerate)
 assert {key(x) for x in audit['degenerate_alignments']}==set(degenerate)
 for row in audit['degenerate_alignments']:
  g=degenerate[key(row)];assert int(g['aligned_length'])==row['aligned_length'] and g['rmsd_status']==row['rmsd_status']
 result=dict(status='verified_completed_full_background_alignment_geometry',alignments=len(seen),degenerate_alignments=len(degenerate),counts=dict(counts),geometry_by_rmsd_status=dict(cross),rmsd_counts=dict(rmsd),readback_sha256=sha(proof),geometry_receipt_sha256=sha(rp),diagnostic_receipt_sha256=sha(diag/'receipt.json'),script_sha256=sha(__file__),terminal_states=states,scope='Full successful-alignment membership, geometry counts, unchanged RMSD classifications and degenerate audit records verified. Coordinate reconstruction and independent quaternion verification performed by completed source audit. No biological qualification or prediction-uncertainty claim; RMSD discrepancies remain flagged.')
 with Path('metadata/background_alignment_geometry_completed_20260928.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
