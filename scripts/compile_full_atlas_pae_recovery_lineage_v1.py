#!/usr/bin/env python3
"""Compile completed PAE recovery/retry receipts into a clustering-gate lineage receipt."""
import argparse, hashlib, json, tempfile
from datetime import datetime, timezone
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('preparation','recovery','base-readback','retry','retry-readback','receipt'):p.add_argument('--'+n,type=Path)
 p.add_argument('--self-test',action='store_true');a=p.parse_args()
 if a.self_test:
  with tempfile.TemporaryDirectory() as d:
   root=Path(d); values={'preparation':dict(status='prepared_interrupted_full_atlas_pae_recovery_queue'),'recovery':dict(status='completed_full_atlas_missing_pae_recovery_attempts_pending_independent_readback'),'base-readback':dict(status='passed_order_independent_full_atlas_pae_recovery_readback',dispositions=dict(preserved_retrieval_verified=4,preserved_retrieval_failed=1,recovery_retrieval_verified=4,recovery_retrieval_failed=1)),'retry':dict(status='completed_full_atlas_pae_explicit_failure_retry_pending_independent_readback',expected_models=2,counts=dict(verified=2,failed=0)),'retry-readback':dict(status='passed_combined_full_atlas_pae_retry_closure',expected_models=10,retained_base_failures=2,retry_dispositions=dict(verified=2,failed=0),final_explicit_failures=0)}
   paths={}
   for n,v in values.items(): paths[n]=root/(n+'.json');paths[n].write_text(json.dumps(v))
   out=root/'out.json';main_args=['--preparation',str(paths['preparation']),'--recovery',str(paths['recovery']),'--base-readback',str(paths['base-readback']),'--retry',str(paths['retry']),'--retry-readback',str(paths['retry-readback']),'--receipt',str(out)]
   import sys; old=sys.argv;sys.argv=[old[0],*main_args];main();sys.argv=old
   assert json.loads(out.read_text())['final_dispositions']['retrieval_verified']==10
  print(json.dumps(dict(status='passed_pae_lineage_composition_literal_control'),indent=2));return
 if any(getattr(a,n.replace('-','_')) is None for n in ('preparation','recovery','base-readback','retry','retry-readback','receipt')):p.error('all receipt paths required')
 if a.receipt.exists():raise FileExistsError('Fresh lineage receipt required')
 prep,recovery,base,retry,closure=(json.loads(getattr(a,n.replace('-','_')).read_text()) for n in ('preparation','recovery','base-readback','retry','retry-readback'))
 if prep.get('status')!='prepared_interrupted_full_atlas_pae_recovery_queue':raise ValueError('Recovery preparation receipt required')
 if recovery.get('status')!='completed_full_atlas_missing_pae_recovery_attempts_pending_independent_readback':raise ValueError('Completed recovery receipt required')
 if base.get('status')!='passed_order_independent_full_atlas_pae_recovery_readback':raise ValueError('Passed base readback required')
 if retry.get('status')!='completed_full_atlas_pae_explicit_failure_retry_pending_independent_readback':raise ValueError('Completed retry receipt required')
 if closure.get('status')!='passed_combined_full_atlas_pae_retry_closure':raise ValueError('Passed retry closure required')
 base_d=base['dispositions'];base_verified=sum(v for k,v in base_d.items() if k.endswith('_retrieval_verified'));base_failed=sum(v for k,v in base_d.items() if k.endswith('_retrieval_failed'))
 if base_failed!=closure['retained_base_failures'] or base_failed!=retry['expected_models']:raise ValueError('Base-failure census differs across independent stages')
 retry_d=closure['retry_dispositions']
 if retry_d['verified']!=retry['counts']['verified'] or retry_d['failed']!=retry['counts']['failed'] or retry_d['failed']!=closure['final_explicit_failures']:raise ValueError('Retry census differs across independent stages')
 expected=closure['expected_models']; final=dict(retrieval_verified=base_verified+retry_d['verified'],retrieval_failed=retry_d['failed'])
 if sum(final.values())!=expected:raise ValueError('Final lineage census differs from source queue')
 inputs=[a.preparation,a.recovery,a.base_readback,a.retry,a.retry_readback]
 result=dict(status='completed_full_atlas_pae_recovery_retry_lineage_for_clustering_gate',checked_utc=datetime.now(timezone.utc).isoformat(),expected_models=expected,final_dispositions=final,base_dispositions=base_d,retry_dispositions=retry_d,final_explicit_failures=retry_d['failed'],source_hashes={str(x):sha(x) for x in inputs},scope='Receipt-only provenance composition after independent base and retry closure checks. It records source-bound PAE availability/dispositions for the interrupted-queue lineage, not confidence calibration, structural clustering, homology, or biological inference.',scientific_eligibility=False,gpu=False,new_predictions=0)
 a.receipt.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
