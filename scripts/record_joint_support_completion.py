"""Bind terminal joint-support stages and all declared artifacts into one record."""
import json
from pathlib import Path
import subprocess
from screen_duplication_domain_alignment_coverage import sha

specs=[('fungal-full-joint-covariate-support-20260927.service','results/model_validation/full-joint-covariate-support-20260927-v1','complete_full_joint_covariate_support_pending_readback'),('fungal-full-joint-support-readback-20260927.service','results/model_validation/full-joint-support-readback-20260927-v1','passed_full_joint_covariate_support_certificate_readback'),('fungal-joint-support-failure-summary-20260927.service','results/model_validation/joint-support-failure-summary-20260927-v1','complete_full_joint_support_failure_summary'),('fungal-joint-support-projection-20260927.service','results/model_validation/joint-support-projection-20260927-v1','complete_unresolved_joint_support_projection')]
receipts=[];stages=[]
for unit,folder,status in specs:
    raw=subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True)
    state=dict(l.split('=',1) for l in raw.splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),(unit,state)
    root=Path(folder);rp=root/'receipt.json';r=json.loads(rp.read_text());assert r['status']==status
    for name,h in r['artifacts'].items():assert sha(root/name)==h,name
    receipts.append(r);stages.append(dict(unit=unit,terminal_state=state,receipt=str(rp),receipt_sha256=sha(rp),artifacts_checked=len(r['artifacts'])))
p,b,d,c=receipts
assert b['source_support_receipt_sha256']==stages[0]['receipt_sha256']
assert d['source_receipt_sha256']==c['source_readback_receipt_sha256']==stages[1]['receipt_sha256']
assert b['unique_inputs']==d['unique_inputs']==c['original_inputs']==28808
assert b['full_settings']==c['full_settings']==82944
assert p['classification_counts']==b['classification_counts']==d['classification_counts']
assert d['unresolved_inputs']==c['unresolved_inputs']==sum(n for k,n in b['classification_counts'].items() if k.startswith('unresolved'))
assert sum(c['classification_counts'].values())==c['unresolved_inputs']
result=dict(status='complete_verified_joint_support_and_constructive_followup',stages=stages,original_classification_counts=b['classification_counts'],original_setting_counts=b['setting_classification_counts'],failure_reason_counts=d['failure_reason_counts'],projection_classification_counts=c['classification_counts'],projection_setting_counts=c['setting_classification_counts'],script_sha256=sha(__file__),scope='Full original certificate readback, all-settings mapping, failure diagnosis and original-matrix constructive follow-up completed. Original classifications retained. Numerical hull support does not establish interior overlap, dense sampling, valid uncertainty, adequate covariance or a biological effect.')
out=Path('metadata/joint_support_completed_20260927.json')
with out.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
