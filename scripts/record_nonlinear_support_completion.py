#!/usr/bin/env python3
"""Bind complete nonlinear support checks and separately preserved constructive certificates."""
import json,subprocess
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def main():
    specs=[('fungal-nonlinear-joint-support-20260927.service','nonlinear-joint-support-20260927-v1','complete_full_nonlinear_joint_support_pending_readback'),('fungal-nonlinear-joint-support-readback-20260927.service','nonlinear-joint-support-readback-20260927-v1','passed_full_nonlinear_joint_support_certificate_readback'),('fungal-nonlinear-joint-support-failure-summary-20260927.service','nonlinear-joint-support-failure-summary-20260927-v1','complete_full_nonlinear_joint_support_failure_summary'),('fungal-nonlinear-joint-support-projection-20260927.service','nonlinear-joint-support-projection-20260927-v1','complete_nonlinear_unresolved_joint_support_projection')]
    receipts=[];stages=[]
    for unit,folder,status in specs:
        state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),(unit,state)
        root=Path('results/model_validation')/folder;rp=root/'receipt.json';r=json.loads(rp.read_text());assert r['status']==status
        for name,h in r['artifacts'].items():assert sha(root/name)==h,name
        receipts.append(r);stages.append(dict(unit=unit,terminal_state=state,receipt=str(rp),receipt_sha256=sha(rp),artifacts_checked=len(r['artifacts'])))
    p,b,d,c=receipts
    assert b['source_support_receipt_sha256']==stages[0]['receipt_sha256']
    assert d['source_receipt_sha256']==c['source_support_receipt_sha256']==stages[1]['receipt_sha256']
    assert p['unique_inputs']==b['unique_inputs']==d['unique_inputs']==c['original_inputs']==57616
    assert b['full_settings']==c['full_settings']==165888
    assert p['classification_counts']==b['classification_counts']==d['classification_counts']
    assert b['setting_classification_counts']==c['setting_classification_counts']
    assert d['unresolved_inputs']==c['unresolved_inputs']==sum(v for k,v in b['classification_counts'].items() if k.startswith('unresolved'))
    assert sum(c['classification_counts'].values())==c['unresolved_inputs']
    failures={x['fit_input_id']:x for x in map(json.loads,(Path(stages[2]['receipt']).parent/'unresolved_inputs.jsonl').read_text().splitlines())}
    projected={x['fit_input_id']:x for x in map(json.loads,(Path(stages[3]['receipt']).parent/'projected_certificates.jsonl').read_text().splitlines())}
    assert set(failures)==set(projected) and len(failures)==d['unresolved_inputs']
    for key in failures:
        assert failures[key]['classification']==projected[key]['original_classification']
        assert failures[key]['polynomial_degree']==projected[key]['polynomial_degree']
    supported=[x for x in projected.values() if x['classification']=='supported_by_projected_nonnegative_weights']
    result=dict(status='complete_verified_nonlinear_joint_support_and_constructive_followup',stages=stages,unique_inputs=57616,full_settings=165888,original_classification_counts=b['classification_counts'],original_setting_counts=b['setting_classification_counts'],degree_classification_counts=d['degree_classification_counts'],failure_reason_counts=d['failure_reason_counts'],maximum_original_negative_weight_mass=d['maximum_negative_mass'],projection_classification_counts=c['classification_counts'],projection_setting_counts=c['projected_setting_classification_counts'],maximum_supported_projected_distance=max([x['primal_distance'] for x in supported],default=None),script_sha256=sha(__file__),scope='Full nonlinear certificate readback, complete setting mapping, failure diagnosis and original-matrix constructive follow-up verified. Original solver classifications retained. Nonnegative projected certificates establish numerical hull support only; no interior-overlap, dense-sampling, covariance-adequacy, calibrated-uncertainty or biological-effect claim.')
    with Path('metadata/nonlinear_joint_support_completed_20260927.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
