#!/usr/bin/env python3
"""Check the new readback against completed production records and tampering."""
import copy
import json
from pathlib import Path
from audit_matched_reml_analytic_refinement import audit_payload,payload_digest
from ancestral_chain_attempt import sha,write_json


def main():
    plan_path=Path('metadata/matched_reml_analytic_refinement_plan_20260928.json');plan=json.loads(plan_path.read_text())
    paths=sorted(Path(plan['output']).glob('*.json'));checked={};payload=None;old=None
    for path in paths:
        r=json.loads(path.read_text())
        if 'payload' not in r:continue
        assert r['plan_sha256']==sha(plan_path) and r['payload_sha256']==payload_digest(r['payload'])
        source=r['source'];assert sha(source['source_fit'])==source['source_fit_sha256']
        old=json.loads(Path(source['source_fit']).read_text())['payload'];payload=r['payload']
        if payload['status']=='refinement_error_requires_review':continue
        audit_payload(payload,old);checked[str(path)]=sha(path)
    assert checked and payload is not None
    for mutation in ['missing_candidate','false_gradient_flag','wrong_raw_coefficient']:
        bad=copy.deepcopy(payload)
        if mutation=='missing_candidate':bad['candidates'].pop()
        elif mutation=='false_gradient_flag':bad['checks']['projected_gradient_pass']=not bad['checks']['projected_gradient_pass']
        else:bad['raw_unit_beta'][0]+=1
        try:audit_payload(bad,old)
        except AssertionError:pass
        else:raise AssertionError('Accepted '+mutation)
    out=Path('metadata/matched_refinement_readback_checks_20260928.json');assert not out.exists()
    write_json(out,dict(status='completed_record_readback_and_three_tamper_checks_passed',records=len(checked),source_output_hashes=checked,
        pins={str(p):sha(p) for p in [plan_path,Path(__file__),Path('scripts/audit_matched_reml_analytic_refinement.py')]},
        scope='Frozen currently completed subset plus missing-candidate, inverted-gradient-flag and changed-coefficient rejection. Full658 audit remains queued; no whole-grid claim.'))
    print('Verified',len(checked),'completed source-bound records and rejected three altered payloads.')


if __name__=='__main__':main()
