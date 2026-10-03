#!/usr/bin/env python3
"""Close full recovery accounting using both exact original completion journals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    args=parser.parse_args(); plan=json.loads(args.plan.read_text()); bindings=dict(plan['pins']); bind(bindings,args.plan)
    audit_plan=json.loads(Path(plan['audit_plan']).read_text()); rp=Path(audit_plan['output'])/'receipt.json'
    audit=json.loads(rp.read_text()); bind(bindings,rp)
    assert audit['status']=='passed_full_baliphy_recovery_integrity_and_original_grid_readback'
    assert audit['plan_sha256']==sha(plan['audit_plan'])
    retry=json.loads(Path(audit['recovery_receipt']).read_text())
    assert sha(audit['recovery_receipt'])==audit['recovery_receipt_sha256']
    assert audit['recovery_plan_sha256']==retry['plan_sha256']==sha(audit_plan['recovery_plan'])
    overlay=json.loads(Path(audit['overlay']).read_text()); bind(bindings,audit['overlay'],audit['overlay_sha256'])
    assert overlay['status']=='verified_full_original_grid_with_whole_recovery_attempt_overlay'
    assert overlay['original_samples_concatenated'] is False and overlay['scientific_eligibility'] is False
    assert all(audit[k]==v for k,v in overlay['summary'].items())
    assert len(overlay['rows'])==overlay['summary']['full_native_chains']==1620
    assert len(overlay['quartets'])==overlay['summary']['full_quartets']==405
    for p,d in audit['source_hashes'].items(): bind(bindings,p,d)
    verify(bindings)
    root=Path(plan['output']); root.mkdir(exist_ok=False); archive=root/'completion_archive.json'; inner=root/'closure_plan.json'
    spec=dict(output=str(archive), completed_status='complete_verified_baliphy_memory_recovery_accounting_archive',
        evidence=dict(recovery=dict(path=audit['recovery_receipt'],expected=dict(status='complete_baliphy_memory_recovery_pending_full_readback',
            plan_sha256=audit['recovery_plan_sha256'])),
            audit=dict(path=str(rp),expected=dict(status=audit['status'],plan_sha256=audit['plan_sha256'],**overlay['summary']))),
        links=[dict(**{'from':'audit'},field='recovery_receipt_sha256',to=audit['recovery_receipt'])],
        launches=plan['launches'],pins=bindings,summary=overlay['summary'],scope=plan['scope'])
    with inner.open('x') as f:json.dump(spec,f,indent=2);f.write('\n')
    subprocess.run([sys.executable,'scripts/record_completed_process_handoffs_v2.py','--plan',str(inner)],check=True)
    proof=json.loads(archive.read_text());assert len(proof['services'])==2
    completed=dict(status='complete_verified_full_baliphy_memory_recovery_accounting',**overlay['summary'],
        source_plan=plan['audit_plan'],source_plan_sha256=sha(plan['audit_plan']),
        overlay=audit['overlay'],overlay_sha256=audit['overlay_sha256'],
        full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),
        bound_source_hashes=len(proof['source_hashes']),exact_process_journals_checked=2,
        scientific_eligibility=False,scope=plan['scope'])
    with Path(plan['completion']).open('x') as f:json.dump(completed,f,indent=2);f.write('\n')
    print(json.dumps(completed,indent=2),flush=True)


if __name__=='__main__':main()
