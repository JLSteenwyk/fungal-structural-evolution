#!/usr/bin/env python3
"""Preserve original timing failure, partial cohorts and unchanged waiting handles."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess
import psutil

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify
from record_project_runtime_checkpoint_v4 import fingerprint,live_record


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    pp=Path('metadata/full_retained_shared_entity_timing_plan_20261003_v1.json')
    plan=json.loads(pp.read_text());verify(plan['pins']);pins=dict(plan['pins']);bind(pins,pp)
    inventory=Path('metadata/full_retained_shared_entity_timing_launches_20261003_v1.json')
    launches=json.loads(inventory.read_text())['launches'];bind(pins,inventory)
    original=json.loads(Path(launches[0]).read_text());bind(pins,launches[0])
    assert sha(original['plan'])==original['plan_sha256'] and fingerprint(original) is None
    rows=[json.loads(l) for l in subprocess.check_output(['journalctl','--user','-u',original['unit'],'-o','json','--no-pager'],text=True).splitlines()]
    exact=[r for r in rows if r.get('_PID')==str(original['pid']) and r.get('_CMDLINE')==' '.join(original['cmdline'])]
    invs={r['_SYSTEMD_INVOCATION_ID'] for r in exact};assert len(invs)==1
    inv=next(iter(invs));rows=[r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    messages='\n'.join(r.get('MESSAGE','') for r in rows)
    assert 'Max absolute difference among violations: 2.56917974' in messages
    assert 'rtol=3e-09, atol=2e-08' in messages and 'backend_guard' in messages
    assert 'returned non-zero exit status 1' in messages
    ends=[r for r in rows if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')]
    assert len(ends)==1 and any('status=1/FAILURE' in r.get('MESSAGE','') for r in rows)
    journal=a.output.with_suffix('');journal.mkdir(exist_ok=False)
    jp=journal/'original-producer-invocation.jsonl'
    with jp.open('x') as f:
        for row in rows:f.write(json.dumps(row,sort_keys=True)+'\n')
    bind(pins,jp)
    waiting=[]
    for lp in launches[1:]:
        bind(pins,lp);r=json.loads(Path(lp).read_text());assert sha(r['plan'])==r['plan_sha256']
        proc=fingerprint(r);assert proc is not None
        waiting.append(live_record(r,proc))
    root=Path(plan['output']);assert root.exists() and not (root/'receipt.json').exists()
    assert not Path(plan['completion']).exists()
    partial=list(sorted(root.rglob('*')))
    for q in partial:
        if q.is_file():bind(pins,q)
    verify(pins)
    result=dict(status='retained_original_full_retained_timing_numerical_guard_failure',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_producer_launch=original,
        invocation_id=inv,original_producer_terminal_exit_code=1,
        manager_completion_records=1,exact_original_wrapper_journal_rows=len(exact),
        current_dependency_waiters=waiting,completed_partial_cohort_receipts=len(list((root/'cohorts').glob('*.receipt.json'))),
        failed_guard='retained_shared_entity_candidate.backend_guard',unchanged_guard_rtol=3e-9,unchanged_guard_atol=2e-8,
        maximum_reported_absolute_difference=2.56917974,maximum_reported_relative_difference=1.52149497e-8,
        reported_mismatched_elements=9,reported_matrix_elements=25,source_hashes=pins,
        original_timing_completion_present=False,production_fit_launched=False,source_or_output_changed=False,
        original_jobs_restarted=False,all_eight_aims_incomplete=True,gpu=False,
        scope='Exact original producer PID/create/command launch matched to invocation-specific native failure '
              'traceback and manager exit1; unchanged source pins and every existing partial artifact hashed. '
              'Both exact original readback/closure controllers are live waiting for a receipt never produced. '
              'No relaxed tolerance, retry, replacement, numerical source acceptance or biological fit. '
              'Separate four-control weighted numerical stage continues and needs its own full closure.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','current_dependency_waiters','original_producer_launch']},indent=2))


if __name__=='__main__':main()
