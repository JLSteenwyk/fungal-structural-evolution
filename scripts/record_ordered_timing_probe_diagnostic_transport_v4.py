#!/usr/bin/env python3
"""Bind the original ordered-probe diagnostic wait, payload and source/artifact verification."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execution', type=Path, required=True)
    p.add_argument('--validation', type=Path, required=True)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--unit', required=True)
    p.add_argument('--actual-tool-session', type=int, required=True)
    p.add_argument('--actual-tool-exit', type=int, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    e = json.loads(a.execution.read_text()); assert not e['timed_out']
    assert e['exit_code'] == a.actual_tool_exit
    pins = {}; validation_status = None; counts = None
    for mapping in [e['source_hashes'], e['artifacts']]:
        for n,h in mapping.items(): bind(pins,n,h)
    if a.actual_tool_exit == 0:
        assert e['status'] == 'exited_zero_with_receipt' and sha(a.validation) == e['receipt_sha256']
        v = json.loads(a.validation.read_text()); validation_status = v['status']
        assert v['status']=='completed_ordered_cohort_timing_probe_diagnostic_v4'
        assert v['biological_fits']==0 and v['planned_groups']==40
        assert 0<v['executed_groups']<=40 and v['executed_groups']+len(v['unattempted_group_ids'])==40
        assert sum(v['probe_status_counts'].values())==v['executed_groups']
        assert v['production_tolerance_changed'] is v['scientific_eligibility'] is v['gpu'] is False
        assert v['original_jobs_restarted'] is False
        assert v['full_original_source_verification_inherited'] is True
        assert v['full_original_source_archive_rehashed'] is False
        assert v['original_prior_ten_cohort_probe_sequence_replayed'] is False
        for n,h in v['source_hashes'].items():bind(pins,n,h)
        for n,h in v['artifacts'].items():bind(pins,n,h)
        counts=dict(executed_groups=v['executed_groups'],probe_status_counts=v['probe_status_counts'],
                    input_mutation=v['detected_input_mutation'],unattempted_group_ids=v['unattempted_group_ids'])
        bind(pins,a.validation)
        status = 'verified_original_ordered_timing_probe_diagnostic_wait_zero'
    else:
        assert a.actual_tool_exit == 1 and e['status'] == 'failed_stage_retained'
        assert e['receipt_sha256'] is None and not a.validation.exists()
        status = 'retained_original_failed_ordered_timing_probe_diagnostic_attempt'
    for path in a.root.rglob('*'):
        if path.is_file(): bind(pins,path)
    rows = [json.loads(l) for l in subprocess.check_output(
        ['journalctl','--user','-u',a.unit,'-o','json','--no-pager'], text=True).splitlines()]
    inv = e['invocation_id']
    rows = [r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    exact = [r for r in rows if r.get('_PID') == str(e['wrapper']['pid']) and r.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])]
    assert len(exact) == 2 and all(r['_SYSTEMD_INVOCATION_ID'] == inv for r in exact)
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=e['wrapper'],invocation_id=inv)
    terminal = {k:v for k,v in e.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and 'Started ' in r.get('MESSAGE','')]
    ends = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    journal = a.execution.with_suffix('')/'original-invocation-journal-line-trace-v4.jsonl'
    with journal.open('x') as f:
        for row in rows: f.write(json.dumps(row,sort_keys=True)+'\n')
    for path in [a.execution,journal,Path(__file__)]: bind(pins,path)
    verify(pins)
    result = dict(status=status,checked_utc=datetime.now(timezone.utc).isoformat(),
        validation_status=validation_status,diagnostic_counts=counts,
        validation_sha256=e['receipt_sha256'],actual_tool_session_id=a.actual_tool_session,
        actual_tool_terminal_exit_code=a.actual_tool_exit,invocation_id=inv,unit=a.unit,
        wrapper=e['wrapper'],entire_terminal_payload_matched=True,
        original_start_records=1,original_completion_records=1,source_hashes=pins,
        full_producer_source_archive_rehashed=False,scientific_eligibility=False,gpu=False,
        source_verification_scope='consumed_export_bindings_only_with_inherited_full_source_verification',
        scope='Exact original ordered-probe diagnostic wait, wrapper PID/create/command/invocation, '
              'complete initial/terminal payload and manager start/end verified. Every consumed '
              'source/code/export binding and produced diagnostic artifact freshly rehashed. Full '
              'original2.38million-source replay is inherited from the separate closed producer. '
              'Probe reviews/rejections/mutations/unattempted groups and any captured factor snapshots/bitwise differences remain explicit. Not a '
              'prior10cohort execution-state replay, covariance repair, tolerance relaxation, '
              'original restart or biological acceptance. Generic wrapper kernel scope is '
              'superseded by its actual ordered-probe command.')
    with a.output.open('x') as f: json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper']},indent=2))


if __name__ == '__main__': main()
