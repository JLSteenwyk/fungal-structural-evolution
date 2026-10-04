#!/usr/bin/env python3
"""Bind original diagnostic tool wait, wrapper payload and invocation journals."""
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
        assert v['status'] in ['completed_exact_failed_cohort_covariance_diagnostic_v1',
                               'completed_extended_precision_failed_cohort_raw_reference_v1',
                               'completed_scoped_failed_cohort_species_raw_reference_v1']
        assert v['biological_fits'] == 0
        assert v['production_tolerance_changed'] is v['scientific_eligibility'] is v['gpu'] is False
        if v['status'] == 'completed_scoped_failed_cohort_species_raw_reference_v1':
            assert v['records'] == 22881 and v['trees_checked'] == 5 and v['original_audits_checked'] == 300
            assert v['raw_entries_per_audit'] == 2 and v['original_guard_reexecuted'] is False
            for n,h in v['source_hashes'].items(): bind(pins,n,h)
            for n,h in v['artifacts'].items(): bind(pins,n,h)
            counts = dict(trees_checked=5,original_audits_checked=300)
        elif 'artifacts' in v:
            assert v['selected_groups'] == 40
            for n,h in v['artifacts'].items(): bind(pins,n,h)
            assert sha(v['full_source_hash_archive']) == v['full_source_hash_archive_sha256']
            counts = v['guard_status_counts']
        else:
            assert v['selected_groups'] == 40
            for n,h in v['source_hashes'].items(): bind(pins,n,h)
            counts = dict(raw_reference_contexts=v['raw_reference_contexts'])
        bind(pins,a.validation)
        status = 'verified_original_retained_covariance_diagnostic_wait_zero'
    else:
        assert a.actual_tool_exit == 1 and e['status'] == 'failed_stage_retained'
        assert e['receipt_sha256'] is None and not a.validation.exists()
        status = 'retained_original_failed_retained_covariance_diagnostic_attempt'
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
    journal = a.execution.with_suffix('')/'original-invocation-journal.jsonl'
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
        scope='Exact original wait, wrapper PID/create/command/invocation, complete initial/terminal '
              'payloads and manager start/end, all diagnostic artifacts freshly hashed. Producer '
              'full-source verification is inherited through its hashed archive, not repeated here. '
              'Reused generic wrapper scope names an old checker; actual command identifies this '
              'diagnostic. A successful diagnostic exit can retain failed original guards; it '
              'does not accept their sources or fit biological effects. Original jobs unchanged.')
    with a.output.open('x') as f: json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper']},indent=2))


if __name__ == '__main__': main()
